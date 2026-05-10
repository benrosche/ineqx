# ============================================================================ #
# Bootstrap standard errors for causal variance decomposition
# ============================================================================ #

# ---------------------------------------------------------------------------- #
# Boot config object
# ---------------------------------------------------------------------------- #

#' Create a bootstrap configuration
#'
#' Groups bootstrap settings into a single object that can be passed to
#' \code{\link{ineqx}} via the \code{se} argument, or to plot methods via
#' the \code{ci} argument. All model arguments (data, formulas, etc.) are
#' inherited from the \code{ineqx()} call automatically.
#'
#' @param B Integer, number of bootstrap replicates. Default 100.
#' @param parallel Logical, use parallel computation. Default FALSE.
#' @param ncores Integer, number of cores for parallel. Default: all but one.
#' @param seed Integer, random seed for reproducibility. Default NULL.
#' @param verbose Logical, print progress messages. Default TRUE.
#'
#' @return An object of class \code{"ineqx_boot_config"}
#'
#' @examples
#' \dontrun{
#' # Bootstrap SEs with 500 replicates
#' ineqx("income", treat = "treat", group = "group", data = mydata,
#'       formula_mu = ~ treat * group,
#'       formula_sigma = ~ treat * group,
#'       se = boot_config(B = 500, seed = 42))
#'
#' # Bootstrap CIs in plot (explicit config or "boot" shorthand)
#' plot(desc_result, ci = boot_config(B = 100))
#' plot(desc_result, ci = "boot")  # equivalent to boot_config()
#' }
#'
#' @export
boot_config <- function(B = 100L, parallel = FALSE, ncores = NULL,
                         seed = NULL, verbose = TRUE) {

  if (B < 2) stop("'B' must be at least 2")

  structure(
    list(
      B = as.integer(B),
      parallel = parallel,
      ncores = ncores,
      seed = seed,
      verbose = verbose
    ),
    class = "ineqx_boot_config"
  )
}


# ---------------------------------------------------------------------------- #
# Main bootstrap function
# ---------------------------------------------------------------------------- #

#' Bootstrap standard errors for causal variance decomposition
#'
#' Computes standard errors by nonparametric bootstrap: resample individuals,
#' re-estimate the GAMLSS, re-extract parameters, and re-compute the
#' decomposition B times. Standard errors are the standard deviation across
#' replicates. Percentile confidence intervals are also computed.
#'
#' For longitudinal data with repeated cross-sections, resampling is performed
#' within each time period independently, preserving the sample size per period.
#'
#' @param data Data.frame, the original individual-level data
#' @param formula_mu Two-sided formula for the mean (mu) equation
#' @param formula_sigma One-sided formula for the log-SD (sigma) equation
#' @param treat Character, treatment variable name
#' @param group Character, grouping variable name
#' @param time Character, time variable name. NULL for cross-sectional.
#' @param post Character, pre/post indicator for DiD. NULL for simple diff.
#' @param ref Numeric, reference time period for longitudinal decomposition
#' @param ystat Character, \code{"Var"} or \code{"CV2"}
#' @param order Character vector of length 3, decomposition ordering
#' @param B Integer, number of bootstrap replicates. Default 100.
#' @param parallel Logical, use parallel computation. Default FALSE.
#' @param ncores Integer, number of cores. Default: all but one.
#' @param seed Integer, random seed. Default NULL.
#' @param verbose Logical, print progress. Default TRUE.
#'
#' @return An object of class \code{"ineqx_boot"} containing:
#' \describe{
#'   \item{se}{List of SEs matching \code{\link{delta_method_se}} output structure}
#'   \item{replicates}{Matrix (B_successful x K) of replicate estimates}
#'   \item{ci}{List of 95\% percentile confidence intervals}
#'   \item{B}{Total requested replicates}
#'   \item{B_successful}{Number of successful replicates}
#'   \item{B_failed}{Number of failed replicates}
#'   \item{type}{\code{"cross"} or \code{"longit"}}
#'   \item{point_estimates}{Named vector from original (non-bootstrap) decomposition}
#'   \item{seed}{The random seed used}
#' }
#'
#' @keywords internal
bootstrap_se <- function(data, formula_mu, formula_sigma,
                          treat, group,
                          time = NULL, post = NULL,
                          ref = NULL, ystat = "Var",
                          order = c("behavioral", "compositional",
                                    "pretreatment"),
                          B = 100L,
                          parallel = FALSE, ncores = NULL,
                          seed = NULL, verbose = TRUE,
                          blend_params = NULL) {

  if (!requireNamespace("gamlss", quietly = TRUE)) {
    stop("Package 'gamlss' is required for bootstrap_se(). ",
         "Install it with install.packages('gamlss')")
  }

  ystat <- match.arg(ystat, c("Var", "CV2"))
  B <- as.integer(B)
  if (B < 2) stop("'B' must be at least 2")

  # Slim `data` to only the columns the bootstrap actually needs. On Windows
  # PSOCK clusters every worker receives a full copy, so dropping unused
  # columns up-front cuts memory by a constant factor of ncores.
  needed <- unique(c(
    all.vars(formula_mu), all.vars(formula_sigma),
    treat, group, time, post
  ))
  needed <- intersect(needed, names(data))
  if (length(needed)) data <- data[, needed, drop = FALSE]

  # Determine type
  is_longit <- !is.null(time)
  type <- if (is_longit) "longit" else "cross"

  # Set time_var for internal use
  if (is.null(time)) {
    data$.time <- 1L
    time_var <- ".time"
  } else {
    time_var <- time
  }

  # Generate bootstrap indices
  boot_indices <- .generate_boot_indices(data, time_var, B, seed)

  # Run original decomposition for point estimates
  orig_result <- .boot_one_replicate(
    data = data,
    formula_mu = formula_mu, formula_sigma = formula_sigma,
    treat = treat, group = group, time = time, post = post,
    ref = ref, ystat = ystat, order = order,
    resample_ids = seq_len(nrow(data)),
    blend_params = blend_params
  )
  if (is.null(orig_result)) {
    stop("Original model estimation failed. Cannot compute bootstrap SEs.")
  }

  # Run bootstrap replicates
  if (parallel) {
    results <- .boot_parallel(
      data = data, formula_mu = formula_mu, formula_sigma = formula_sigma,
      treat = treat, group = group, time = time, post = post,
      ref = ref, ystat = ystat, order = order,
      boot_indices = boot_indices, ncores = ncores,
      blend_params = blend_params
    )
  } else {
    results <- vector("list", B)
    bar_width <- 30L
    report_at <- unique(c(round(seq(0, B, length.out = 11L)), B))
    for (b in seq_len(B)) {
      results[[b]] <- tryCatch(
        .boot_one_replicate(
          data = data,
          formula_mu = formula_mu, formula_sigma = formula_sigma,
          treat = treat, group = group, time = time, post = post,
          ref = ref, ystat = ystat, order = order,
          resample_ids = boot_indices[[b]],
          blend_params = blend_params
        ),
        error = function(e) NULL
      )
      if (verbose && b %in% report_at) {
        filled <- round(bar_width * b / B)
        bar <- paste0(strrep("=", filled), strrep(" ", bar_width - filled))
        message(sprintf("  Bootstrap SE [%s] %d/%d", bar, b, B))
      }
    }
  }

  # Filter out failed replicates
  successful <- !vapply(results, is.null, logical(1))
  B_successful <- sum(successful)
  B_failed <- B - B_successful

  if (B_failed > 0 && B_failed / B > 0.2) {
    warning(B_failed, " of ", B, " bootstrap replicates failed (",
            round(100 * B_failed / B), "%). ",
            "Consider checking model specification or increasing sample size.",
            call. = FALSE)
  }
  min_required <- min(B, 20)
  if (B_successful < min_required) {
    stop("Only ", B_successful, " of ", B, " bootstrap replicates succeeded. ",
         "Cannot compute reliable SEs. Check model specification.")
  }

  # Stack results into matrix
  good_results <- results[successful]
  replicates <- do.call(rbind, good_results)
  colnames(replicates) <- names(good_results[[1]])

  # Compute SEs and CIs
  boot_sds <- apply(replicates, 2, stats::sd)
  boot_ci_lower <- apply(replicates, 2, stats::quantile, probs = 0.025)
  boot_ci_upper <- apply(replicates, 2, stats::quantile, probs = 0.975)

  # Format SEs to match delta_method_se output structure
  se_formatted <- .format_boot_se(boot_sds, type, orig_result)

  # Format CIs
  ci_formatted <- lapply(names(boot_sds), function(nm) {
    c(lower = unname(boot_ci_lower[nm]), upper = unname(boot_ci_upper[nm]))
  })
  names(ci_formatted) <- names(boot_sds)

  structure(
    list(
      se = se_formatted,
      replicates = replicates,
      ci = ci_formatted,
      B = B,
      B_successful = B_successful,
      B_failed = B_failed,
      type = type,
      point_estimates = orig_result,
      seed = seed
    ),
    class = "ineqx_boot"
  )
}


# ---------------------------------------------------------------------------- #
# Attach bootstrap SEs to a decomposition result
# ---------------------------------------------------------------------------- #

#' Attach bootstrap SEs to a causal decomposition result
#'
#' Replaces the standard errors in a decomposition result object with
#' bootstrap SEs. This is useful when you've computed the decomposition
#' and bootstrap SEs separately.
#'
#' @param result An \code{ineqx_causal_cross} or \code{ineqx_causal_longit} object
#' @param boot An \code{ineqx_boot} object
#'
#' @return The \code{result} object with \code{$se} replaced by bootstrap SEs
#'   and additional fields \code{$se_method} and \code{$boot}.
#'
#' @keywords internal
attach_boot_se <- function(result, boot) {
  stopifnot(inherits(boot, "ineqx_boot"))
  stopifnot(inherits(result, "ineqx_causal_cross") ||
            inherits(result, "ineqx_causal_longit"))

  result$se <- boot$se
  result$se_method <- "bootstrap"
  result$boot <- boot
  result
}


# ---------------------------------------------------------------------------- #
# Print method for ineqx_boot
# ---------------------------------------------------------------------------- #

#' @export
print.ineqx_boot <- function(x, ...) {
  cat("Bootstrap standard errors\n")
  cat("  Replicates:", x$B_successful, "of", x$B,
      if (x$B_failed > 0) paste0("(", x$B_failed, " failed)") else "",
      "\n")
  cat("  Type:", x$type, "\n")
  if (!is.null(x$seed)) cat("  Seed:", x$seed, "\n")
  cat("\nStandard errors:\n")

  if (x$type == "cross") {
    se <- x$se
    cat("  SE(tau_B):    ", round(se$se_tau_B, 4), "\n")
    cat("  SE(tau_W):    ", round(se$se_tau_W, 4), "\n")
    cat("  SE(tau_total):", round(se$se_tau_total, 4), "\n")
  } else {
    for (t_name in names(x$se)) {
      cat("\n  Time:", t_name, "\n")
      se_t <- x$se[[t_name]]
      for (nm in names(se_t)) {
        cat("    ", nm, ":", round(se_t[[nm]], 4), "\n")
      }
    }
  }

  invisible(x)
}


# ---------------------------------------------------------------------------- #
# Internal: Single bootstrap replicate
# ---------------------------------------------------------------------------- #

.boot_one_replicate <- function(data, formula_mu, formula_sigma,
                                 treat, group, time, post,
                                 ref, ystat, order, resample_ids,
                                 blend_params = NULL) {

  # Resample
  boot_data <- data[resample_ids, , drop = FALSE]

  # Fit GAMLSS
  model <- gamlss::gamlss(
    formula = formula_mu,
    sigma.formula = formula_sigma,
    data = boot_data,
    trace = FALSE
  )

  # Extract parameters (no vcov needed)
  params <- ineqx_params(
    model = model, data = boot_data,
    treat = treat, group = group,
    time = time, post = post,
    ystat = ystat, vcov = FALSE, verbose = FALSE
  )

  # Free the fitted gamlss and resampled data — decomposition only needs params
  rm(model, boot_data)
  gc(verbose = FALSE)

  # Blend with user params if provided
  if (!is.null(blend_params)) {
    params <- .blend_causal_params(blend_params, params, ref)
  }

  # Helper: append pre-period anchor estimates if available (DiD models)
  .append_did_param <- function(acc, pdata_row, suffix) {
    extra_cols <- c("mu1", "sigma1", "mu0_pre", "mu1_pre",
                    "sigma0_pre", "sigma1_pre")
    for (col in extra_cols) {
      if (col %in% names(pdata_row) && !is.na(pdata_row[[col]])) {
        acc <- c(acc, stats::setNames(pdata_row[[col]],
                                       paste0(col, "_", suffix)))
      }
    }
    acc
  }

  # Compute decomposition
  if (params$type == "cross_sectional") {
    result <- causal_decompose_cross(params, ref = ref)
    # Also store beta/lambda/mu0/sigma0 per group for *.params plot CIs
    pdata <- params$data
    param_ests <- c()
    for (g in sort(unique(pdata$group))) {
      idx <- pdata$group == g
      param_ests <- c(param_ests,
        stats::setNames(pdata$beta[idx],   paste0("beta_",   g)),
        stats::setNames(pdata$lambda[idx], paste0("lambda_", g)),
        stats::setNames(pdata$mu0[idx],    paste0("mu0_",    g)),
        stats::setNames(pdata$sigma0[idx], paste0("sigma0_", g))
      )
      param_ests <- .append_did_param(param_ests, pdata[idx, ], g)
    }
    c(tau_B = result$tau_B,
      tau_W = result$tau_W,
      tau_total = result$tau_total,
      param_ests)
  } else {
    result <- causal_decompose_longit(params, order = order, ref = ref)
    # Also store beta/lambda/mu0/sigma0 per group × time for *.params plot CIs
    pdata <- params$data
    param_ests <- c()
    for (t in sort(unique(pdata$time))) {
      for (g in sort(unique(pdata$group))) {
        idx <- pdata$time == t & pdata$group == g
        if (!any(idx)) next
        param_ests <- c(param_ests,
          stats::setNames(pdata$beta[idx],   paste0("beta_",   g, "_", t)),
          stats::setNames(pdata$lambda[idx], paste0("lambda_", g, "_", t)),
          stats::setNames(pdata$mu0[idx],    paste0("mu0_",    g, "_", t)),
          stats::setNames(pdata$sigma0[idx], paste0("sigma0_", g, "_", t))
        )
        param_ests <- .append_did_param(param_ests, pdata[idx, ],
                                         paste0(g, "_", t))
      }
    }
    # Flatten decomposition into named vector
    estimates <- c()
    for (t_name in names(result$results)) {
      r <- result$results[[t_name]]
      t_ests <- c(
        # Split components (10)
        Delta_beta_B   = r$Delta_beta_B,
        Delta_beta_W   = r$Delta_beta_W,
        Delta_lambda_B = r$Delta_lambda_B,
        Delta_lambda_W = r$Delta_lambda_W,
        Delta_pi_B     = r$Delta_pi_B,
        Delta_pi_W     = r$Delta_pi_W,
        Delta_mu_B     = r$Delta_mu_B,
        Delta_mu_W     = r$Delta_mu_W,
        Delta_sigma_B  = r$Delta_sigma_B,
        Delta_sigma_W  = r$Delta_sigma_W,
        # Aggregate parameter components
        Delta_beta   = r$Delta_beta,
        Delta_lambda = r$Delta_lambda,
        Delta_mu     = r$Delta_mu,
        Delta_sigma  = r$Delta_sigma,
        # Combined
        Delta_behavioral    = r$Delta_behavioral,
        Delta_compositional = r$Delta_compositional,
        Delta_pretreatment  = r$Delta_pretreatment,
        # Totals
        Delta_B      = r$Delta_B,
        Delta_W      = r$Delta_W,
        Delta_total  = r$Delta_total
      )
      names(t_ests) <- paste0(names(t_ests), "_", t_name)
      estimates <- c(estimates, t_ests)
    }

    # Cross-sectional component levels per time (for wibe plot CIs):
    # tau_{B,W,total}_<t>, het_{B,W}_<t>, cov_{B,W}_<t>
    # ref time levels live in components_t0 / tau_*_t0 of any non-ref result
    level_ests <- c()
    if (length(result$results) > 0) {
      ref_t <- as.character(result$ref %||% NA)
      first_r <- result$results[[1]]
      if (!is.na(ref_t)) {
        ct0 <- first_r$components_t0
        level_ests <- c(level_ests,
          stats::setNames(first_r$tau_B_t0,            paste0("tau_B_",     ref_t)),
          stats::setNames(first_r$tau_W_t0,            paste0("tau_W_",     ref_t)),
          stats::setNames(first_r$tau_B_t0 + first_r$tau_W_t0,
                          paste0("tau_total_", ref_t)),
          stats::setNames(ct0$het_B,                   paste0("het_B_",     ref_t)),
          stats::setNames(ct0$cov_B,                   paste0("cov_B_",     ref_t)),
          stats::setNames(ct0$rescale_B,               paste0("rescale_B_", ref_t)),
          stats::setNames(ct0$het_W,                   paste0("het_W_",     ref_t)),
          stats::setNames(ct0$cov_W,                   paste0("cov_W_",     ref_t)),
          stats::setNames(ct0$rescale_W,               paste0("rescale_W_", ref_t))
        )
      }
      for (t_name in names(result$results)) {
        r <- result$results[[t_name]]
        ct <- r$components_t
        level_ests <- c(level_ests,
          stats::setNames(r$tau_B_t,            paste0("tau_B_",     t_name)),
          stats::setNames(r$tau_W_t,            paste0("tau_W_",     t_name)),
          stats::setNames(r$tau_B_t + r$tau_W_t,
                          paste0("tau_total_", t_name)),
          stats::setNames(ct$het_B,             paste0("het_B_",     t_name)),
          stats::setNames(ct$cov_B,             paste0("cov_B_",     t_name)),
          stats::setNames(ct$rescale_B,         paste0("rescale_B_", t_name)),
          stats::setNames(ct$het_W,             paste0("het_W_",     t_name)),
          stats::setNames(ct$cov_W,             paste0("cov_W_",     t_name)),
          stats::setNames(ct$rescale_W,         paste0("rescale_W_", t_name))
        )
      }
    }

    c(estimates, level_ests, param_ests)
  }
}


# ---------------------------------------------------------------------------- #
# Internal: Generate bootstrap resampling indices
# ---------------------------------------------------------------------------- #

.generate_boot_indices <- function(data, time_var, B, seed) {
  if (!is.null(seed)) set.seed(seed)

  n <- nrow(data)
  time_vals <- unique(data[[time_var]])

  if (length(time_vals) <= 1) {
    # Cross-sectional: simple resampling
    lapply(seq_len(B), function(b) sample.int(n, n, replace = TRUE))
  } else {
    # Longitudinal: resample within each time period
    period_idx <- split(seq_len(n), data[[time_var]])
    lapply(seq_len(B), function(b) {
      unlist(lapply(period_idx, function(idx) {
        sample(idx, length(idx), replace = TRUE)
      }), use.names = FALSE)
    })
  }
}


# ---------------------------------------------------------------------------- #
# Internal: Parallel bootstrap execution
# ---------------------------------------------------------------------------- #

.boot_parallel <- function(data, formula_mu, formula_sigma,
                            treat, group, time, post,
                            ref, ystat, order,
                            boot_indices, ncores,
                            blend_params = NULL) {

  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("Package 'parallel' is required for parallel bootstrap")
  }
  if (is.null(ncores)) ncores <- max(1L, parallel::detectCores() - 1L)

  cl <- parallel::makeCluster(ncores)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  # Export data and settings to workers
  parallel::clusterExport(cl, c(
    "data", "formula_mu", "formula_sigma", "treat", "group",
    "time", "post", "ref", "ystat", "order",
    "blend_params"
  ), envir = environment())

  parallel::clusterEvalQ(cl, {
    requireNamespace("ineqx", quietly = TRUE)
    requireNamespace("gamlss", quietly = TRUE)
  })

  parallel::parLapply(cl, boot_indices, function(ids) {
    tryCatch(
      ineqx:::.boot_one_replicate(
        data = data, formula_mu = formula_mu,
        formula_sigma = formula_sigma,
        treat = treat, group = group, time = time,
        post = post, ref = ref, ystat = ystat,
        order = order, resample_ids = ids,
        blend_params = blend_params
      ),
      error = function(e) NULL
    )
  })
}


# ---------------------------------------------------------------------------- #
# Internal: Format bootstrap SEs to match delta_method_se output
# ---------------------------------------------------------------------------- #

.format_boot_se <- function(boot_sds, type, orig_result) {
  if (type == "cross") {
    list(
      se_tau_B = unname(boot_sds["tau_B"]),
      se_tau_W = unname(boot_sds["tau_W"]),
      se_tau_total = unname(boot_sds["tau_total"])
    )
  } else {
    # Longitudinal: group by time period
    # Names are like "Delta_beta_2010", "Delta_lambda_2010", etc.
    all_names <- names(boot_sds)

    # Extract unique time suffixes
    time_suffixes <- unique(sub("^.*_([^_]+)$", "\\1", all_names))

    se_list <- list()
    for (t_name in time_suffixes) {
      pattern <- paste0("_", t_name, "$")
      t_names <- grep(pattern, all_names, value = TRUE)
      t_ses <- boot_sds[t_names]

      # Strip the time suffix and add "se_" prefix
      base_names <- sub(pattern, "", names(t_ses))
      se_names <- paste0("se_", base_names)

      se_t <- as.list(unname(t_ses))
      names(se_t) <- se_names
      se_list[[t_name]] <- se_t
    }
    se_list
  }
}
