# ============================================================================ #
# Bootstrap standard errors for causal variance decomposition
# ============================================================================ #

# ---------------------------------------------------------------------------- #
# Boot config object
# ---------------------------------------------------------------------------- #

#' Create a bootstrap configuration
#'
#' Groups all bootstrap settings into a single object that can be passed to
#' \code{\link{ineqx}} via the \code{boot} argument.
#'
#' @param data Data.frame, the original individual-level data used to fit
#'   the GAMLSS model
#' @param formula_mu Formula for the mean (mu) equation
#' @param formula_sigma Formula for the log-SD (sigma) equation (one-sided,
#'   e.g., \code{~treat*group})
#' @param treat Character, name of the treatment variable (coded 0/1)
#' @param group Character, name of the grouping variable. Required when
#'   using \code{boot_config} via \code{\link{ineqx}}.
#' @param time Character, name of the time variable. NULL for cross-sectional.
#' @param post Character, name of the pre/post indicator for DiD designs.
#'   NULL for simple difference estimator.
#' @param B Integer, number of bootstrap replicates. Default 200.
#' @param family A gamlss.family object (e.g., \code{gamlss.dist::NO()}).
#'   If NULL, uses the normal distribution.
#' @param parallel Logical, use parallel computation. Default FALSE.
#' @param ncores Integer, number of cores for parallel. Default: all but one.
#' @param seed Integer, random seed for reproducibility. Default NULL.
#' @param verbose Logical, print progress messages. Default TRUE.
#' @param gamlss_args Named list of additional arguments passed to
#'   \code{gamlss::gamlss()}.
#'
#' @return An object of class \code{"ineqx_boot_config"}
#'
#' @examples
#' \dontrun{
#' bc <- boot_config(
#'   data = mydata,
#'   formula_mu = y ~ treat * group,
#'   formula_sigma = ~ treat * group,
#'   treat = "treat",
#'   B = 500,
#'   seed = 42
#' )
#' result <- ineqx(params, se_method = "bootstrap", boot = bc)
#' }
#'
#' @export
boot_config <- function(data, formula_mu, formula_sigma, treat,
                         group = NULL, time = NULL,
                         post = NULL, B = 200L, family = NULL,
                         parallel = FALSE, ncores = NULL,
                         seed = NULL, verbose = TRUE,
                         gamlss_args = list()) {

  if (!is.data.frame(data)) stop("'data' must be a data.frame")
  if (!inherits(formula_mu, "formula")) stop("'formula_mu' must be a formula")
  if (!inherits(formula_sigma, "formula")) stop("'formula_sigma' must be a formula")
  if (!is.character(treat) || length(treat) != 1) stop("'treat' must be a single character string")
  if (!treat %in% names(data)) stop("'", treat, "' not found in data")
  if (B < 2) stop("'B' must be at least 2")

  structure(
    list(
      data = data,
      formula_mu = formula_mu,
      formula_sigma = formula_sigma,
      treat = treat,
      group = group,
      time = time,
      post = post,
      B = as.integer(B),
      family = family,
      parallel = parallel,
      ncores = ncores,
      seed = seed,
      verbose = verbose,
      gamlss_args = gamlss_args
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
#' @param formula_mu Formula for the mean (mu) equation
#' @param formula_sigma Formula for the log-SD (sigma) equation
#' @param treat Character, treatment variable name
#' @param group Character, grouping variable name
#' @param time Character, time variable name. NULL for cross-sectional.
#' @param post Character, pre/post indicator for DiD. NULL for simple diff.
#' @param ref Numeric, reference time period for longitudinal decomposition
#' @param ystat Character, \code{"Var"} or \code{"CV2"}
#' @param order Character vector of length 3, decomposition ordering
#' @param B Integer, number of bootstrap replicates. Default 200.
#' @param family A gamlss.family object. NULL uses normal distribution.
#' @param parallel Logical, use parallel computation. Default FALSE.
#' @param ncores Integer, number of cores. Default: all but one.
#' @param seed Integer, random seed. Default NULL.
#' @param verbose Logical, print progress. Default TRUE.
#' @param gamlss_args Named list of additional arguments for gamlss().
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
                          B = 200L, family = NULL,
                          parallel = FALSE, ncores = NULL,
                          seed = NULL, verbose = TRUE,
                          gamlss_args = list()) {

  if (!requireNamespace("gamlss", quietly = TRUE)) {
    stop("Package 'gamlss' is required for bootstrap_se(). ",
         "Install it with install.packages('gamlss')")
  }

  ystat <- match.arg(ystat, c("Var", "CV2"))
  B <- as.integer(B)
  if (B < 2) stop("'B' must be at least 2")

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
    ref = ref, ystat = ystat, order = order, family = family,
    gamlss_args = gamlss_args, resample_ids = seq_len(nrow(data))
  )
  if (is.null(orig_result)) {
    stop("Original model estimation failed. Cannot compute bootstrap SEs.")
  }

  # Run bootstrap replicates
  if (parallel) {
    results <- .boot_parallel(
      data = data, formula_mu = formula_mu, formula_sigma = formula_sigma,
      treat = treat, group = group, time = time, post = post,
      ref = ref, ystat = ystat, order = order, family = family,
      gamlss_args = gamlss_args, boot_indices = boot_indices, ncores = ncores
    )
  } else {
    results <- vector("list", B)
    report_interval <- max(1L, B %/% 10L)
    for (b in seq_len(B)) {
      if (verbose && b %% report_interval == 0) {
        message(sprintf("Bootstrap: %d%% complete (%d/%d)",
                        round(100 * b / B), b, B))
      }
      results[[b]] <- tryCatch(
        .boot_one_replicate(
          data = data,
          formula_mu = formula_mu, formula_sigma = formula_sigma,
          treat = treat, group = group, time = time, post = post,
          ref = ref, ystat = ystat, order = order, family = family,
          gamlss_args = gamlss_args, resample_ids = boot_indices[[b]]
        ),
        error = function(e) NULL
      )
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
    cat("  SE(delta_B):    ", round(se$se_delta_B, 4), "\n")
    cat("  SE(delta_W):    ", round(se$se_delta_W, 4), "\n")
    cat("  SE(delta_total):", round(se$se_delta_total, 4), "\n")
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
                                 ref, ystat, order, family,
                                 gamlss_args, resample_ids) {

  # Resample
  boot_data <- data[resample_ids, , drop = FALSE]

  # Fit GAMLSS
  fit_args <- c(
    list(
      formula = formula_mu,
      sigma.formula = formula_sigma,
      data = boot_data,
      trace = FALSE
    ),
    gamlss_args
  )
  if (!is.null(family)) fit_args$family <- family

  model <- do.call(gamlss::gamlss, fit_args)

  # Extract parameters (no vcov needed)
  params <- extract_params(
    model = model, treat = treat, group = group,
    time = time, post = post, data = boot_data,
    ref = ref, ystat = ystat, vcov = FALSE
  )

  # Compute decomposition
  if (params$type == "cross_sectional") {
    result <- causal_decompose_cross(params)
    c(delta_B = result$delta_B,
      delta_W = result$delta_W,
      delta_total = result$delta_total)
  } else {
    result <- causal_decompose_longit(params, order = order)
    # Flatten into named vector
    estimates <- c()
    for (t_name in names(result$results)) {
      r <- result$results[[t_name]]
      t_ests <- c(
        Delta_beta   = r$Delta_beta,
        Delta_lambda = r$Delta_lambda,
        Delta_pi_B   = r$Delta_pi_B,
        Delta_pi_W   = r$Delta_pi_W,
        Delta_mu     = r$Delta_mu,
        Delta_sigma  = r$Delta_sigma,
        Delta_behavioral    = r$Delta_behavioral,
        Delta_compositional = r$Delta_compositional,
        Delta_pretreatment  = r$Delta_pretreatment,
        Delta_B      = r$Delta_B,
        Delta_W      = r$Delta_W,
        Delta_total  = r$Delta_total
      )
      names(t_ests) <- paste0(names(t_ests), "_", t_name)
      estimates <- c(estimates, t_ests)
    }
    estimates
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
                            ref, ystat, order, family,
                            gamlss_args, boot_indices, ncores) {

  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("Package 'parallel' is required for parallel bootstrap")
  }
  if (is.null(ncores)) ncores <- max(1L, parallel::detectCores() - 1L)

  cl <- parallel::makeCluster(ncores)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  # Export data and settings to workers
  parallel::clusterExport(cl, c(
    "data", "formula_mu", "formula_sigma", "treat", "group",
    "time", "post", "ref", "ystat", "order", "family", "gamlss_args"
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
        order = order, family = family,
        gamlss_args = gamlss_args, resample_ids = ids
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
      se_delta_B = unname(boot_sds["delta_B"]),
      se_delta_W = unname(boot_sds["delta_W"]),
      se_delta_total = unname(boot_sds["delta_total"])
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
