# ============================================================================ #
# ineqx: Unified variance decomposition (descriptive + causal)
# ============================================================================ #

#' Variance decomposition
#'
#' A unified function for both descriptive and causal variance decomposition.
#' If \code{treat} is not specified (and no \code{params} are provided), a
#' descriptive within/between decomposition is performed. If \code{treat} is
#' specified or \code{params} are provided, a causal decomposition of the
#' treatment effect on inequality is performed.
#'
#' There are five usage modes:
#' \describe{
#'   \item{Descriptive (raw data)}{Provide \code{y} but omit \code{treat} and
#'     \code{params}. Decomposes total inequality into within- and between-group
#'     components. For longitudinal data with a reference period, further
#'     decomposes the change over time into contributions from changing means,
#'     dispersions, and group composition.}
#'   \item{Descriptive (counterfactual reference)}{Provide \code{y}, \code{data},
#'     and an \code{ineqx_desc_params} object via \code{params}. The params
#'     define a counterfactual baseline; observed group-level statistics are
#'     computed from the data. The decomposition shows how each parameter
#'     contributes to the change from counterfactual to observed inequality.}
#'   \item{Causal (integrated estimation)}{Provide \code{treat}, \code{formula_mu},
#'     \code{formula_sigma}, and \code{data}. The function fits a GAMLSS model,
#'     extracts parameters, and performs the causal decomposition.}
#'   \item{Causal (counterfactual reference)}{Provide an \code{ineqx_params}
#'     object via \code{params} together with \code{treat}, \code{formula_mu},
#'     and \code{formula_sigma}. The params define a counterfactual baseline
#'     (e.g., zero treatment effects), and a GAMLSS model is fitted to estimate
#'     parameters for observed periods. The two are blended for a longitudinal
#'     decomposition showing how treatment effects changed from the baseline.
#'     Delta method SEs are not available; use \code{boot_config()} for
#'     bootstrap SEs.}
#'   \item{Causal (externally estimated model)}{Provide an \code{ineqx_params}
#'     object via \code{params} (without \code{treat}/formulas). Use
#'     \code{\link{ineqx_params}} to create one manually or from a fitted
#'     gamlss model.}
#' }
#'
#' @param y Character, name of the outcome variable in \code{data}.
#' @param ystat Character, one of \code{"Var"} (default), \code{"CV2"}, or
#'   \code{"VL"} (variance of log). \code{"VL"} is supported for descriptive
#'   decomposition and for integrated causal estimation (when \code{treat},
#'   \code{formula_mu}, and \code{formula_sigma} are supplied with raw
#'   \code{y}/\code{data}). It is implemented by log-transforming \code{y}
#'   and running the standard \code{"Var"} decomposition; output is labelled
#'   as \code{"VL"}. \code{"VL"} is not supported when \code{params} is
#'   supplied; in that case, fit your model on \code{log(y)} yourself and
#'   pass \code{ystat = "Var"}.
#' @param treat Character, name of the treatment variable in \code{data}.
#'   Coded 0/1. If NULL, a descriptive decomposition is performed.
#' @param post Character, pre/post indicator for DiD designs. NULL for
#'   simple difference estimator. Only used in integrated estimation mode.
#'   For \code{ystat = "CV2"}, keep \code{y} on the outcome's level scale
#'   and use \code{post} for the DiD contrast; applying \code{"CV2"} to
#'   first-differenced outcomes is not recommended because it targets
#'   relative dispersion in changes and can be unstable when mean changes
#'   are near zero.
#' @param group Character, name of the grouping variable in \code{data}.
#' @param time Character, name of the time variable in \code{data}. If NULL,
#'   a single cross-section is assumed.
#' @param ref Numeric, reference time period. Required for longitudinal
#'   decomposition (both descriptive and causal).
#' @param order Decomposition ordering. Either:
#'   \itemize{
#'     \item \code{"shapley"} (default): averages across all 6 possible orderings.
#'     \item For descriptive: a permutation of \code{c("mu", "sigma", "pi")}.
#'     \item For causal: a permutation of
#'       \code{c("behavioral", "compositional", "pretreatment")}.
#'   }
#'   Only relevant for longitudinal data with a reference period.
#' @param formula_mu One-sided formula for the mean equation (integrated
#'   estimation mode). E.g., \code{~ treat * group + controls}. The outcome
#'   \code{y} is prepended internally.
#' @param formula_sigma One-sided formula for the log-SD equation (integrated
#'   estimation mode). E.g., \code{~ treat * group + controls}.
#' @param estimand Character, either \code{"marginal"} (default) or
#'   \code{"residual"}. Selects whether within-group dispersion is the marginal
#'   counterfactual variance (law of total variance over the covariate
#'   distribution; controls contribute to within-group inequality) or the
#'   residual/conditional scale (paper Appendix B.7). See
#'   \code{\link{ineqx_params}}. Only used in the integrated- and
#'   blending-estimation paths.
#' @param params A parameter object created by \code{\link{ineqx_params}}.
#'   Either an \code{ineqx_desc_params} (descriptive counterfactual reference)
#'   or an \code{ineqx_params} (causal decomposition). When an
#'   \code{ineqx_params} is provided together with \code{treat},
#'   \code{formula_mu}, and \code{formula_sigma}, blending mode is used:
#'   the params define a counterfactual baseline, and a GAMLSS model is
#'   fitted to estimate parameters for observed periods. See Details.
#' @param weights Character, name of the weight variable in \code{data}.
#'   If NULL, equal weights are used.
#' @param se Standard error method. One of:
#'   \itemize{
#'     \item \code{"delta"} or \code{TRUE} (default): delta method SEs. For
#'       causal, requires vcov in params. For descriptive, uses sampling-based
#'       covariance.
#'     \item \code{"none"} or \code{FALSE}: skip SE computation
#'     \item \code{"boot"}: bootstrap SEs with default settings (equivalent to
#'       \code{se = boot_config()})
#'     \item A \code{\link{boot_config}} object: bootstrap SEs with custom
#'       settings (causal only)
#'   }
#' @param data A data.frame containing the variables. For manual descriptive
#'   mode (when \code{y = NULL}), must contain columns \code{group}, \code{pi},
#'   \code{mu}, \code{sigma}, and optionally a time column.
#'
#' @return For descriptive decomposition: an \code{ineqx_desc} object.
#'   For cross-sectional causal: an \code{ineqx_causal_cross} object.
#'   For longitudinal causal (both specific orderings and Shapley): an
#'   \code{ineqx_causal_longit} object. Shapley results include additional
#'   fields \code{$shapley}, \code{$all_orderings}, and \code{$ranges}.
#'
#' @examples
#' data(incdat)
#'
#' # Descriptive decomposition from raw data
#' ineqx("inc", group = "group", time = "year", data = incdat, ref = 1)
#'
#' # Descriptive with counterfactual reference
#' ref <- ineqx_params(data = data.frame(
#'   group = 1:3, pi = 1/3, mu = 0, sigma = 0
#' ))
#' ineqx("inc", group = "group", time = "year",
#'       params = ref, data = incdat)
#'
#' # Causal: from manual params (externally estimated model)
#' params <- ineqx_params(
#'   data = data.frame(
#'     group = c("workers", "managers"),
#'     pi = c(0.5, 0.5),
#'     mu0 = c(500, 1000),
#'     sigma0 = c(200, 400),
#'     beta = c(60, 100),
#'     lambda = c(-0.1, -0.2)
#'   )
#' )
#' ineqx("inc", group = "group", data = incdat, params = params, se = "none")
#'
#' \dontrun{
#' # Causal: integrated estimation
#' ineqx("inc", treat = "x", group = "group", data = incdat,
#'       formula_mu = ~ x * factor(group),
#'       formula_sigma = ~ x * factor(group),
#'       se = "delta")
#'
#' # Causal: counterfactual reference (blending)
#' cf_ref <- ineqx_params(data = data.frame(
#'   group = 1:3, pi = 1/3,
#'   mu0 = 500, sigma0 = 100, beta = 0, lambda = 0
#' ))
#' ineqx("inc", treat = "x", group = "group",
#'       params = cf_ref,
#'       formula_mu = ~ x * factor(group),
#'       formula_sigma = ~ x * factor(group),
#'       se = "none", data = incdat)
#' }
#'
#' @param ... Additional arguments passed to \code{gamlss::gamlss()} for the
#'   integrated- and blending-estimation paths. Useful for bumping the
#'   iteration limit on saturated models (\code{n.cyc = 100}) or selecting a
#'   non-default family.
#' @export
ineqx <- function(y = NULL, ystat = "Var", treat = NULL, post = NULL,
                  group = NULL, time = NULL, ref = NULL, order = "shapley",
                  formula_mu = NULL, formula_sigma = NULL,
                  estimand = c("marginal", "residual"),
                  params = NULL, weights = NULL, se = "delta", data = NULL,
                  ...) {

  ystat <- match.arg(ystat, c("Var", "CV2", "VL"))
  estimand <- match.arg(estimand)

  # -------------------------------------------------------------------- #
  # VL: variance of log(y). Implemented by running the standard Var
  # decomposition on a log-scale fit; output is labelled "VL".
  #
  # Two code paths produce VL:
  #   1. Integrated path (y + data supplied): we log-transform data[[y]]
  #      here and proceed with ystat = "Var" internals.
  #   2. Split path (params supplied): the params object must come from a
  #      fit on log(y) -- detected via params$transform == "log", which is
  #      set by fit_ineqx_model(transform = "log") and propagated through
  #      ineqx_params(). In that case the parameters already live on the
  #      log scale, so we just switch ystat to "Var" and relabel on output.
  # -------------------------------------------------------------------- #
  ystat_label <- ystat
  if (ystat == "VL") {
    if (!is.null(params)) {
      params_transform <- params$transform %||% "identity"
      if (!identical(params_transform, "log")) {
        stop("ystat = 'VL' with params requires a params object built from ",
             "a log(y) fit. The supplied params has transform = '",
             params_transform, "'. To compute V_L with the split-step ",
             "workflow, fit with fit_ineqx_model(..., transform = 'log'), ",
             "extract params via ineqx_params(model = ...), then call ",
             "ineqx(params = ..., ystat = 'VL'). Note: because V_L ",
             "measures dispersion in log earnings rather than in earnings ",
             "themselves, its decomposition can diverge from income-scale ",
             "changes in inequality (see Rosche 2026).",
             call. = FALSE)
      }
      warning("Because V_L measures dispersion in log earnings rather than in ",
              "earnings themselves, its decomposition can diverge from ",
              "income-scale changes in inequality (see Rosche 2026).",
              call. = FALSE)
      ystat <- "Var"   # params already on log scale; run internals as Var
    } else {
      if (is.null(y)) {
        stop("ystat = 'VL' requires 'y' so it can be log-transformed.")
      }
      if (missing(data) || is.null(data) || !y %in% names(data)) {
        stop("ystat = 'VL' requires 'data' containing column '", y, "'.")
      }
      yvals <- data[[y]]
      bad <- !is.na(yvals) & yvals <= 0
      if (any(bad)) {
        stop("ystat = 'VL' requires strictly positive '", y, "'; found ",
             sum(bad), " non-positive value(s).")
      }
      warning("Because V_L measures dispersion in log earnings rather than in ",
              "earnings themselves, its decomposition can diverge from ",
              "income-scale changes in inequality (see Rosche 2026).",
              call. = FALSE)
      data[[y]] <- log(yvals)
      ystat <- "Var"   # run internals on Var of log(y)
    }
  }

  # -------------------------------------------------------------------- #
  # Parse se argument (shared across descriptive and causal)
  # -------------------------------------------------------------------- #

  parsed_se <- .parse_se(se)
  se_method <- parsed_se$method
  boot <- parsed_se$boot

  # -------------------------------------------------------------------- #
  # Dispatch: Descriptive vs Causal
  # -------------------------------------------------------------------- #

  if (!is.null(params) && inherits(params, "ineqx_desc_params")) {
    # ================================================================== #
    # DESCRIPTIVE WITH COUNTERFACTUAL REFERENCE
    # ================================================================== #

    if (!is.null(treat)) {
      stop("Cannot use descriptive params with 'treat'. ",
           "For causal decomposition, provide causal params ",
           "(mu0, sigma0, beta, lambda).")
    }
    if (is.null(y)) {
      stop("'y' is required when using descriptive params with observed data.")
    }

    # Auto-default ref for single cross-section params
    if (is.null(params$times) && is.null(ref)) {
      ref <- 0L
    }
    # Multi-period params require explicit ref
    if (!is.null(params$times) && is.null(ref)) {
      stop("'ref' is required when params has multiple time periods. ",
           "Available times in params: ",
           paste(params$times, collapse = ", "))
    }

    message("Computing descriptive decomposition...")
    result <- .ineq_descriptive_with_ref(
      y = y, group = group, time = time, weights = weights,
      data = data, ref_params = params, ref = ref,
      ystat = ystat, order = order)

    # Compute delta method SEs for descriptive
    if (se_method == "delta") {
      result$se <- delta_method_desc_se(result$wibe, result$totals,
                                         result$deltas, ystat, ref = ref,
                                         order = order)
      result$se_method <- "delta"
    } else {
      result$se <- NULL
      result$se_method <- "none"
    }

    result$ystat <- ystat_label
    message("Finished.\n")
    return(result)

  } else if (is.null(treat) && is.null(params)) {
    # ================================================================== #
    # DESCRIPTIVE DECOMPOSITION (from raw data)
    # ================================================================== #

    if (is.null(y)) {
      stop("'y' is required for descriptive decomposition. ",
           "To use pre-computed parameters, pass them via ",
           "params = ineqx_params(data = ...).")
    }

    message("Computing descriptive decomposition...")
    result <- .ineq_descriptive(y = y, group = group, time = time,
                                weights = weights, data = data, ref = ref,
                                ystat = ystat, order = order)

    # Compute delta method SEs for descriptive
    if (se_method == "delta") {
      result$se <- delta_method_desc_se(result$wibe, result$totals,
                                         result$deltas, ystat, ref = ref,
                                         order = order)
      result$se_method <- "delta"
    } else {
      result$se <- NULL
      result$se_method <- "none"
    }

    result$ystat <- ystat_label
    message("Finished.\n")
    return(result)
  }

  # y is required for causal decomposition, *except* for the params-only path
  # where the user supplies pre-computed (mu0, sigma0, beta, lambda) and no
  # GAMLSS needs to be fit. In that case y/data are not used downstream.
  params_only <- !is.null(params) && inherits(params, "ineqx_params") &&
                  is.null(treat) && is.null(formula_mu) && is.null(formula_sigma)
  if (is.null(y) && !params_only) {
    stop("'y' is required for causal decomposition.")
  }

  # ================================================================== #
  # CAUSAL DECOMPOSITION
  # ================================================================== #

  # -------------------------------------------------------------------- #
  # Obtain params
  # -------------------------------------------------------------------- #

  blending <- FALSE
  user_params <- NULL

  if (!is.null(params) && inherits(params, "ineqx_params") &&
      !is.null(treat) && !is.null(formula_mu) && !is.null(formula_sigma)) {
    # ================================================================ #
    # BLENDING MODE: user params + model estimation
    # ================================================================ #
    blending <- TRUE
    user_params <- params

    # Auto-default ref for single cross-section user params
    if (is.null(user_params$times) && is.null(ref)) {
      ref <- 0L
    }
    # Multi-period user params require explicit ref
    if (!is.null(user_params$times) && is.null(ref)) {
      stop("'ref' is required when params has multiple time periods. ",
           "Available times in params: ",
           paste(user_params$times, collapse = ", "))
    }

    # Fit GAMLSS model on data
    message("Fitting GAMLSS model...")
    full_formula_mu <- .make_two_sided(formula_mu, y)
    data <- .complete_model_data(data, full_formula_mu, formula_sigma,
                                 y, group, treat, post, time, weights)
    model <- fit_ineqx_model(
      formula_mu = full_formula_mu,
      formula_sigma = formula_sigma,
      data = data,
      weights = weights,
      ...
    )

    # Extract model params (no vcov needed — delta method unavailable for blended)
    # Use data's first time period as ref for extraction (the real ref is set
    # on the blended params; this just needs to be a valid time in the data)
    message("Extracting parameters...")
    model_params <- ineqx_params(
      model = model, data = data,
      treat = treat, group = group,
      time = time, post = post,
      ystat = ystat, estimand = estimand,
      vcov = FALSE
    )

    # Free the fitted gamlss — bootstrap (if any) refits per replicate
    rm(model)
    gc(verbose = FALSE)

    # Blend user params with model params
    params <- .blend_causal_params(user_params, model_params, ref)

    # Override ystat from ineqx() call (user_params may have defaulted to "Var")
    params$ystat <- ystat

    # Delta method SE not available for blended params
    if (is.character(se) && se == "delta") {
      message("Delta method SEs not available for blended params; ",
              "falling back to se = 'none'. Use boot_config() for bootstrap SEs.")
      se <- "none"
    }

  } else if (!is.null(params)) {
    # Externally estimated model
    stopifnot(inherits(params, "ineqx_params"))

    # The explicit ystat= argument is the single source of truth. If the
    # supplied params object carries a different ystat (most often because
    # ineqx_params() defaulted to "Var" while the user wants "CV2"), override
    # it and surface the change so the user notices.
    if (!is.null(params$ystat) && !identical(params$ystat, ystat)) {
      message("Overriding params$ystat ('", params$ystat,
              "') with ineqx() argument ystat = '", ystat, "'.")
    }
    params$ystat <- ystat

  } else {
    # Integrated estimation
    if (is.null(formula_mu) || is.null(data) || is.null(group)) {
      stop("For integrated estimation, 'formula_mu', 'data', and 'group' ",
           "are required.\n",
           "Alternatively, pass pre-computed parameters via 'params'.")
    }
    if (is.null(formula_sigma)) {
      stop("'formula_sigma' is required for GAMLSS fitting.")
    }

    # Construct two-sided formula from y + one-sided formula_mu
    full_formula_mu <- .make_two_sided(formula_mu, y)

    # Drop rows with NAs in model variables
    data <- .complete_model_data(data, full_formula_mu, formula_sigma,
                                 y, group, treat, post, time, weights)

    # Fit GAMLSS
    message("Fitting GAMLSS model...")
    model <- fit_ineqx_model(
      formula_mu = full_formula_mu,
      formula_sigma = formula_sigma,
      data = data,
      weights = weights,
      ...
    )

    # Determine whether to extract vcov (needed for delta method SEs)
    extract_vcov <- se_method == "delta"

    # Extract params
    message("Extracting parameters...")
    params <- ineqx_params(
      model = model, data = data,
      treat = treat, group = group,
      time = time, post = post,
      ystat = ystat, estimand = estimand,
      vcov = extract_vcov
    )

    # Free the fitted gamlss — bootstrap (if any) refits per replicate
    rm(model)
    gc(verbose = FALSE)
  }

  # -------------------------------------------------------------------- #
  # Causal se: check vcov availability for delta method
  # -------------------------------------------------------------------- #

  if (se_method == "delta" && is.null(params$vcov)) {
    message("vcov not available in params; falling back to se = 'none'.")
    se_method <- "none"
  }

  # Suppress vcov-based SEs in decomposition if not using delta method
  params_for_decomp <- params
  if (se_method != "delta") {
    params_for_decomp$vcov <- NULL
  }

  # -------------------------------------------------------------------- #
  # Parse order argument
  # -------------------------------------------------------------------- #

  use_shapley <- is.character(order) && length(order) == 1 && order == "shapley"

  if (!use_shapley) {
    valid_components <- c("behavioral", "compositional", "pretreatment")
    order <- match.arg(order, valid_components, several.ok = TRUE)
    if (length(order) != 3 || !setequal(order, valid_components)) {
      stop("'order' must be 'shapley' or a permutation of ",
           "c('behavioral', 'compositional', 'pretreatment')")
    }
  }

  # -------------------------------------------------------------------- #
  # Dispatch decomposition
  # -------------------------------------------------------------------- #

  message("Computing decomposition...")
  if (params$type == "cross_sectional") {
    # Cross-sectional: order is ignored
    result <- causal_decompose_cross(params_for_decomp, ref = ref)

  } else if (use_shapley) {
    # Longitudinal with Shapley averaging
    result <- causal_shapley(params_for_decomp, ref = ref)

  } else {
    # Longitudinal with specific ordering
    result <- causal_decompose_longit(params_for_decomp, order = order, ref = ref)
  }

  # -------------------------------------------------------------------- #
  # Bootstrap SEs (if requested)
  # -------------------------------------------------------------------- #

  if (se_method == "bootstrap") {
    message("Computing bootstrap standard errors...")

    boot_result <- bootstrap_se(
      data = data,
      formula_mu = full_formula_mu,
      formula_sigma = formula_sigma,
      treat = treat,
      group = group,
      time = time,
      post = post,
      ref = ref,
      ystat = params$ystat,
      estimand = estimand,
      order = if (use_shapley) c("behavioral", "compositional", "pretreatment") else order,
      B = boot$B,
      parallel = boot$parallel,
      ncores = boot$ncores,
      seed = boot$seed,
      verbose = boot$verbose,
      cl_type = boot$cl_type,
      blend_params = user_params
    )

    result$se <- boot_result$se
    result$se_method <- "bootstrap"
    result$boot <- boot_result
  } else if (se_method == "delta") {
    result$se_method <- "delta"
  } else {
    result$se <- NULL
    result$se_method <- "none"
  }

  # Restore original params (with vcov) so that plot(type="params", ci=TRUE)

  # can compute parameter-level CIs even when decomposition SEs used bootstrap
  # (params_for_decomp may have had vcov stripped to suppress delta-method SEs)
  result$params <- params

  # Relabel ystat for VL: internals ran on Var of log(y); surface "VL" to callers.
  result$ystat <- ystat_label

  message("Finished.\n")
  result
}


# ---------------------------------------------------------------------------- #
# Internal: Blend user-provided and model-estimated causal params
# ---------------------------------------------------------------------------- #

.blend_causal_params <- function(user_params, model_params, ref) {

  # Validate group alignment
  if (!setequal(user_params$groups, model_params$groups)) {
    stop("Groups in 'params' (", paste(user_params$groups, collapse = ", "),
         ") do not match groups in data (",
         paste(model_params$groups, collapse = ", "), ")")
  }

  # Determine user params times
  u_data <- user_params$data
  if (is.null(user_params$times)) {
    # Single cross-section: assign time=0
    u_data$time <- 0L
    user_times <- 0L
  } else {
    user_times <- user_params$times
  }

  # Model params times
  m_data <- model_params$data
  model_times <- model_params$times
  if (is.null(model_times)) {
    # Single CS model (no time column) — shouldn't happen typically, but handle
    m_data$time <- 1L
    model_times <- 1L
  }

  # Merge: user periods override model periods on overlap
  model_only_times <- setdiff(model_times, user_times)
  m_data_keep <- m_data[m_data$time %in% model_only_times, ]

  # Align columns before rbind: user-supplied params may not include the
  # DiD-specific or extended columns (mu1, sigma1, mu0_pre, ...). Pad missing
  # columns on either side with NA so rbind works regardless of which path
  # produced each side.
  all_cols <- union(names(u_data), names(m_data_keep))
  for (col in setdiff(all_cols, names(u_data))) u_data[[col]] <- NA_real_
  for (col in setdiff(all_cols, names(m_data_keep))) m_data_keep[[col]] <- NA_real_
  u_data      <- u_data[, all_cols, drop = FALSE]
  m_data_keep <- m_data_keep[, all_cols, drop = FALSE]

  merged_data <- rbind(u_data, m_data_keep)
  merged_data <- merged_data[order(merged_data$time, merged_data$group), ]
  rownames(merged_data) <- NULL

  all_times <- sort(unique(merged_data$time))

  # Merge vcov: keep model vcov for model-only periods
  merged_vcov <- NULL
  if (!is.null(model_params$vcov)) {
    if (is.list(model_params$vcov)) {
      # Longitudinal: named list of matrices per time
      merged_vcov <- model_params$vcov[as.character(model_only_times)]
      if (length(merged_vcov) == 0) merged_vcov <- NULL
    } else if (is.matrix(model_params$vcov) && length(model_only_times) == 1) {
      # Single model period remaining
      merged_vcov <- list()
      merged_vcov[[as.character(model_only_times)]] <- model_params$vcov
    }
  }

  # Validate ref against merged times
  if (!(ref %in% all_times)) {
    stop("'ref' = ", ref, " not found in merged time periods. ",
         "Available: ", paste(all_times, collapse = ", "))
  }

  structure(
    list(
      data = merged_data,
      ystat = user_params$ystat,
      estimand = model_params$estimand %||% "marginal",
      type = "longitudinal",
      vcov = merged_vcov,
      groups = sort(unique(merged_data$group)),
      times = all_times,
      n_groups = length(unique(merged_data$group)),
      n_times = length(all_times),
      blended = TRUE,
      user_times = user_times,
      model_times = model_only_times,
      is_did = isTRUE(model_params$is_did),
      post   = model_params$post
    ),
    class = "ineqx_params"
  )
}


# ---------------------------------------------------------------------------- #
# Internal: Construct two-sided formula from y + one-sided formula
# ---------------------------------------------------------------------------- #

.make_two_sided <- function(formula, y) {
  if (length(formula) == 3) {
    # Already two-sided, return as-is
    return(formula)
  }
  as.formula(paste(y, "~", paste(deparse(formula[[2]]), collapse = " ")))
}


# ---------------------------------------------------------------------------- #
# Internal: Drop rows with NAs in model-relevant variables only
# ---------------------------------------------------------------------------- #

.complete_model_data <- function(data, formula_mu, formula_sigma,
                                 y, group, treat = NULL, post = NULL,
                                 time = NULL, weights = NULL) {
  vars <- unique(c(
    all.vars(formula_mu), all.vars(formula_sigma),
    y, group,
    if (!is.null(treat)) treat,
    if (!is.null(post)) post,
    if (!is.null(time)) time,
    if (!is.null(weights) && is.character(weights)) weights
  ))
  vars <- intersect(vars, names(data))
  data <- data[, vars, drop = FALSE]
  complete <- complete.cases(data)
  n_dropped <- sum(!complete)
  if (n_dropped > 0) {
    data <- data[complete, , drop = FALSE]
    message("N = ", nrow(data),
            " (", n_dropped, " rows omitted because of missing values).")
  } else {
    message("N = ", nrow(data), ".")
  }
  data
}


# ---------------------------------------------------------------------------- #
# Internal: Parse se argument (shared by descriptive and causal paths)
# ---------------------------------------------------------------------------- #

.parse_se <- function(se) {
  if (isFALSE(se)) {
    return(list(method = "none", boot = NULL))
  }
  if (isTRUE(se)) {
    return(list(method = "delta", boot = NULL))
  }
  if (inherits(se, "ineqx_boot_config")) {
    return(list(method = "bootstrap", boot = se))
  }
  if (is.character(se)) {
    if (identical(se, "boot")) {
      return(list(method = "bootstrap", boot = boot_config()))
    }
    se <- match.arg(se, c("delta", "none"))
    return(list(method = se, boot = NULL))
  }
  stop("'se' must be TRUE, 'delta', 'none', 'boot', or a boot_config() object.")
}
