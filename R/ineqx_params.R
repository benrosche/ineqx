# ============================================================================ #
# ineqx_params: Standardized input for causal variance decomposition
# ============================================================================ #

#' Create an ineqx_params object
#'
#' Constructs a standardized parameter object that serves as the interface
#' between model estimation and variance decomposition. This decouples the
#' two stages: users can construct \code{ineqx_params} from any estimation
#' method (GAMLSS, Bayesian models, simulation, or manual specification).
#'
#' There are three usage modes:
#' \describe{
#'   \item{Descriptive (manual)}{Pass a data.frame with columns
#'     \code{group}, \code{pi}, \code{mu}, \code{sigma}. Returns an
#'     \code{ineqx_desc_params} object representing a counterfactual
#'     reference for descriptive decomposition.}
#'   \item{Causal (manual)}{Pass a data.frame with columns
#'     \code{group}, \code{pi}, \code{mu0}, \code{sigma0}, \code{beta},
#'     \code{lambda} (and optionally \code{time}). Returns an
#'     \code{ineqx_params} object.}
#'   \item{Causal (from fitted gamlss)}{Pass a fitted \code{gamlss} object
#'     via \code{model}, along with the raw data and variable names. The
#'     function extracts treatment effects via counterfactual predictions.}
#' }
#'
#' The mode is auto-detected from the columns in \code{data}: if \code{mu}
#' and \code{sigma} are present (without \code{mu0}, \code{sigma0}), the
#' descriptive path is used; otherwise the causal path is used.
#'
#' @param data For descriptive mode: a data.frame with columns
#'   \code{group}, \code{pi}, \code{mu}, \code{sigma}.
#'   For causal manual mode: a data.frame with columns
#'   \code{group}, \code{pi}, \code{mu0}, \code{sigma0}, \code{beta},
#'   \code{lambda} (and optionally \code{time}).
#'   For model mode: the data.frame used to fit the model.
#' @param model A fitted \code{gamlss} object. If provided, parameters are
#'   extracted from the model via counterfactual predictions. Default: NULL
#'   (manual specification).
#' @param treat Character, name of the treatment variable in \code{data}.
#'   Required when \code{model} is provided. Must be coded as 0 (untreated)
#'   and 1 (treated), or with multiple non-zero treatment levels.
#' @param group Character, name of the grouping variable in \code{data}.
#'   Required when \code{model} is provided.
#' @param time Character, name of the time variable. NULL for cross-sectional.
#'   Used in both manual mode (if \code{data} contains a \code{time} column)
#'   and model mode.
#' @param post Character, name of the pre/post indicator for DiD designs.
#'   Only used in model mode. NULL for simple difference estimator. For
#'   \code{ystat = "CV2"}, keep the outcome in levels and use \code{post}
#'   to construct the DiD contrast; applying \code{"CV2"} to first
#'   differences is not recommended.
#' @param ystat Character, either \code{"Var"} (variance) or \code{"CV2"}
#'   (squared coefficient of variation). Default: \code{"Var"}.
#' @param estimand Character, either \code{"marginal"} (default) or
#'   \code{"residual"}. The \code{"marginal"} estimand computes each group's
#'   variance by the law of total variance over the covariate distribution of
#'   the treated (paper eqs 10, 37-38): the average conditional residual
#'   variance plus the dispersion of predicted conditional means across
#'   covariate profiles. Controls therefore contribute to within-group
#'   inequality. The \code{"residual"} variant (paper Appendix B.7) uses the
#'   conditional scale parameter directly and omits the predicted-mean
#'   dispersion term, answering a residual-inequality question instead. In the
#'   manual (data.frame) mode this only tags the object; the supplied
#'   \code{sigma0}/\code{lambda} are taken as marginal \eqn{S}/\eqn{\Lambda}
#'   (or residual \eqn{\sigma}/\eqn{\lambda} when \code{"residual"}).
#' @param vcov For manual mode: an optional variance-covariance matrix (or
#'   named list of matrices for longitudinal data).
#'   For model mode: logical, whether to extract the vcov via numerical
#'   Jacobian (default: TRUE). Set to FALSE for faster computation without SEs.
#' @param na.action A function that handles NAs in \code{data}. Only used in
#'   model mode. Applied to the subset of \code{data} containing the variables
#'   referenced by the model's mu/sigma formulas plus \code{treat}, \code{group},
#'   \code{time}, and \code{post}. Defaults to \code{\link[stats]{na.omit}},
#'   matching the integrated \code{\link{ineqx}()} path. Pass \code{na.fail} to
#'   error on NAs, or \code{na.pass} to keep them (and let \code{predict()}
#'   propagate them).
#'
#' @return An object of class \code{"ineqx_desc_params"} (descriptive mode) or
#'   \code{"ineqx_params"} (causal mode)
#'
#' @details
#' When \code{model} is provided, the function generates predictions under
#' counterfactual treatment assignments. For each observation, it predicts
#' the outcome under treat=0 (baseline) and treat=1 (treated), for both the
#' mean (mu) and standard deviation (sigma) equations. The average marginal
#' effects are computed within each group-time cell.
#'
#' For the simple difference estimator, the ATT is:
#' \deqn{\beta_{D,j} = E[\hat{\mu}(D=1) - \hat{\mu}(D=0) | G=j]}
#'
#' @examples
#' # Descriptive counterfactual reference (equal groups, no inequality)
#' desc_ref <- ineqx_params(
#'   data = data.frame(
#'     group = c("workers", "managers"),
#'     pi = c(0.5, 0.5),
#'     mu = c(0, 0),
#'     sigma = c(0, 0)
#'   )
#' )
#'
#' # Causal manual specification for a two-group example
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
#'
#' \dontrun{
#' # From a fitted gamlss model
#' params <- ineqx_params(
#'   model = my_gamlss, data = mydata,
#'   treat = "x", group = "SES",
#'   ystat = "Var", vcov = TRUE
#' )
#'
#' # Stepwise DiD: pass `post`, the indicator for the post-treatment
#' # observation in a paired pre/post panel. ineqx_params() then runs the
#' # 4-prediction DiD extraction (parallel-trends adjusted) instead of the
#' # 2-prediction simple-difference extraction.
#' did_params <- ineqx_params(
#'   model = my_did_gamlss, data = panel_data,
#'   treat = "treat", group = "edu", time = "year5",
#'   post = "post01", ystat = "Var", vcov = TRUE
#' )
#'
#' # V_L (variance of log earnings) via the split-step workflow:
#' log_model  <- fit_ineqx_model(
#'   formula_mu    = y ~ treat * group * time,
#'   formula_sigma =   ~ treat * group * time,
#'   data = mydata, transform = "log"
#' )
#' log_params <- ineqx_params(model = log_model, data = mydata,
#'                            treat = "treat", group = "group", time = "time",
#'                            ystat = "Var", vcov = TRUE)
#' fit_vl <- ineqx(params = log_params, ystat = "VL", ref = 1980)
#' }
#'
#' @export
ineqx_params <- function(data, model = NULL, treat = NULL, group = NULL,
                          time = NULL, post = NULL,
                          ystat = "Var", estimand = c("marginal", "residual"),
                          vcov = NULL,
                          na.action = stats::na.omit,
                          verbose = TRUE) {

  estimand <- match.arg(estimand)

  if (identical(ystat, "VL")) {
    stop("ystat = 'VL' is not supported in ineqx_params(). To compute V_L ",
         "with the split-step workflow, fit your model with ",
         "fit_ineqx_model(..., transform = 'log'), then call ",
         "ineqx_params(model = ..., ystat = 'Var') and pass the result to ",
         "ineqx(params = ..., ystat = 'VL'). The transform tag is ",
         "propagated automatically so ineqx() recognises the params as ",
         "living on the log scale. Note: because V_L measures dispersion ",
         "in log earnings rather than in earnings themselves, its ",
         "decomposition can diverge from income-scale changes in ",
         "inequality (see Rosche 2026).",
         call. = FALSE)
  }
  ystat <- match.arg(ystat, c("Var", "CV2"))

  if (!is.null(model)) {
    # -------------------------------------------------------------------- #
    # MODEL EXTRACTION PATH
    # -------------------------------------------------------------------- #
    .ineqx_params_from_model(
      model = model, data = data,
      treat = treat, group = group,
      time = time, post = post,
      ystat = ystat, estimand = estimand, vcov = vcov,
      na.action = na.action
    )
  } else {
    # -------------------------------------------------------------------- #
    # MANUAL SPECIFICATION PATH
    # -------------------------------------------------------------------- #

    # Auto-detect descriptive vs causal based on columns
    causal_cols <- c("mu0", "sigma0", "beta", "lambda")
    desc_cols <- c("mu", "sigma")

    if (is.data.frame(data) &&
        all(desc_cols %in% names(data)) &&
        !any(causal_cols %in% names(data))) {
      .ineqx_desc_params_manual(data = data, ystat = ystat)
    } else {
      .ineqx_params_manual(data = data, ystat = ystat, vcov = vcov,
                           estimand = estimand)
    }
  }
}


# ---------------------------------------------------------------------------- #
# Internal: Manual specification path
# ---------------------------------------------------------------------------- #

.ineqx_params_manual <- function(data, ystat, vcov,
                                 estimand = "marginal") {

  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame")
  }

  # Check required columns
  required_cols <- c("group", "pi", "mu0", "sigma0", "beta", "lambda")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in 'data': ",
         paste(missing_cols, collapse = ", "),
         "\nRequired columns: group, pi, mu0, sigma0, beta, lambda",
         "\n(and 'time' for longitudinal decomposition)",
         "\nDid you mean to pass a fitted model via model=?")
  }

  # Determine type
  has_time <- "time" %in% names(data)
  n_times <- if (has_time) length(unique(data$time)) else 1L
  type <- if (has_time && n_times > 1) "longitudinal" else "cross_sectional"

  # Validate numeric columns
  numeric_cols <- c("pi", "mu0", "sigma0", "beta", "lambda")
  for (col in numeric_cols) {
    if (!is.numeric(data[[col]])) {
      stop("Column '", col, "' must be numeric")
    }
    if (any(is.na(data[[col]]))) {
      stop("Column '", col, "' contains NA values")
    }
  }

  # Validate sigma0 >= 0
  if (any(data$sigma0 < 0)) {
    stop("'sigma0' must be non-negative (it represents the baseline SD)")
  }

  # Validate pi sums to 1 within each time period
  if (has_time) {
    times <- unique(data$time)
    for (t in times) {
      pi_sum <- sum(data$pi[data$time == t])
      if (abs(pi_sum - 1) > 1e-6) {
        stop("'pi' must sum to 1 within each time period. ",
             "At time = ", t, ", sum(pi) = ", round(pi_sum, 6))
      }
    }
  } else {
    pi_sum <- sum(data$pi)
    if (abs(pi_sum - 1) > 1e-6) {
      stop("'pi' must sum to 1. Currently sum(pi) = ", round(pi_sum, 6))
    }
  }

  # Validate no duplicate group-time combinations
  if (has_time) {
    gt <- paste(data$group, data$time, sep = "_")
  } else {
    gt <- as.character(data$group)
  }
  if (any(duplicated(gt))) {
    stop("Duplicate group-time combinations found in 'data'")
  }

  # Sort by time then group for consistency
  if (has_time) {
    data <- data[order(data$time, data$group), ]
  } else {
    data <- data[order(data$group), ]
  }

  # Validate vcov if provided
  n_groups <- length(unique(data$group))
  if (!is.null(vcov)) {
    expected_dim <- 4 * n_groups
    if (is.matrix(vcov)) {
      if (nrow(vcov) != expected_dim || ncol(vcov) != expected_dim) {
        stop("vcov matrix must be ", expected_dim, " x ", expected_dim,
             " (4 * J where J = ", n_groups, " groups). ",
             "Got ", nrow(vcov), " x ", ncol(vcov))
      }
    } else if (is.list(vcov)) {
      for (nm in names(vcov)) {
        V <- vcov[[nm]]
        if (!is.matrix(V) || nrow(V) != expected_dim || ncol(V) != expected_dim) {
          stop("vcov[['", nm, "']] must be a ", expected_dim, " x ", expected_dim,
               " matrix. Got ", nrow(V), " x ", ncol(V))
        }
      }
    } else {
      stop("vcov must be a matrix or a named list of matrices")
    }
  }

  .make_ineqx_params(data, ystat, type, vcov, has_time, estimand)
}


# ---------------------------------------------------------------------------- #
# Internal: Model extraction path
# ---------------------------------------------------------------------------- #

.ineqx_params_from_model <- function(model, data, treat, group,
                                      time, post, ystat,
                                      estimand = "marginal", vcov,
                                      na.action = stats::na.omit,
                                      verbose = TRUE) {

  if (!requireNamespace("gamlss", quietly = TRUE)) {
    stop("Package 'gamlss' is required when model is provided. ",
         "Install it with install.packages('gamlss')")
  }

  # Validate required args for model path
  if (is.null(treat)) stop("'treat' is required when model is provided")
  if (is.null(group)) stop("'group' is required when model is provided")
  if (!is.data.frame(data)) stop("'data' must be a data.frame")
  if (!treat %in% names(data)) stop("'", treat, "' not found in data")
  if (!group %in% names(data)) stop("'", group, "' not found in data")
  if (!is.null(time) && !time %in% names(data)) stop("'", time, "' not found in data")
  if (!is.null(post) && !post %in% names(data)) stop("'", post, "' not found in data")

  # Apply na.action to model-relevant columns. Mirrors the integrated ineqx()
  # path's .complete_model_data() so the manual two-step entry point
  # (fit_ineqx_model() + ineqx_params(model = ...)) handles NAs identically.
  vars <- unique(c(
    all.vars(formula(model, what = "mu")),
    all.vars(formula(model, what = "sigma")),
    treat, group, time, post
  ))
  vars <- intersect(vars, names(data))
  n_before <- nrow(data)
  data <- na.action(data[, vars, drop = FALSE])
  n_after <- nrow(data)
  if (verbose && n_after < n_before) {
    message("  N = ", n_after, " (", n_before - n_after,
            " rows omitted by na.action).")
  }

  # vcov defaults to TRUE for model path
  if (is.null(vcov)) vcov <- TRUE
  if (!is.logical(vcov)) {
    stop("'vcov' must be TRUE or FALSE when model is provided")
  }

  treat_vals <- sort(unique(data[[treat]]))
  if (!(0 %in% treat_vals)) {
    stop("Treatment variable '", treat, "' must contain 0 (untreated)")
  }
  treat_levels <- treat_vals[treat_vals != 0]

  if (length(treat_levels) > 1) {
    message("Multiple treatment levels detected: ",
            paste(treat_levels, collapse = ", "),
            ". Effect computed as weighted average.")
  }

  # Create time variable if missing
  if (is.null(time)) {
    data$.time <- 1L
    time_var <- ".time"
  } else {
    time_var <- time
  }

  group_levels <- sort(unique(data[[group]]))
  time_levels <- sort(unique(data[[time_var]]))

  is_did <- !is.null(post)

  # -------------------------------------------------------------------- #
  # Generate predictions under counterfactual treatment (and post, for DiD)
  #
  # Simple difference: 2 predictions per row (D=0, D=1).
  # DiD:               4 predictions per row (D x P combinations).
  # -------------------------------------------------------------------- #
  if (verbose) message("  Generating counterfactual predictions...")

  .pred_at <- function(treat_val, post_val = NULL) {
    nd <- data
    nd[[treat]] <- treat_val
    if (!is.null(post_val)) nd[[post]] <- post_val
    list(
      mu    = predict(model, what = "mu",    type = "response",
                      newdata = nd, data = data),
      sigma = predict(model, what = "sigma", type = "response",
                      newdata = nd, data = data)
    )
  }

  .pred_at_treat1 <- function(post_val = NULL) {
    if (length(treat_levels) == 1) {
      return(.pred_at(treat_levels[1], post_val))
    }
    mu_acc    <- rep(0, nrow(data))
    sigma_acc <- rep(0, nrow(data))
    for (tl in treat_levels) {
      p <- .pred_at(tl, post_val)
      freq <- sum(data[[treat]] == tl) / sum(data[[treat]] != 0)
      mu_acc    <- mu_acc    + freq * p$mu
      sigma_acc <- sigma_acc + freq * p$sigma
    }
    list(mu = mu_acc, sigma = sigma_acc)
  }

  if (is_did) {
    p_D0P0 <- .pred_at(0,        0)
    p_D0P1 <- .pred_at(0,        1)
    p_D1P0 <- .pred_at_treat1(0)
    p_D1P1 <- .pred_at_treat1(1)
  } else {
    p_D0 <- .pred_at(0)
    p_D1 <- .pred_at_treat1()
  }

  # -------------------------------------------------------------------- #
  # Per-observation counterfactual potential outcomes m_hat_i(d), s_hat_i(d)
  # (used by the marginal estimand). For DiD these are constructed at the
  # individual level via the parallel-trends identity (paper eqs 45-46)
  # BEFORE any aggregation, which the marginal law-of-total-variance
  # aggregation (paper eqs 10, 37-38) requires.
  # -------------------------------------------------------------------- #
  if (estimand == "marginal") {
    if (is_did) {
      po_m0 <- p_D1P0$mu + (p_D0P1$mu - p_D0P0$mu)
      po_m1 <- p_D1P1$mu
      po_s0 <- p_D1P0$sigma * p_D0P1$sigma / p_D0P0$sigma
      po_s1 <- p_D1P1$sigma
      # Pre-period anchors (P=0): treated (D=1) vs untreated (D=0) outcomes,
      # for the Figure-3 / pretrends panels.
      pre_m1 <- p_D1P0$mu; pre_s1 <- p_D1P0$sigma
      pre_m0 <- p_D0P0$mu; pre_s0 <- p_D0P0$sigma
    } else {
      po_m0 <- p_D0$mu; po_m1 <- p_D1$mu
      po_s0 <- p_D0$sigma; po_s1 <- p_D1$sigma
    }
  }

  # -------------------------------------------------------------------- #
  # Compute per group-time quantities
  #
  # For DiD, the canonical estimands are computed on the treated post-period
  # subpopulation (D=1, P=1) within each (group, time) cell. The
  # counterfactual untreated mean for that subpopulation is constructed via
  # the parallel-trends identity:
  #   mu0 = E[mu_hat(D=1, P=0)] + (E[mu_hat(D=0, P=1)] - E[mu_hat(D=0, P=0)])
  # so that beta = mu1 - mu0 = beta_DP (the ATT). Lambda is analogous on
  # the log-SD scale.
  # -------------------------------------------------------------------- #

  result_list <- list()
  for (t in time_levels) {
    for (g in group_levels) {
      idx_cell <- data[[time_var]] == t & data[[group]] == g
      if (sum(idx_cell) == 0) next

      # Target sub-population: treated post-period (DiD) or whole cell (simple).
      if (is_did) {
        idx_tp <- idx_cell & data[[treat]] != 0 & data[[post]] == 1
        if (sum(idx_tp) == 0) {
          # No treated post-period observations in this cell — fall back to
          # whole-cell averaging so the row is still produced.
          idx_tp <- idx_cell
        }
        idx_use <- which(idx_tp)
      } else {
        idx_use <- which(idx_cell)
      }
      n_treated <- sum(data[[treat]][idx_use] != 0)

      if (estimand == "marginal") {
        # ---- Marginal: law of total variance over the covariate dist. ----
        gm <- .marginal_group_moments(po_m0, po_m1, po_s0, po_s1, idx_use)
        mu0_g    <- gm$mu0;    mu1_g    <- gm$mu1
        sigma0_g <- gm$sigma0; sigma1_g <- gm$sigma1
        beta_g   <- gm$beta;   lambda_g <- gm$lambda

        if (is_did) {
          an <- .marginal_group_moments(pre_m0, pre_m1, pre_s0, pre_s1, idx_use)
          mu0_pre_g    <- an$mu0;    mu1_pre_g    <- an$mu1
          sigma0_pre_g <- an$sigma0; sigma1_pre_g <- an$sigma1
        } else {
          mu1_pre_g    <- NA_real_
          mu0_pre_g    <- NA_real_
          sigma1_pre_g <- NA_real_
          sigma0_pre_g <- NA_real_
        }
      } else if (is_did) {
        # ---- Residual (conditional) variant, DiD ----
        mu_D0P0 <- mean(p_D0P0$mu[idx_use])
        mu_D0P1 <- mean(p_D0P1$mu[idx_use])
        mu_D1P0 <- mean(p_D1P0$mu[idx_use])
        mu_D1P1 <- mean(p_D1P1$mu[idx_use])

        lsig_D0P0 <- mean(log(p_D0P0$sigma[idx_use]))
        lsig_D0P1 <- mean(log(p_D0P1$sigma[idx_use]))
        lsig_D1P0 <- mean(log(p_D1P0$sigma[idx_use]))
        lsig_D1P1 <- mean(log(p_D1P1$sigma[idx_use]))

        mu1_g    <- mu_D1P1
        mu0_g    <- mu_D1P0 + (mu_D0P1 - mu_D0P0)  # DiD-implied counterfactual
        beta_g   <- mu1_g - mu0_g                  # = beta_DP (ATT)

        lsig1     <- lsig_D1P1
        lsig0_DiD <- lsig_D1P0 + (lsig_D0P1 - lsig_D0P0)
        sigma1_g  <- exp(lsig1)
        sigma0_g  <- exp(lsig0_DiD)
        lambda_g  <- lsig1 - lsig0_DiD             # = lambda_DP

        # Pre-period anchors (used by outcome.params and pretrends plots)
        mu1_pre_g    <- mu_D1P0
        mu0_pre_g    <- mu_D0P0
        sigma1_pre_g <- exp(lsig_D1P0)
        sigma0_pre_g <- exp(lsig_D0P0)
      } else {
        # ---- Residual (conditional) variant, simple difference ----
        mu0_g     <- mean(p_D0$mu[idx_use])
        sigma0_g  <- mean(p_D0$sigma[idx_use])
        beta_g    <- mean(p_D1$mu[idx_use] - p_D0$mu[idx_use])
        lambda_g  <- mean(log(p_D1$sigma[idx_use]) -
                          log(p_D0$sigma[idx_use]))

        mu1_g    <- mu0_g + beta_g
        sigma1_g <- sigma0_g * exp(lambda_g)

        mu1_pre_g    <- NA_real_
        mu0_pre_g    <- NA_real_
        sigma1_pre_g <- NA_real_
        sigma0_pre_g <- NA_real_
      }

      result_list[[length(result_list) + 1]] <- data.frame(
        group = g, time = t, n_treated = n_treated,
        mu0 = mu0_g, sigma0 = sigma0_g,
        mu1 = mu1_g, sigma1 = sigma1_g,
        beta = beta_g, lambda = lambda_g,
        mu0_pre = mu0_pre_g, mu1_pre = mu1_pre_g,
        sigma0_pre = sigma0_pre_g, sigma1_pre = sigma1_pre_g,
        stringsAsFactors = FALSE
      )
    }
  }

  result_df <- do.call(rbind, result_list)

  # Normalize pi within each time period (share of treated [post-period, in
  # the DiD case] subpopulation in each group)
  result_df$pi <- NA_real_
  for (t in time_levels) {
    idx <- result_df$time == t
    total_treated <- sum(result_df$n_treated[idx])
    if (total_treated > 0) {
      result_df$pi[idx] <- result_df$n_treated[idx] / total_treated
    } else {
      result_df$pi[idx] <- 1 / sum(idx)
    }
  }

  base_cols <- c("group", "time", "pi",
                 "mu0", "sigma0", "mu1", "sigma1", "beta", "lambda")
  did_cols  <- c("mu0_pre", "mu1_pre", "sigma0_pre", "sigma1_pre")
  out_cols  <- if (is_did) c(base_cols, did_cols) else base_cols
  out_df    <- result_df[, out_cols]

  # Extract vcov if requested
  vcov_mat <- NULL
  if (vcov) {
    message("  Computing variance-covariance matrix...")
    vcov_mat <- .extract_vcov_numerical(
      model = model, treat = treat, group = group,
      time_var = time_var, data = data,
      group_levels = group_levels, time_levels = time_levels,
      out_df = out_df, estimand = estimand,
      post = if (is_did) post else NULL
    )
  }

  # Remove time column if cross-sectional
  if (is.null(time)) {
    out_df$time <- NULL
  }

  # Use manual path for validation and S3 construction, then attach DiD flags
  obj <- .ineqx_params_manual(data = out_df, ystat = ystat, vcov = vcov_mat,
                              estimand = estimand)
  obj$is_did <- is_did
  obj$post   <- if (is_did) post else NULL
  # Propagate the response-scale tag set by fit_ineqx_model(). "log" tells
  # downstream ineqx(params, ystat = "VL") that it can treat this fit as a
  # variance of log(y) decomposition; "identity" or NULL means VL is not
  # available through the split-step path.
  obj$transform <- attr(model, "ineqx_transform")
  if (is.null(obj$transform)) obj$transform <- "identity"
  obj
}


# ---------------------------------------------------------------------------- #
# Internal: Construct the S3 object
# ---------------------------------------------------------------------------- #

.make_ineqx_params <- function(data, ystat, type, vcov, has_time,
                               estimand = "marginal") {
  structure(
    list(
      data = data,
      ystat = ystat,
      estimand = estimand,
      type = type,
      vcov = vcov,
      groups = sort(unique(data$group)),
      times = if (has_time) sort(unique(data$time)) else NULL,
      n_groups = length(unique(data$group)),
      n_times = if (has_time) length(unique(data$time)) else 1L
    ),
    class = "ineqx_params"
  )
}


# ---------------------------------------------------------------------------- #
# Internal: Marginal group moments via the law of total variance
#
# Aggregates per-observation counterfactual predictions into the marginal
# group-level moments used by the decomposition (paper eqs 10, 37-38):
#   mu_j(d)  = mean_i m_hat_i(d)
#   S2_j(d)  = mean_i [ s_hat_i(d)^2 + m_hat_i(d)^2 ] - mu_j(d)^2
# i.e. the average conditional residual variance PLUS the dispersion of
# predicted conditional means across covariate profiles within the group.
#
#   m0, m1 : predicted conditional means  m_hat_i(0), m_hat_i(1)  (full vectors)
#   s0, s1 : predicted conditional SDs    s_hat_i(0), s_hat_i(1)  (full vectors)
#   idx    : integer row indices of the target sub-population in the cell
# Returns list(mu0, mu1, sigma0, sigma1, beta, lambda).
# ---------------------------------------------------------------------------- #

.marginal_group_moments <- function(m0, m1, s0, s1, idx) {
  m0 <- m0[idx]; m1 <- m1[idx]; s0 <- s0[idx]; s1 <- s1[idx]
  mu0  <- mean(m0)
  mu1  <- mean(m1)
  S2_0 <- mean(s0^2 + m0^2) - mu0^2
  S2_1 <- mean(s1^2 + m1^2) - mu1^2
  list(
    mu0    = mu0,
    mu1    = mu1,
    sigma0 = sqrt(S2_0),
    sigma1 = sqrt(S2_1),
    beta   = mu1 - mu0,
    lambda = 0.5 * log(S2_1 / S2_0)
  )
}


# ---------------------------------------------------------------------------- #
# Internal: Extract raw coefficient vcov from QR decompositions
# ---------------------------------------------------------------------------- #

#' Extract block-diagonal coefficient vcov from a gamlss model
#'
#' Bypasses \code{gamlss::vcov.gamlss()} / \code{gen.likelihood()} which use
#' \code{get()} to resolve \code{data} from the model's call — this fails when
#' the model was created inside \code{fit_ineqx_model()} because the call
#' references local variable names that no longer exist.
#'
#' Instead, we extract the per-equation vcov directly from the QR
#' decompositions stored on the model (\code{model$mu.qr},
#' \code{model$sigma.qr}) and combine them into a block-diagonal matrix.
#'
#' @param model A fitted gamlss object
#' @return A (p_mu + p_sigma) x (p_mu + p_sigma) block-diagonal vcov matrix
#' @keywords internal
.extract_raw_vcov <- function(model) {
  mu_coefs <- coef(model, what = "mu")
  sigma_coefs <- coef(model, what = "sigma")
  p_mu <- length(mu_coefs)
  p_sigma <- length(sigma_coefs)

  # Extract per-equation vcov from QR decompositions
  R_mu <- qr.R(model$mu.qr)[1:p_mu, 1:p_mu, drop = FALSE]
  vcov_mu <- chol2inv(R_mu)

  R_sigma <- qr.R(model$sigma.qr)[1:p_sigma, 1:p_sigma, drop = FALSE]
  vcov_sigma <- chol2inv(R_sigma)

  # Combine into block-diagonal matrix
  p <- p_mu + p_sigma
  raw_vcov <- matrix(0, p, p)
  raw_vcov[1:p_mu, 1:p_mu] <- vcov_mu
  raw_vcov[(p_mu + 1):p, (p_mu + 1):p] <- vcov_sigma

  # Transfer names
  nms <- c(names(mu_coefs), names(sigma_coefs))
  dimnames(raw_vcov) <- list(nms, nms)

  raw_vcov
}


# ---------------------------------------------------------------------------- #
# Internal: Extract vcov via numerical Jacobian
# ---------------------------------------------------------------------------- #

#' @keywords internal
.extract_vcov_numerical <- function(model, treat, group, time_var, data,
                                     group_levels, time_levels, out_df,
                                     estimand = "marginal", post = NULL) {

  # Extract the raw coefficient vcov from the model's QR decompositions.
  raw_vcov <- tryCatch(
    .extract_raw_vcov(model),
    error = function(e) {
      warning("Could not extract vcov from model: ", e$message,
              ". Returning NULL.")
      return(NULL)
    }
  )
  if (is.null(raw_vcov)) return(NULL)

  mu_coefs <- coef(model, what = "mu")
  sigma_coefs <- coef(model, what = "sigma")
  raw_coefs <- c(mu_coefs, sigma_coefs)
  p <- length(raw_coefs)
  p_mu <- length(mu_coefs)
  p_sigma <- length(sigma_coefs)

  if (nrow(raw_vcov) != p || ncol(raw_vcov) != p) {
    vcov_mu <- tryCatch(stats::vcov(model, what = "mu"), error = function(e) NULL)
    vcov_sigma <- tryCatch(stats::vcov(model, what = "sigma"), error = function(e) NULL)
    if (!is.null(vcov_mu) && !is.null(vcov_sigma)) {
      raw_vcov <- matrix(0, p, p)
      raw_vcov[1:p_mu, 1:p_mu] <- vcov_mu
      raw_vcov[(p_mu + 1):p, (p_mu + 1):p] <- vcov_sigma
    } else {
      warning("Could not construct full vcov matrix. Returning NULL.")
      return(NULL)
    }
  }

  # --- Fast analytical Jacobian via model matrices ---
  # Build model matrices for counterfactual predictions (once, not per coef).
  # For simple-difference: 2 matrices each (D=0, D=1).
  # For DiD: 4 matrices each (D x P combinations) — matches the corrected
  # extraction in .ineqx_params_from_model.
  is_did <- !is.null(post)

  mu_form <- formula(model, what = "mu")
  sigma_form <- formula(model, what = "sigma")
  if (length(mu_form) == 3) mu_form[2] <- NULL
  if (length(sigma_form) == 3) sigma_form[2] <- NULL

  .mm <- function(form, dat, what) {
    if (what == "mu") {
      model.matrix(form, data = dat, xlev = model$mu.xlevels,
                   contrasts.arg = model$contrasts)
    } else {
      model.matrix(form, data = dat, xlev = model$sigma.xlevels,
                   contrasts.arg = model$contrasts)
    }
  }

  if (is_did) {
    .with_dp <- function(d_val, p_val) {
      dat <- data; dat[[treat]] <- d_val; dat[[post]] <- p_val; dat
    }
    X_mu_D0P0  <- .mm(mu_form,    .with_dp(0, 0), "mu")
    X_mu_D0P1  <- .mm(mu_form,    .with_dp(0, 1), "mu")
    X_mu_D1P0  <- .mm(mu_form,    .with_dp(1, 0), "mu")
    X_mu_D1P1  <- .mm(mu_form,    .with_dp(1, 1), "mu")
    X_sig_D0P0 <- .mm(sigma_form, .with_dp(0, 0), "sigma")
    X_sig_D0P1 <- .mm(sigma_form, .with_dp(0, 1), "sigma")
    X_sig_D1P0 <- .mm(sigma_form, .with_dp(1, 0), "sigma")
    X_sig_D1P1 <- .mm(sigma_form, .with_dp(1, 1), "sigma")

    # DiD-implied linear combinations:
    #   beta_j(t)   = mean_idx_tp(X_mu_D1P1 - X_mu_D1P0 - X_mu_D0P1 + X_mu_D0P0) %*% mu_coefs
    #   mu0_j(t)    = mean_idx_tp(X_mu_D1P0 + X_mu_D0P1 - X_mu_D0P0) %*% mu_coefs
    #   lambda_j(t) = mean_idx_tp(X_sig_D1P1 - X_sig_D1P0 - X_sig_D0P1 + X_sig_D0P0) %*% sigma_coefs
    #   log_sigma0_j(t) = mean_idx_tp(X_sig_D1P0 + X_sig_D0P1 - X_sig_D0P0) %*% sigma_coefs
    X_mu_beta  <- X_mu_D1P1  - X_mu_D1P0  - X_mu_D0P1  + X_mu_D0P0
    X_mu_mu0   <- X_mu_D1P0  + X_mu_D0P1  - X_mu_D0P0
    X_sig_lam  <- X_sig_D1P1 - X_sig_D1P0 - X_sig_D0P1 + X_sig_D0P0
    X_sig_lsig <- X_sig_D1P0 + X_sig_D0P1 - X_sig_D0P0
  } else {
    data_cf0 <- data; data_cf0[[treat]] <- 0
    data_cf1 <- data; data_cf1[[treat]] <- 1
    X_mu_0  <- .mm(mu_form,    data_cf0, "mu")
    X_mu_1  <- .mm(mu_form,    data_cf1, "mu")
    X_sig_0 <- .mm(sigma_form, data_cf0, "sigma")
    X_sig_1 <- .mm(sigma_form, data_cf1, "sigma")

    X_mu_beta  <- X_mu_1 - X_mu_0
    X_mu_mu0   <- X_mu_0
    X_sig_lam  <- X_sig_1 - X_sig_0
    X_sig_lsig <- X_sig_0
  }

  # Theta vector per time-group, ordered (within each time block) as
  # [beta_1..J, mu0_1..J, lambda_1..J, log_sigma0_1..J] to match the layout
  # consumed by delta_method_se().

  J <- length(group_levels)
  n_t <- length(time_levels)
  n_theta <- 4 * J * n_t

  # Precompute the per-(time, group) target index set: treated post-period
  # for DiD, whole cell otherwise. This is the identical averaging set used
  # by the point estimate in .ineqx_params_from_model().
  cell_idx <- vector("list", n_t * J)
  for (ti in seq_along(time_levels)) {
    for (ji in seq_along(group_levels)) {
      idx_cell <- which(data[[time_var]] == time_levels[ti] &
                          data[[group]] == group_levels[ji])
      if (is_did && length(idx_cell) > 0) {
        idx_tp <- idx_cell[data[[treat]][idx_cell] != 0 &
                            data[[post]][idx_cell] == 1]
        if (length(idx_tp) > 0) idx_cell <- idx_tp
      }
      cell_idx[[(ti - 1) * J + ji]] <- idx_cell
    }
  }

  if (estimand == "residual") {
    # --- Fast analytical Jacobian: theta is linear in the raw coefficients ---
    jacobian <- matrix(0, n_theta, p)
    for (ti in seq_along(time_levels)) {
      for (ji in seq_along(group_levels)) {
        idx <- cell_idx[[(ti - 1) * J + ji]]
        if (length(idx) == 0) next
        offset <- (ti - 1) * 4 * J
        jacobian[offset + ji,         1:p_mu]       <- colMeans(X_mu_beta[idx, , drop = FALSE])
        jacobian[offset + J + ji,     1:p_mu]       <- colMeans(X_mu_mu0[idx, , drop = FALSE])
        jacobian[offset + 2 * J + ji, (p_mu + 1):p] <- colMeans(X_sig_lam[idx, , drop = FALSE])
        jacobian[offset + 3 * J + ji, (p_mu + 1):p] <- colMeans(X_sig_lsig[idx, , drop = FALSE])
      }
    }
    full_vcov <- jacobian %*% raw_vcov %*% t(jacobian)
  } else {
    # --- Numerical Jacobian of the full marginal g-computation routine ---
    # The marginal S0/Lambda are nonlinear in BOTH the mu- and sigma-equation
    # coefficients (the law of total variance folds the predicted conditional
    # means into the group variance), so the Jacobian gamma -> theta is taken by
    # forward finite differences over the prediction + g-computation map (paper
    # Appendix C/H). Predictions are reconstructed from the prebuilt model
    # matrices (no predict() calls), and -- because perturbing coefficient k
    # shifts only one equation's linear predictor by h * X[, k] -- the loop
    # updates just the affected side incrementally (a rank-1 update of the base
    # predictions) rather than recomputing all matvecs. This makes the FD loop
    # O(N * p) instead of O(N * p^2).
    mu_linkinv  <- stats::make.link(model$mu.link)$linkinv
    sig_linkinv <- stats::make.link(model$sigma.link)$linkinv

    # Per-state design matrices, keyed so a perturbation can update one state set.
    if (is_did) {
      Xmu  <- list(D0P0 = X_mu_D0P0,  D0P1 = X_mu_D0P1,
                   D1P0 = X_mu_D1P0,  D1P1 = X_mu_D1P1)
      Xsig <- list(D0P0 = X_sig_D0P0, D0P1 = X_sig_D0P1,
                   D1P0 = X_sig_D1P0, D1P1 = X_sig_D1P1)
    } else {
      Xmu  <- list(`0` = X_mu_0,  `1` = X_mu_1)
      Xsig <- list(`0` = X_sig_0, `1` = X_sig_1)
    }

    mu_c0  <- raw_coefs[1:p_mu]
    sig_c0 <- raw_coefs[(p_mu + 1):p]

    # Base linear predictors per state (the only matvecs; computed once).
    eta_mu0  <- lapply(Xmu,  function(X) as.vector(X %*% mu_c0))
    eta_sig0 <- lapply(Xsig, function(X) as.vector(X %*% sig_c0))

    # theta from per-state mean/sd PREDICTION lists. Builds the potential
    # outcomes (same po formulas as the point estimate) and aggregates via the
    # shared .marginal_group_moments(), producing the 4J * n_t theta vector in
    # the [beta, mu0, lambda, log_sigma0] layout per time block.
    .theta_from_pred <- function(m, s) {
      if (is_did) {
        po_m0 <- m$D1P0 + (m$D0P1 - m$D0P0); po_m1 <- m$D1P1
        po_s0 <- s$D1P0 * s$D0P1 / s$D0P0;   po_s1 <- s$D1P1
      } else {
        po_m0 <- m[["0"]]; po_m1 <- m[["1"]]
        po_s0 <- s[["0"]]; po_s1 <- s[["1"]]
      }
      out <- numeric(n_theta)
      for (ti in seq_along(time_levels)) {
        offset <- (ti - 1) * 4 * J
        for (ji in seq_along(group_levels)) {
          idx <- cell_idx[[(ti - 1) * J + ji]]
          if (length(idx) == 0) next
          gm <- .marginal_group_moments(po_m0, po_m1, po_s0, po_s1, idx)
          out[offset + ji]         <- gm$beta
          out[offset + J + ji]     <- gm$mu0
          out[offset + 2 * J + ji] <- gm$lambda
          out[offset + 3 * J + ji] <- log(gm$sigma0)
        }
      }
      out
    }

    m_pred0 <- lapply(eta_mu0,  mu_linkinv)
    s_pred0 <- lapply(eta_sig0, sig_linkinv)
    theta0  <- .theta_from_pred(m_pred0, s_pred0)

    jacobian <- matrix(0, n_theta, p)
    h_vec <- pmax(abs(raw_coefs) * 1e-6, 1e-7)
    for (k in seq_len(p)) {
      h <- h_vec[k]
      if (k <= p_mu) {
        # Mean coefficient: only the mu predictions shift; reuse base s_pred0.
        m_p   <- Map(function(eta, X) mu_linkinv(eta + h * X[, k]), eta_mu0, Xmu)
        th    <- .theta_from_pred(m_p, s_pred0)
      } else {
        # Scale coefficient: only the sigma predictions shift; reuse base m_pred0.
        j     <- k - p_mu
        s_p   <- Map(function(eta, X) sig_linkinv(eta + h * X[, j]), eta_sig0, Xsig)
        th    <- .theta_from_pred(m_pred0, s_p)
      }
      jacobian[, k] <- (th - theta0) / h
    }
    full_vcov <- jacobian %*% raw_vcov %*% t(jacobian)
  }

  if (n_t == 1) {
    return(full_vcov)
  }

  vcov_list <- list()
  for (i in seq_along(time_levels)) {
    idx <- ((i - 1) * 4 * J + 1):(i * 4 * J)
    vcov_list[[as.character(time_levels[i])]] <- full_vcov[idx, idx, drop = FALSE]
  }
  vcov_list
}


# ---------------------------------------------------------------------------- #
# Print method
# ---------------------------------------------------------------------------- #

#' @export
print.ineqx_params <- function(x, ...) {
  cat("ineqx_params object\n")
  cat("  Type:", x$type, "\n")
  cat("  Inequality measure:", x$ystat, "\n")
  cat("  Estimand:", x$estimand %||% "marginal", "\n")
  cat("  Groups (", x$n_groups, "): ",
      paste(x$groups, collapse = ", "), "\n", sep = "")
  if (!is.null(x$times)) {
    cat("  Time periods (", x$n_times, "): ",
        paste(x$times, collapse = ", "), "\n", sep = "")
  }

  d <- x$data
  show_cols <- c("group", "pi", "mu0", "sigma0", "beta", "lambda")

  if (!is.null(x$times)) {
    cat("\nGroup-level parameters:\n")
    for (t in x$times) {
      cat(sprintf("\n  time %s:\n", t))
      sub <- d[d$time == t, show_cols, drop = FALSE]
      rownames(sub) <- NULL
      .print_indented_df(sub, indent = 4, digits = 4)
    }
  } else {
    cat("\nGroup-level parameters:\n")
    sub <- d[, show_cols, drop = FALSE]
    rownames(sub) <- NULL
    .print_indented_df(sub, indent = 2, digits = 4)
  }

  if (!is.null(x$vcov)) {
    if (is.matrix(x$vcov)) {
      cat(sprintf("\nvcov: %dx%d\n", nrow(x$vcov), ncol(x$vcov)))
    } else {
      cat(sprintf("\nvcov: %d time periods\n", length(x$vcov)))
    }
  }
  invisible(x)
}


# ---------------------------------------------------------------------------- #
# Internal: Descriptive params manual specification
# ---------------------------------------------------------------------------- #

.ineqx_desc_params_manual <- function(data, ystat) {

  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame")
  }

  # Check required columns
  required_cols <- c("group", "pi", "mu", "sigma")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in 'data': ",
         paste(missing_cols, collapse = ", "),
         "\nFor descriptive params, required columns: group, pi, mu, sigma",
         "\nFor causal params, required columns: group, pi, mu0, sigma0, beta, lambda")
  }

  # Validate numeric columns
  numeric_cols <- c("pi", "mu", "sigma")
  for (col in numeric_cols) {
    if (!is.numeric(data[[col]])) {
      stop("Column '", col, "' must be numeric")
    }
    if (any(is.na(data[[col]]))) {
      stop("Column '", col, "' contains NA values")
    }
  }

  # Validate sigma >= 0
  if (any(data$sigma < 0)) {
    stop("'sigma' must be non-negative")
  }

  # Determine if longitudinal (optional time column)
  has_time <- "time" %in% names(data)

  # Validate pi sums to 1 (within each time period if longitudinal)
  if (has_time) {
    times <- unique(data$time)
    for (t in times) {
      pi_sum <- sum(data$pi[data$time == t])
      if (abs(pi_sum - 1) > 1e-6) {
        stop("'pi' must sum to 1 within each time period. ",
             "At time = ", t, ", sum(pi) = ", round(pi_sum, 6))
      }
    }
  } else {
    pi_sum <- sum(data$pi)
    if (abs(pi_sum - 1) > 1e-6) {
      stop("'pi' must sum to 1. Currently sum(pi) = ", round(pi_sum, 6))
    }
  }

  # Validate no duplicate group(-time) combinations
  if (has_time) {
    gt <- paste(data$group, data$time, sep = "_")
  } else {
    gt <- as.character(data$group)
  }
  if (any(duplicated(gt))) {
    stop("Duplicate group", if (has_time) "-time", " combinations found in 'data'")
  }

  # Sort for consistency
  if (has_time) {
    data <- data[order(data$time, data$group), ]
  } else {
    data <- data[order(data$group), ]
  }

  structure(
    list(
      data = data,
      ystat = ystat,
      groups = sort(unique(data$group)),
      n_groups = length(unique(data$group)),
      times = if (has_time) sort(unique(data$time)) else NULL,
      n_times = if (has_time) length(unique(data$time)) else NULL
    ),
    class = "ineqx_desc_params"
  )
}


# ---------------------------------------------------------------------------- #
# Print method for descriptive params
# ---------------------------------------------------------------------------- #

#' @export
print.ineqx_desc_params <- function(x, ...) {
  cat("ineqx_desc_params object (descriptive counterfactual reference)\n")
  cat("  Inequality measure:", x$ystat, "\n")
  cat("  Groups (", x$n_groups, "): ",
      paste(x$groups, collapse = ", "), "\n", sep = "")
  if (!is.null(x$times)) {
    cat("  Time periods (", x$n_times, "): ",
        paste(x$times, collapse = ", "), "\n", sep = "")
  }

  d <- x$data
  show_cols <- c("group", "pi", "mu", "sigma")

  if (!is.null(x$times)) {
    cat("\nGroup-level parameters:\n")
    for (t in x$times) {
      cat(sprintf("\n  time %s:\n", t))
      sub <- d[d$time == t, show_cols, drop = FALSE]
      rownames(sub) <- NULL
      .print_indented_df(sub, indent = 4, digits = 4)
    }
  } else {
    cat("\nGroup-level parameters:\n")
    sub <- d[, show_cols, drop = FALSE]
    rownames(sub) <- NULL
    .print_indented_df(sub, indent = 2, digits = 4)
  }
  invisible(x)
}
