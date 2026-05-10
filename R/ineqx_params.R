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
#'   Only used in model mode. NULL for simple difference estimator.
#' @param ystat Character, either \code{"Var"} (variance) or \code{"CV2"}
#'   (squared coefficient of variation). Default: \code{"Var"}.
#' @param vcov For manual mode: an optional variance-covariance matrix (or
#'   named list of matrices for longitudinal data).
#'   For model mode: logical, whether to extract the vcov via numerical
#'   Jacobian (default: TRUE). Set to FALSE for faster computation without SEs.
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
#' }
#'
#' @export
ineqx_params <- function(data, model = NULL, treat = NULL, group = NULL,
                          time = NULL, post = NULL,
                          ystat = "Var", vcov = NULL, verbose = TRUE) {

  ystat <- match.arg(ystat, c("Var", "CV2"))

  if (!is.null(model)) {
    # -------------------------------------------------------------------- #
    # MODEL EXTRACTION PATH
    # -------------------------------------------------------------------- #
    .ineqx_params_from_model(
      model = model, data = data,
      treat = treat, group = group,
      time = time, post = post,
      ystat = ystat, vcov = vcov
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
      .ineqx_params_manual(data = data, ystat = ystat, vcov = vcov)
    }
  }
}


# ---------------------------------------------------------------------------- #
# Internal: Manual specification path
# ---------------------------------------------------------------------------- #

.ineqx_params_manual <- function(data, ystat, vcov) {

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

  .make_ineqx_params(data, ystat, type, vcov, has_time)
}


# ---------------------------------------------------------------------------- #
# Internal: Model extraction path
# ---------------------------------------------------------------------------- #

.ineqx_params_from_model <- function(model, data, treat, group,
                                      time, post, ystat, vcov, verbose = TRUE) {

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

      if (is_did) {
        idx_tp <- idx_cell & data[[treat]] != 0 & data[[post]] == 1
        if (sum(idx_tp) == 0) {
          # No treated post-period observations in this cell — fall back to
          # whole-cell averaging so the row is still produced (with a warning
          # surfaced via NA n_treated_post).
          idx_tp <- idx_cell
        }
        mu_D0P0 <- mean(p_D0P0$mu[idx_tp])
        mu_D0P1 <- mean(p_D0P1$mu[idx_tp])
        mu_D1P0 <- mean(p_D1P0$mu[idx_tp])
        mu_D1P1 <- mean(p_D1P1$mu[idx_tp])

        lsig_D0P0 <- mean(log(p_D0P0$sigma[idx_tp]))
        lsig_D0P1 <- mean(log(p_D0P1$sigma[idx_tp]))
        lsig_D1P0 <- mean(log(p_D1P0$sigma[idx_tp]))
        lsig_D1P1 <- mean(log(p_D1P1$sigma[idx_tp]))

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

        n_treated <- sum(data[[treat]][idx_tp] != 0)
      } else {
        mu0_g     <- mean(p_D0$mu[idx_cell])
        sigma0_g  <- mean(p_D0$sigma[idx_cell])
        beta_g    <- mean(p_D1$mu[idx_cell] - p_D0$mu[idx_cell])
        lambda_g  <- mean(log(p_D1$sigma[idx_cell]) -
                          log(p_D0$sigma[idx_cell]))

        mu1_g    <- mu0_g + beta_g
        sigma1_g <- sigma0_g * exp(lambda_g)

        mu1_pre_g    <- NA_real_
        mu0_pre_g    <- NA_real_
        sigma1_pre_g <- NA_real_
        sigma0_pre_g <- NA_real_

        n_treated <- sum(data[[treat]][idx_cell] != 0)
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
      out_df = out_df,
      post = if (is_did) post else NULL
    )
  }

  # Remove time column if cross-sectional
  if (is.null(time)) {
    out_df$time <- NULL
  }

  # Use manual path for validation and S3 construction, then attach DiD flags
  obj <- .ineqx_params_manual(data = out_df, ystat = ystat, vcov = vcov_mat)
  obj$is_did <- is_did
  obj$post   <- if (is_did) post else NULL
  obj
}


# ---------------------------------------------------------------------------- #
# Internal: Construct the S3 object
# ---------------------------------------------------------------------------- #

.make_ineqx_params <- function(data, ystat, type, vcov, has_time) {
  structure(
    list(
      data = data,
      ystat = ystat,
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
                                     post = NULL) {

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

  # Theta vector per time-group: [beta_1..J, mu0_1..J, lambda_1..J, log_sigma0_1..J]
  # All quantities are linear in the raw coefficients, so the Jacobian is analytic.
  # The averaging set differs by mode: simple-difference averages over the cell;
  # DiD averages over the treated post-period subset of the cell (idx_tp).

  J <- length(group_levels)
  n_t <- length(time_levels)
  n_theta <- 4 * J * n_t
  jacobian <- matrix(0, n_theta, p)

  for (ti in seq_along(time_levels)) {
    t_val <- time_levels[ti]
    for (ji in seq_along(group_levels)) {
      g <- group_levels[ji]
      idx_cell <- which(data[[time_var]] == t_val & data[[group]] == g)
      if (length(idx_cell) == 0) next

      if (is_did) {
        idx <- idx_cell[data[[treat]][idx_cell] != 0 &
                         data[[post]][idx_cell] == 1]
        if (length(idx) == 0) idx <- idx_cell
      } else {
        idx <- idx_cell
      }

      offset <- (ti - 1) * 4 * J
      jacobian[offset + ji,           1:p_mu]      <- colMeans(X_mu_beta[idx, , drop = FALSE])
      jacobian[offset + J + ji,       1:p_mu]      <- colMeans(X_mu_mu0[idx, , drop = FALSE])
      jacobian[offset + 2 * J + ji,   (p_mu + 1):p] <- colMeans(X_sig_lam[idx, , drop = FALSE])
      jacobian[offset + 3 * J + ji,   (p_mu + 1):p] <- colMeans(X_sig_lsig[idx, , drop = FALSE])
    }
  }

  full_vcov <- jacobian %*% raw_vcov %*% t(jacobian)

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
