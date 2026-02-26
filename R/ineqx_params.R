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
#' There are two usage modes:
#' \describe{
#'   \item{Manual specification}{Pass a data.frame with columns
#'     \code{group}, \code{pi}, \code{mu0}, \code{sigma0}, \code{beta},
#'     \code{lambda} (and optionally \code{time}).}
#'   \item{From a fitted gamlss model}{Pass a fitted \code{gamlss} object
#'     via \code{model}, along with the raw data and variable names. The
#'     function extracts treatment effects via counterfactual predictions.}
#' }
#'
#' @param data For manual mode: a data.frame with required columns
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
#' @param ref Numeric, reference time period for longitudinal decomposition.
#'   Required when data contains multiple time periods.
#' @param vcov For manual mode: an optional variance-covariance matrix (or
#'   named list of matrices for longitudinal data).
#'   For model mode: logical, whether to extract the vcov via numerical
#'   Jacobian (default: TRUE). Set to FALSE for faster computation without SEs.
#'
#' @return An object of class \code{"ineqx_params"}
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
#' # Manual specification for a two-group example
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
                          ystat = "Var", ref = NULL, vcov = NULL) {

  ystat <- match.arg(ystat, c("Var", "CV2"))

  if (!is.null(model)) {
    # -------------------------------------------------------------------- #
    # MODEL EXTRACTION PATH
    # -------------------------------------------------------------------- #
    .ineqx_params_from_model(
      model = model, data = data,
      treat = treat, group = group,
      time = time, post = post,
      ystat = ystat, ref = ref, vcov = vcov
    )
  } else {
    # -------------------------------------------------------------------- #
    # MANUAL SPECIFICATION PATH
    # -------------------------------------------------------------------- #
    .ineqx_params_manual(data = data, ystat = ystat, ref = ref, vcov = vcov)
  }
}


# ---------------------------------------------------------------------------- #
# Internal: Manual specification path
# ---------------------------------------------------------------------------- #

.ineqx_params_manual <- function(data, ystat, ref, vcov) {

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

  # Validate ref for longitudinal
  if (type == "longitudinal" && is.null(ref)) {
    stop("'ref' (reference time period) is required for longitudinal decomposition ",
         "(data has ", n_times, " time periods)")
  }
  if (type == "longitudinal" && !is.null(ref) && !(ref %in% unique(data$time))) {
    stop("'ref' = ", ref, " not found in data$time. ",
         "Available time periods: ", paste(sort(unique(data$time)), collapse = ", "))
  }

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

  # Validate sigma0 > 0
  if (any(data$sigma0 <= 0)) {
    stop("'sigma0' must be strictly positive (it represents the baseline SD)")
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

  .make_ineqx_params(data, ystat, type, ref, vcov, has_time)
}


# ---------------------------------------------------------------------------- #
# Internal: Model extraction path
# ---------------------------------------------------------------------------- #

.ineqx_params_from_model <- function(model, data, treat, group,
                                      time, post, ystat, ref, vcov) {

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

  # -------------------------------------------------------------------- #
  # Generate predictions under counterfactual treatment assignments
  # -------------------------------------------------------------------- #

  data_cf0 <- data
  data_cf0[[treat]] <- 0
  pred_mu_0 <- predict(model, what = "mu", type = "response",
                       newdata = data_cf0, data = data)
  pred_sigma_0 <- predict(model, what = "sigma", type = "response",
                          newdata = data_cf0, data = data)

  if (length(treat_levels) == 1) {
    data_cf1 <- data
    data_cf1[[treat]] <- treat_levels[1]
    pred_mu_1 <- predict(model, what = "mu", type = "response",
                         newdata = data_cf1, data = data)
    pred_sigma_1 <- predict(model, what = "sigma", type = "response",
                            newdata = data_cf1, data = data)
  } else {
    pred_mu_1 <- rep(0, nrow(data))
    pred_sigma_1 <- rep(0, nrow(data))
    for (tl in treat_levels) {
      data_cf <- data
      data_cf[[treat]] <- tl
      p_mu <- predict(model, what = "mu", type = "response",
                      newdata = data_cf, data = data)
      p_sigma <- predict(model, what = "sigma", type = "response",
                         newdata = data_cf, data = data)
      freq <- sum(data[[treat]] == tl) / sum(data[[treat]] != 0)
      pred_mu_1 <- pred_mu_1 + freq * p_mu
      pred_sigma_1 <- pred_sigma_1 + freq * p_sigma
    }
  }

  # -------------------------------------------------------------------- #
  # Compute per group-time quantities
  # -------------------------------------------------------------------- #

  result_list <- list()
  for (t in time_levels) {
    for (g in group_levels) {
      idx <- data[[time_var]] == t & data[[group]] == g
      if (sum(idx) == 0) next

      mu0_g  <- mean(pred_mu_0[idx])
      sigma0_g <- mean(pred_sigma_0[idx])
      beta_g  <- mean(pred_mu_1[idx] - pred_mu_0[idx])
      lambda_g <- mean(log(pred_sigma_1[idx]) - log(pred_sigma_0[idx]))
      n_treated <- sum(data[[treat]][idx] != 0)

      result_list[[length(result_list) + 1]] <- data.frame(
        group = g, time = t, n_treated = n_treated,
        mu0 = mu0_g, sigma0 = sigma0_g,
        beta = beta_g, lambda = lambda_g,
        stringsAsFactors = FALSE
      )
    }
  }

  result_df <- do.call(rbind, result_list)

  # Normalize pi within each time period
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

  out_df <- result_df[, c("group", "time", "pi", "mu0", "sigma0", "beta", "lambda")]

  # Extract vcov if requested
  vcov_mat <- NULL
  if (vcov) {
    vcov_mat <- .extract_vcov_numerical(
      model = model, treat = treat, group = group,
      time_var = time_var, data = data,
      group_levels = group_levels, time_levels = time_levels,
      out_df = out_df
    )
  }

  # Remove time column if cross-sectional
  if (is.null(time)) {
    out_df$time <- NULL
  }

  # Use manual path for validation and S3 construction
  .ineqx_params_manual(data = out_df, ystat = ystat, ref = ref, vcov = vcov_mat)
}


# ---------------------------------------------------------------------------- #
# Internal: Construct the S3 object
# ---------------------------------------------------------------------------- #

.make_ineqx_params <- function(data, ystat, type, ref, vcov, has_time) {
  structure(
    list(
      data = data,
      ystat = ystat,
      type = type,
      ref = ref,
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
# Internal: Extract vcov via numerical Jacobian
# ---------------------------------------------------------------------------- #

#' @keywords internal
.extract_vcov_numerical <- function(model, treat, group, time_var, data,
                                     group_levels, time_levels, out_df) {

  raw_vcov <- tryCatch(
    stats::vcov(model),
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

  if (nrow(raw_vcov) != p || ncol(raw_vcov) != p) {
    vcov_mu <- tryCatch(stats::vcov(model, what = "mu"), error = function(e) NULL)
    vcov_sigma <- tryCatch(stats::vcov(model, what = "sigma"), error = function(e) NULL)
    if (!is.null(vcov_mu) && !is.null(vcov_sigma)) {
      p_mu <- length(mu_coefs)
      p_sigma <- length(sigma_coefs)
      raw_vcov <- matrix(0, p, p)
      raw_vcov[1:p_mu, 1:p_mu] <- vcov_mu
      raw_vcov[(p_mu + 1):p, (p_mu + 1):p] <- vcov_sigma
    } else {
      warning("Could not construct full vcov matrix. Returning NULL.")
      return(NULL)
    }
  }

  theta_base <- .coefs_to_theta(model, raw_coefs, treat, group, time_var,
                                 data, group_levels, time_levels)

  eps <- 1e-5
  n_theta <- length(theta_base)
  jacobian <- matrix(0, n_theta, p)

  for (k in seq_len(p)) {
    coefs_plus <- raw_coefs
    coefs_plus[k] <- coefs_plus[k] + eps
    theta_plus <- .coefs_to_theta(model, coefs_plus, treat, group, time_var,
                                   data, group_levels, time_levels)
    jacobian[, k] <- (theta_plus - theta_base) / eps
  }

  full_vcov <- jacobian %*% raw_vcov %*% t(jacobian)

  J <- length(group_levels)
  n_t <- length(time_levels)

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
# Internal: Map raw model coefficients to ineqx theta vector
# ---------------------------------------------------------------------------- #

#' @keywords internal
.coefs_to_theta <- function(model, coefs, treat, group, time_var, data,
                             group_levels, time_levels) {

  p_mu <- length(coef(model, what = "mu"))
  mu_coefs <- coefs[1:p_mu]
  sigma_coefs <- coefs[(p_mu + 1):length(coefs)]

  model_mod <- model
  model_mod$mu.coefficients <- mu_coefs
  model_mod$sigma.coefficients <- sigma_coefs

  data_cf0 <- data
  data_cf0[[treat]] <- 0
  pred_mu_0 <- predict(model_mod, what = "mu", type = "response",
                       newdata = data_cf0, data = data)
  pred_sigma_0 <- predict(model_mod, what = "sigma", type = "response",
                          newdata = data_cf0, data = data)

  data_cf1 <- data
  data_cf1[[treat]] <- 1
  pred_mu_1 <- predict(model_mod, what = "mu", type = "response",
                       newdata = data_cf1, data = data)
  pred_sigma_1 <- predict(model_mod, what = "sigma", type = "response",
                          newdata = data_cf1, data = data)

  theta <- c()
  for (t in time_levels) {
    betas <- mu0s <- lambdas <- log_sigma0s <- numeric(length(group_levels))
    for (j in seq_along(group_levels)) {
      g <- group_levels[j]
      idx <- data[[time_var]] == t & data[[group]] == g
      if (sum(idx) == 0) next

      mu0s[j] <- mean(pred_mu_0[idx])
      betas[j] <- mean(pred_mu_1[idx] - pred_mu_0[idx])
      lambdas[j] <- mean(log(pred_sigma_1[idx]) - log(pred_sigma_0[idx]))
      log_sigma0s[j] <- mean(log(pred_sigma_0[idx]))
    }
    theta <- c(theta, betas, mu0s, lambdas, log_sigma0s)
  }

  theta
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
    cat("  Reference period:", x$ref, "\n")
  }
  cat("\nParameter ranges:\n")
  d <- x$data
  cat("  pi:     [", round(min(d$pi), 4), ", ", round(max(d$pi), 4), "]\n", sep = "")
  cat("  mu0:    [", round(min(d$mu0), 2), ", ", round(max(d$mu0), 2), "]\n", sep = "")
  cat("  sigma0: [", round(min(d$sigma0), 2), ", ", round(max(d$sigma0), 2), "]\n", sep = "")
  cat("  beta:   [", round(min(d$beta), 4), ", ", round(max(d$beta), 4), "]\n", sep = "")
  cat("  lambda: [", round(min(d$lambda), 4), ", ", round(max(d$lambda), 4), "]\n", sep = "")
  if (!is.null(x$vcov)) {
    if (is.matrix(x$vcov)) {
      cat("\nVariance-covariance matrix: ", nrow(x$vcov), "x", ncol(x$vcov), "\n")
    } else {
      cat("\nVariance-covariance matrices: ", length(x$vcov), " time periods\n")
    }
  }
  invisible(x)
}
