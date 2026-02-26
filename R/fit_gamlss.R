# ============================================================================ #
# Convenience wrapper for fitting GAMLSS models
# ============================================================================ #

#' Fit a GAMLSS model for variance decomposition
#'
#' A convenience wrapper around \code{gamlss::gamlss()} for users who want
#' to fit a GAMLSS model without worrying about the syntax. The model
#' simultaneously estimates the conditional mean and conditional (log)
#' standard deviation, which are the inputs to the causal variance
#' decomposition.
#'
#' Users who want more control over the model specification should fit
#' their own \code{gamlss} model directly and use \code{\link{ineqx_params}}
#' to extract the decomposition parameters.
#'
#' @param formula_mu Formula for the mean equation (e.g.,
#'   \code{y ~ treat * group + controls})
#' @param formula_sigma Formula for the log-SD equation (e.g.,
#'   \code{~ treat * group + controls}). Note: one-sided formula.
#' @param data A data.frame
#' @param weights Optional character string naming the weight variable in
#'   \code{data}, or a numeric vector of weights
#' @param family A \code{gamlss.family} object specifying the distribution.
#'   Default: \code{gamlss.dist::NO()} (normal distribution).
#' @param ... Additional arguments passed to \code{gamlss::gamlss()}
#'
#' @return A fitted \code{gamlss} object
#'
#' @details
#' The GAMLSS framework models both the conditional mean and conditional
#' variance of the outcome distribution. The mean equation captures how
#' treatment affects average outcomes within each group. The sigma equation
#' captures how treatment affects the dispersion of outcomes within each group.
#'
#' For the simple difference estimator, the model is:
#' \deqn{\mu_{ij} = \alpha_j + \beta_{D,j} D_i + Z_i \gamma}
#' \deqn{\log(\sigma_{ij}) = \alpha_j^{(\sigma)} + \lambda_{D,j} D_i + Z_i \gamma^{(\sigma)}}
#'
#' For the DiD estimator, add pre/post interactions:
#' \deqn{\mu_{ij} = \alpha_j + \beta_{D,j} D_i + \beta_P P_i + \beta_{DP,j} D_i P_i + Z_i \gamma}
#'
#' @seealso \code{\link{ineqx_params}} for extracting decomposition parameters
#'   from the fitted model
#'
#' @keywords internal
fit_ineqx_model <- function(formula_mu, formula_sigma, data,
                             weights = NULL, family = NULL, ...) {

  if (!requireNamespace("gamlss", quietly = TRUE)) {
    stop("Package 'gamlss' is required for fit_ineqx_model(). ",
         "Install it with install.packages('gamlss')")
  }

  # Default family
  if (is.null(family)) {
    family <- gamlss.dist::NO()
  }

  # Handle weights
  w <- NULL
  if (!is.null(weights)) {
    if (is.character(weights)) {
      if (!weights %in% names(data)) {
        stop("Weight variable '", weights, "' not found in data")
      }
      w <- data[[weights]]
    } else {
      w <- weights
    }
  }

  # Fit the model
  if (is.null(w)) {
    gamlss::gamlss(
      formula = formula_mu,
      sigma.formula = formula_sigma,
      data = data,
      family = family,
      ...
    )
  } else {
    gamlss::gamlss(
      formula = formula_mu,
      sigma.formula = formula_sigma,
      weights = w,
      data = data,
      family = family,
      ...
    )
  }
}
