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
#' @param na.action A function that handles NAs in \code{data}. Applied to the
#'   subset of \code{data} containing the variables referenced by
#'   \code{formula_mu}, \code{formula_sigma}, and \code{weights} before fitting.
#'   Defaults to \code{\link[stats]{na.omit}}, matching \code{\link{ineqx_params}}
#'   and the integrated \code{\link{ineqx}()} path so the split-step workflow
#'   keeps both steps on identical rows. Pass \code{na.fail} to error on NAs.
#' @param transform Either \code{"identity"} (default) or \code{"log"}.
#'   When \code{"log"}, the response column on the LHS of \code{formula_mu}
#'   is log-transformed in a local copy of \code{data} before fitting, and
#'   the returned model is tagged with \code{attr(model, "ineqx_transform")
#'   = "log"}. Downstream \code{\link{ineqx_params}} and \code{\link{ineqx}}
#'   read this tag, which is what lets the split-step workflow compute
#'   \eqn{V_L} (variance of log earnings) via \code{ystat = "VL"}. Requires
#'   the LHS of \code{formula_mu} to be a simple variable name (e.g.,
#'   \code{y}, not \code{log(y)} or \code{I(y + 1)}) and that column to be
#'   strictly positive.
#' @param ... Additional arguments passed to \code{gamlss::gamlss()}
#'
#' @return A fitted \code{gamlss} object. When \code{transform = "log"}, the
#'   object carries \code{attr(., "ineqx_transform") = "log"} so downstream
#'   functions know the fit lives on the log scale.
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
#' @section Split-step workflow:
#' \code{fit_ineqx_model()} is the entry point to the split-step workflow,
#' where the GAMLSS fit is cached once and many decompositions (varying
#' \code{ystat} and \code{ref}) are derived cheaply on top:
#' \preformatted{
#' model  <- fit_ineqx_model(y ~ treat * group * time + controls,
#'                                ~ treat * group * time + controls,
#'                            data = d)
#' params <- ineqx_params(model = model, data = d,
#'                        treat = "treat", group = "group", time = "time",
#'                        ystat = "Var", vcov = TRUE)
#' fit_var <- ineqx(params = params, ystat = "Var", ref = 1980)
#' fit_cv2 <- ineqx(params = params, ystat = "CV2", ref = 1980)
#' }
#' For the stepwise DiD, pass \code{post = "post01"} (or whatever your
#' pre/post column is) to \code{ineqx_params()}. For \eqn{V_L}, fit with
#' \code{transform = "log"} and pass \code{ystat = "VL"} to \code{ineqx()}.
#'
#' @section Single-level grouping variables:
#' A factor or character predictor with only one observed level cannot enter a
#' model formula (\code{model.matrix()} errors with "contrasts can be applied
#' only to factors with 2 or more levels"). Any such term is automatically
#' dropped before fitting, with a message, since the between-group component it
#' would identify is zero by construction. The rest of the model is estimated
#' normally and the decomposition returns \code{VarB = 0} (within equals total).
#'
#' @seealso \code{\link{ineqx_params}} for extracting decomposition parameters
#'   from the fitted model; \code{\link{ineqx}} for running the decomposition.
#'
#' @export
fit_ineqx_model <- function(formula_mu, formula_sigma, data,
                             weights = NULL, family = NULL,
                             transform = c("identity", "log"),
                             na.action = stats::na.omit, ...) {

  if (!requireNamespace("gamlss", quietly = TRUE)) {
    stop("Package 'gamlss' is required for fit_ineqx_model(). ",
         "Install it with install.packages('gamlss')")
  }

  transform <- match.arg(transform)

  # ------------------------------------------------------------------ #
  # Materialize weights into the data as `.ineqx_w` BEFORE applying
  # na.action, so that row-dropping stays aligned when a numeric weight
  # vector is supplied (the vector becomes a column and is subset with the
  # rest of the data).
  # ------------------------------------------------------------------ #
  if (!is.null(weights)) {
    if (is.character(weights)) {
      if (!weights %in% names(data)) {
        stop("Weight variable '", weights, "' not found in data")
      }
      data$.ineqx_w <- data[[weights]]
    } else {
      if (length(weights) != nrow(data)) {
        stop("Numeric 'weights' must have length nrow(data) (",
             length(weights), " != ", nrow(data), ").")
      }
      data$.ineqx_w <- weights
    }
  }

  # ------------------------------------------------------------------ #
  # Apply na.action to model-relevant columns. Mirrors ineqx_params() so
  # the two split-step entry points handle NAs identically and gamlss does
  # not error on incomplete rows.
  # ------------------------------------------------------------------ #
  vars <- unique(c(all.vars(formula_mu), all.vars(formula_sigma)))
  if (!is.null(weights)) vars <- c(vars, ".ineqx_w")
  vars <- intersect(vars, names(data))
  data <- na.action(data[, vars, drop = FALSE])

  # ------------------------------------------------------------------ #
  # Drop terms built on single-level factor/character variables. A
  # one-level factor cannot enter the model (model.matrix() errors with
  # "contrasts can be applied only to factors with 2 or more levels") and
  # the between-group component it would identify is zero by construction,
  # so we strip those terms and let the rest of the model estimate.
  # ------------------------------------------------------------------ #
  drop_mu    <- .drop_constant_factor_terms(formula_mu, data)
  drop_sigma <- .drop_constant_factor_terms(formula_sigma, data)
  formula_mu    <- drop_mu$formula
  formula_sigma <- drop_sigma$formula
  dropped_vars  <- unique(c(drop_mu$dropped, drop_sigma$dropped))
  if (length(dropped_vars) > 0L) {
    message("The grouping variable ('", paste(dropped_vars, collapse = "', '"),
            "') has only one level. The between-group component is therefore ",
            "zero by construction, and the group-level term will be dropped ",
            "from the model.")
  }

  # ------------------------------------------------------------------ #
  # transform = "log": log-transform the LHS column of formula_mu in a
  # local copy of data. Mirrors the auto-log path inside ineqx() for the
  # integrated VL workflow (see ineqx.R:170-201), so downstream
  # ineqx_params() + ineqx(params, ystat = "VL") yields identical numbers
  # to the integrated ineqx(y, data, ystat = "VL") call.
  # ------------------------------------------------------------------ #
  if (transform == "log") {
    if (length(formula_mu) < 3L) {
      stop("transform = 'log' requires a two-sided formula_mu with a ",
           "response on the LHS (e.g. y ~ treat * group + controls).")
    }
    lhs <- formula_mu[[2L]]
    if (!is.name(lhs)) {
      stop("transform = 'log' requires the LHS of formula_mu to be a ",
           "simple variable name (e.g. y), not an expression like ",
           "log(y) or I(y + 1). Got: ", deparse(lhs))
    }
    y_var <- as.character(lhs)
    if (!y_var %in% names(data)) {
      stop("transform = 'log': response variable '", y_var,
           "' not found in data.")
    }
    yvals <- data[[y_var]]
    bad <- !is.na(yvals) & yvals <= 0
    if (any(bad)) {
      stop("transform = 'log' requires strictly positive '", y_var,
           "'; found ", sum(bad), " non-positive value(s).")
    }
    data[[y_var]] <- log(yvals)
  }

  # Default family
  if (is.null(family)) {
    family <- gamlss.dist::NO()
  }

  # Weights were materialized into `data$.ineqx_w` above (before na.action).
  # gamlss resolves the `.ineqx_w` symbol from the data frame via
  # match.call() + model.frame().

  # Fit the model
  if (is.null(weights)) {
    model <- gamlss::gamlss(
      formula = formula_mu,
      sigma.formula = formula_sigma,
      data = data,
      family = family,
      ...
    )
  } else {
    # `.ineqx_w` must exist both as a column in `data` (for model.frame()
    # resolution in predict()) AND as a local variable (for R's standard
    # argument evaluation when gamlss() is called).
    .ineqx_w <- data$.ineqx_w  # nolint
    model <- gamlss::gamlss(
      formula = formula_mu,
      sigma.formula = formula_sigma,
      weights = .ineqx_w,
      data = data,
      family = family,
      ...
    )
  }

  # Tag the model with the scale on which the response was modeled, so
  # that ineqx_params() can propagate this to params and ineqx() knows
  # whether ystat = "VL" is valid for a params-supplied call.
  attr(model, "ineqx_transform") <- transform
  model
}


# ---------------------------------------------------------------------------- #
# Internal: Drop formula terms built on single-level factor/character variables
#
# A factor/character/logical variable with only one observed level cannot enter
# a model formula -- model.matrix() errors with "contrasts can be applied only
# to factors with 2 or more levels". Such a term also carries no information
# (the between-group component it would identify is zero by construction), so we
# remove every term that references it and keep the remaining structure.
# Returns list(formula, dropped) where `dropped` names the offending variables.
# ---------------------------------------------------------------------------- #

.drop_constant_factor_terms <- function(formula, data) {
  tt <- stats::terms(formula)
  term_labels <- attr(tt, "term.labels")
  if (length(term_labels) == 0L) {
    return(list(formula = formula, dropped = character(0)))
  }

  is_single <- function(col) {
    (is.factor(col) || is.character(col) || is.logical(col)) &&
      length(unique(col[!is.na(col)])) < 2L
  }
  single_vars <- names(data)[vapply(data, is_single, logical(1))]
  if (length(single_vars) == 0L) {
    return(list(formula = formula, dropped = character(0)))
  }

  term_vars <- lapply(term_labels, function(tl) all.vars(stats::reformulate(tl)))
  drop_term <- vapply(term_vars, function(v) any(v %in% single_vars), logical(1))
  if (!any(drop_term)) {
    return(list(formula = formula, dropped = character(0)))
  }

  dropped <- intersect(unique(unlist(term_vars[drop_term])), single_vars)
  kept    <- term_labels[!drop_term]

  has_intercept <- attr(tt, "intercept") == 1L
  resp <- if (attr(tt, "response") == 1L) deparse(formula[[2L]]) else NULL

  rhs <- kept
  if (!has_intercept) rhs <- c(rhs, "0")
  if (length(rhs) == 0L) rhs <- if (has_intercept) "1" else "0"

  new_f <- stats::reformulate(rhs, response = resp)
  environment(new_f) <- environment(formula)
  list(formula = new_f, dropped = dropped)
}
