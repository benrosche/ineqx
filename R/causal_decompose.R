# ============================================================================ #
# Causal variance decomposition
# Cross-sectional and longitudinal with configurable ordering
# ============================================================================ #

# ---------------------------------------------------------------------------- #
# Cross-sectional causal decomposition
# ---------------------------------------------------------------------------- #

#' Cross-sectional causal variance decomposition
#'
#' Decomposes the treatment effect on total variance (or CV^2) into
#' between-group and within-group components at a single timepoint.
#'
#' The between-group treatment effect is:
#' \deqn{\delta_B^D = \text{Var}_\pi(\beta_D) + 2\,\text{Cov}_\pi(\mu(0), \beta_D)}
#'
#' The within-group treatment effect is:
#' \deqn{\delta_W^D = \sum_j \pi_j \sigma_j^2(0) [\exp(2\lambda_{D,j}) - 1]}
#'
#' @param params An \code{ineqx_params} object. For cross-sectional
#'   decomposition, this should contain a single time period (or no time column).
#'
#' @return An object of class \code{"ineqx_causal_cross"} containing:
#' \describe{
#'   \item{delta_B}{Scalar, between-group treatment effect on inequality}
#'   \item{delta_W}{Scalar, within-group treatment effect on inequality}
#'   \item{delta_total}{Scalar, total treatment effect (delta_B + delta_W)}
#'   \item{components}{List of sub-components for interpretive analysis}
#'   \item{by_group}{data.frame of group-level contributions}
#'   \item{params}{The input ineqx_params object}
#' }
#'
#' @keywords internal
causal_decompose_cross <- function(params) {

  stopifnot(inherits(params, "ineqx_params"))

  d <- params$data
  ystat <- params$ystat

  # If longitudinal, use only the first time period (or ref)
  if (params$type == "longitudinal") {
    t_use <- if (!is.null(params$ref)) params$ref else params$times[1]
    d <- d[d$time == t_use, ]
  }

  pi_j     <- d$pi
  mu0_j    <- d$mu0
  sigma0_j <- d$sigma0
  beta_j   <- d$beta
  lambda_j <- d$lambda
  groups   <- d$group

  # Compute delta_B and delta_W using exact plug-in formulas
  wb <- compute_delta_WB(pi_j, mu0_j, sigma0_j, beta_j, lambda_j, ystat)

  # Sub-components for interpretation (Var decomposition)
  mu_bar <- sum(pi_j * mu0_j)
  beta_bar <- sum(pi_j * beta_j)
  var_pi_beta <- sum(pi_j * (beta_j - beta_bar)^2)
  cov_pi_mu_beta <- sum(pi_j * (mu0_j - mu_bar) * (beta_j - beta_bar))

  f_j <- exp(2 * lambda_j) - 1
  sigma2_bar <- sum(pi_j * sigma0_j^2)
  f_bar <- sum(pi_j * f_j)
  cov_pi_sigma2_f <- sum(pi_j * (sigma0_j^2 - sigma2_bar) * (f_j - f_bar))

  # Group-level contributions
  by_group <- data.frame(
    group = groups,
    pi = pi_j,
    mu0 = mu0_j,
    sigma0 = sigma0_j,
    beta = beta_j,
    lambda = lambda_j,
    f = f_j,
    contrib_W = pi_j * sigma0_j^2 * f_j,
    stringsAsFactors = FALSE
  )

  # Compute delta method SEs if vcov is available
  se_list <- NULL
  if (!is.null(params$vcov)) {
    se_list <- delta_method_se(params, type = "cross")
  }

  structure(
    list(
      delta_B = wb$delta_B,
      delta_W = wb$delta_W,
      delta_total = wb$delta_total,
      se = se_list,
      components = list(
        Var_pi_beta = var_pi_beta,
        Cov_pi_mu_beta = cov_pi_mu_beta,
        mean_sigma2_0 = sigma2_bar,
        mean_f = f_bar,
        Cov_pi_sigma2_f = cov_pi_sigma2_f
      ),
      by_group = by_group,
      ystat = ystat,
      params = params
    ),
    class = "ineqx_causal_cross"
  )
}


# ---------------------------------------------------------------------------- #
# Longitudinal causal decomposition
# ---------------------------------------------------------------------------- #

#' Longitudinal causal variance decomposition
#'
#' Decomposes the change in treatment effect on inequality from a reference
#' period to each subsequent period into behavioral, compositional, and
#' pre-treatment components, using sequential parameter switching.
#'
#' The six components are:
#' \describe{
#'   \item{Delta_beta}{Changing treatment effects on means (behavioral, between)}
#'   \item{Delta_lambda}{Changing treatment effects on variances (behavioral, within)}
#'   \item{Delta_pi_B}{Changing group composition (compositional, between)}
#'   \item{Delta_pi_W}{Changing group composition (compositional, within)}
#'   \item{Delta_mu}{Changing pre-treatment means (pre-treatment, between)}
#'   \item{Delta_sigma}{Changing pre-treatment variances (pre-treatment, within)}
#' }
#'
#' These group into four combined components:
#' behavioral = Delta_beta + Delta_lambda,
#' compositional = Delta_pi_B + Delta_pi_W,
#' pre-treatment = Delta_mu + Delta_sigma.
#'
#' @param params An \code{ineqx_params} object with multiple time periods
#' @param order Character vector of length 3, a permutation of
#'   \code{c("behavioral", "compositional", "pretreatment")} specifying
#'   the order in which parameters are switched from baseline to time-t values.
#'   This ordering is applied in parallel to both the between-group
#'   (beta, pi, mu) and within-group (lambda, pi, sigma) sub-decompositions.
#'   Default: \code{c("behavioral", "compositional", "pretreatment")}.
#'
#' @return An object of class \code{"ineqx_causal_longit"} containing:
#' \describe{
#'   \item{results}{List keyed by time period, each containing the 6 components
#'     plus 3 combined components and the total}
#'   \item{order}{The ordering used}
#'   \item{ystat}{The inequality measure}
#'   \item{params}{The input ineqx_params object}
#' }
#'
#' @keywords internal
causal_decompose_longit <- function(params,
                                     order = c("behavioral", "compositional",
                                               "pretreatment")) {

  stopifnot(inherits(params, "ineqx_params"))
  stopifnot(params$type == "longitudinal")

  # Validate ordering
  valid_components <- c("behavioral", "compositional", "pretreatment")
  order <- match.arg(order, valid_components, several.ok = TRUE)
  if (length(order) != 3 || !setequal(order, valid_components)) {
    stop("'order' must be a permutation of c('behavioral', 'compositional', 'pretreatment')")
  }

  ref_time <- params$ref
  ystat <- params$ystat
  d <- params$data

  # Reference period data
  d0 <- d[d$time == ref_time, ]

  results <- list()
  for (t in setdiff(params$times, ref_time)) {
    dt <- d[d$time == t, ]

    # Between-group sub-decomposition
    between <- .decompose_between_sequential(d0, dt, order, ystat)

    # Within-group sub-decomposition
    within <- .decompose_within_sequential(d0, dt, order, ystat)

    # Cross-sectional decomposition at t0 and t for context
    wb_t0 <- compute_delta_WB(d0$pi, d0$mu0, d0$sigma0, d0$beta, d0$lambda, ystat)
    wb_t  <- compute_delta_WB(dt$pi, dt$mu0, dt$sigma0, dt$beta, dt$lambda, ystat)

    result_t <- list(
      time = t,
      # 6 detailed components
      Delta_beta   = between$behavioral,
      Delta_lambda = within$behavioral,
      Delta_pi_B   = between$compositional,
      Delta_pi_W   = within$compositional,
      Delta_mu     = between$pretreatment,
      Delta_sigma  = within$pretreatment,
      # 3 combined components (paper Eq. 18)
      Delta_behavioral    = between$behavioral + within$behavioral,
      Delta_compositional = between$compositional + within$compositional,
      Delta_pretreatment  = between$pretreatment + within$pretreatment,
      # Between and within totals
      Delta_B = between$delta_total,
      Delta_W = within$delta_total,
      # Grand total
      Delta_total = between$delta_total + within$delta_total,
      # Cross-sectional at each time for reference
      delta_B_t0 = wb_t0$delta_B,
      delta_W_t0 = wb_t0$delta_W,
      delta_B_t  = wb_t$delta_B,
      delta_W_t  = wb_t$delta_W
    )

    results[[as.character(t)]] <- result_t
  }

  # Compute delta method SEs if vcov is available
  se_list <- NULL
  if (!is.null(params$vcov)) {
    se_list <- delta_method_se(params, type = "longit", order = order)
  }

  structure(
    list(
      results = results,
      se = se_list,
      order = order,
      ystat = ystat,
      ref = ref_time,
      params = params
    ),
    class = "ineqx_causal_longit"
  )
}


# ---------------------------------------------------------------------------- #
# Internal: Sequential parameter switching for between-group
# ---------------------------------------------------------------------------- #

#' @keywords internal
.decompose_between_sequential <- function(d0, dt, order, ystat) {

  # Map abstract component names to between-group parameter names
  param_map <- list(
    behavioral    = list(param = "beta", from = d0$beta,  to = dt$beta),
    compositional = list(param = "pi",   from = d0$pi,    to = dt$pi),
    pretreatment  = list(param = "mu",   from = d0$mu0,   to = dt$mu0)
  )

  # Start with all parameters at t0 values
  current <- list(beta = d0$beta, pi = d0$pi, mu = d0$mu0)

  # Compute delta_B^D = B(pi, mu+beta) - B(pi, mu) at current parameter values
  .eval_delta_B <- function(pars) {
    compute_delta_B(pars$pi, pars$mu, pars$beta, ystat)
  }

  prev_value <- .eval_delta_B(current)
  baseline_value <- prev_value
  components <- list()

  for (step in order) {
    p <- param_map[[step]]$param
    current[[p]] <- param_map[[step]]$to
    new_value <- .eval_delta_B(current)
    components[[step]] <- new_value - prev_value
    prev_value <- new_value
  }

  components$delta_total <- prev_value - baseline_value
  components
}


# ---------------------------------------------------------------------------- #
# Internal: Sequential parameter switching for within-group
# ---------------------------------------------------------------------------- #

#' @keywords internal
.decompose_within_sequential <- function(d0, dt, order, ystat) {

  # Map abstract component names to within-group parameter names
  param_map <- list(
    behavioral    = list(param = "lambda", from = d0$lambda, to = dt$lambda),
    compositional = list(param = "pi",     from = d0$pi,     to = dt$pi),
    pretreatment  = list(param = "sigma",  from = d0$sigma0, to = dt$sigma0)
  )

  # For CV2, we also need mu — but within-group CV2 depends on total means
  # which in turn depend on pi and mu0. For the within-group switching,
  # mu0 is not being switched (that's in between-group). So we need to
  # decide which mu0 to use. Since this is the within-group sub-decomposition,
  # we evaluate W(d=1) - W(d=0) using the group means.
  # For CV2: we need mu0 and beta to compute the grand mean.
  # We use mu0 and beta as they are at the current switching state.
  # Since mu0 is NOT switched in the within sub-decomposition,
  # we need to handle it specially for CV2.

  # Start with all parameters at t0 values
  current <- list(lambda = d0$lambda, pi = d0$pi, sigma = d0$sigma0)

  # For CV2, we need mu0 and beta. These change in the between-group
  # sub-decomposition but are held constant here.
  # Use time-t values for mu0 and beta (since between is decomposed separately)
  # Actually, for the sequential switching, we need the mu0 and beta
  # that correspond to the same switching state. Since we switch in the
  # same order for both sub-decompositions, and within only touches
  # lambda/pi/sigma, we should use:
  # - mu0 at the current time's switching state for the between decomposition
  # But since between and within are independent sub-decompositions,
  # we use the FACTUAL values at each time to evaluate the cross-sectional
  # within-group treatment effect.

  # For Var: delta_W = sum(pi * sigma^2 * f) where f = exp(2*lambda) - 1
  # For CV2: delta_W = [sum(pi * (sigma*exp(lambda))^2) / gmu1^2] -
  #                    [sum(pi * sigma^2) / gmu0^2]
  # where gmu1 = sum(pi*(mu0+beta)), gmu0 = sum(pi*mu0)
  # For within sub-decomposition, mu0 and beta are NOT switched,
  # so we use the corresponding values from dt for the "destination"

  .eval_delta_W <- function(pars) {
    f <- exp(2 * pars$lambda) - 1
    if (ystat == "Var") {
      # For Var, delta_W = sum(pi * sigma^2 * f) is independent of mu/beta,
      # so the between and within sub-decompositions are truly independent.
      sum(pars$pi * pars$sigma^2 * f)
    } else {
      # For CV2, delta_W depends on grand means (via mu0 and beta).
      # Since mu0 and beta are switched in the between-group sub-decomposition
      # but NOT here, we use time-t values as the "destination" state.
      # This ensures the telescoping sums work correctly when between and
      # within totals are added, though individual CV2 components should be
      # interpreted with care.
      mu0 <- dt$mu0
      beta <- dt$beta
      gmu1 <- sum(pars$pi * (mu0 + beta))
      gmu0 <- sum(pars$pi * mu0)
      sigma1 <- pars$sigma * exp(pars$lambda)
      sum(pars$pi * sigma1^2) / gmu1^2 - sum(pars$pi * pars$sigma^2) / gmu0^2
    }
  }

  prev_value <- .eval_delta_W(current)
  baseline_value <- prev_value
  components <- list()

  for (step in order) {
    p <- param_map[[step]]$param
    current[[p]] <- param_map[[step]]$to
    new_value <- .eval_delta_W(current)
    components[[step]] <- new_value - prev_value
    prev_value <- new_value
  }

  components$delta_total <- prev_value - baseline_value
  components
}
