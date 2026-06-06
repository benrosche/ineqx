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
#'   \item{tau_B}{Scalar, between-group treatment effect on inequality}
#'   \item{tau_W}{Scalar, within-group treatment effect on inequality}
#'   \item{tau_total}{Scalar, total treatment effect (tau_B + tau_W)}
#'   \item{components}{List of sub-components for interpretive analysis}
#'   \item{by_group}{data.frame of group-level contributions}
#'   \item{params}{The input ineqx_params object}
#' }
#'
#' @keywords internal
causal_decompose_cross <- function(params, ref = NULL) {

  stopifnot(inherits(params, "ineqx_params"))

  d <- params$data
  ystat <- params$ystat

  # If longitudinal, use only the specified ref period (or first time)
  if (params$type == "longitudinal") {
    t_use <- if (!is.null(ref)) ref else params$times[1]
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

  # Sub-components for interpretation
  comps <- .compute_cross_components(pi_j, mu0_j, sigma0_j, beta_j, lambda_j, ystat)

  # Group-level contributions (these sum exactly to delta_W and delta_B)
  mu1_j <- mu0_j + beta_j
  sigma1_j <- sigma0_j * exp(lambda_j)
  f_j <- exp(2 * lambda_j) - 1

  if (ystat == "Var") {
    mu0_bar <- sum(pi_j * mu0_j)
    mu1_bar <- sum(pi_j * mu1_j)
    contrib_B_j <- pi_j * (mu1_j - mu1_bar)^2 - pi_j * (mu0_j - mu0_bar)^2
    contrib_W_j <- pi_j * (sigma1_j^2 - sigma0_j^2)
  } else {
    # CV2: B1_j - B0_j and W1_j - W0_j per group
    mu0_bar <- sum(pi_j * mu0_j)
    mu1_bar <- sum(pi_j * mu1_j)
    contrib_B_j <- pi_j * ((mu1_j / mu1_bar - 1)^2 - (mu0_j / mu0_bar - 1)^2)
    contrib_W_j <- pi_j * ((sigma1_j / mu1_bar)^2 - (sigma0_j / mu0_bar)^2)
  }

  by_group <- data.frame(
    group = groups,
    pi = pi_j,
    mu0 = mu0_j,
    sigma0 = sigma0_j,
    beta = beta_j,
    lambda = lambda_j,
    f = f_j,
    contrib_B = contrib_B_j,
    contrib_W = contrib_W_j,
    stringsAsFactors = FALSE
  )

  # Compute delta method SEs if vcov is available
  se_list <- NULL
  lambda_test <- NULL
  beta_test <- NULL
  if (!is.null(params$vcov)) {
    se_list <- delta_method_se(params, type = "cross", ref = ref)

    # Joint Wald test: H0: all lambda_j = 0
    J <- length(lambda_j)
    V <- .get_vcov_for_time(params,
           if (params$type == "longitudinal") ref else NULL)
    lambda_idx <- (2 * J + 1):(3 * J)
    V_lambda <- V[lambda_idx, lambda_idx, drop = FALSE]
    lambda_test <- tryCatch({
      W_stat <- as.numeric(t(lambda_j) %*% solve(V_lambda) %*% lambda_j)
      p_value <- pchisq(W_stat, df = J, lower.tail = FALSE)
      list(statistic = W_stat, df = J, p_value = p_value)
    }, error = function(e) NULL)

    # Joint Wald test: H0: beta_j homogeneous across groups (mean effect equal).
    # Scale-aware -- identity fit tests a constant ABSOLUTE effect; log fit tests
    # a constant PROPORTIONAL (geometric-mean) effect. beta occupies the first J
    # slots of the parameter vector (beta, mu0, lambda, log sigma0).
    beta_scale <- if (!is.null(params$transform)) params$transform else "identity"
    beta_test <- tryCatch({
      if (J < 2L) NULL else {
        V_beta <- V[1:J, 1:J, drop = FALSE]
        Cmat   <- cbind(-1, diag(J - 1L))      # (J-1) x J: beta_g - beta_1
        Cb     <- as.numeric(Cmat %*% beta_j)
        W_stat <- as.numeric(t(Cb) %*% solve(Cmat %*% V_beta %*% t(Cmat)) %*% Cb)
        list(statistic = W_stat, df = J - 1L,
             p_value = pchisq(W_stat, df = J - 1L, lower.tail = FALSE),
             scale = beta_scale)
      }
    }, error = function(e) NULL)
  }

  # Total inequality under control vs treatment
  if (ystat == "Var") {
    mu0_bar_all <- sum(pi_j * mu0_j)
    mu1_bar_all <- sum(pi_j * mu1_j)
    ineq_control <- sum(pi_j * (mu0_j^2 + sigma0_j^2)) - mu0_bar_all^2
    ineq_treated <- sum(pi_j * (mu1_j^2 + sigma1_j^2)) - mu1_bar_all^2
  } else {
    ineq_control <- CV2T_pi(pi_j, mu0_j, sigma0_j)
    ineq_treated <- CV2T_pi(pi_j, mu1_j, sigma1_j)
  }

  structure(
    list(
      tau_B = wb$tau_B,
      tau_W = wb$tau_W,
      tau_total = wb$tau_total,
      ineq_control = ineq_control,
      ineq_treated = ineq_treated,
      se = se_list,
      lambda_test = lambda_test,
      beta_test = beta_test,
      components = comps,
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
#' pre-treatment components, using a unified sequential parameter switching
#' scheme that tracks each switch's contribution to both \eqn{\tau_B} and
#' \eqn{\tau_W} simultaneously.
#'
#' For each parameter \eqn{p \in \{\beta, \lambda, \pi, \mu, \sigma\}}, the
#' result reports two split components:
#' \describe{
#'   \item{Delta_<p>_B}{Contribution of switching p (from t0 to t) to the change in \eqn{\tau_B}}
#'   \item{Delta_<p>_W}{Contribution of switching p (from t0 to t) to the change in \eqn{\tau_W}}
#' }
#'
#' For the variance (\code{ystat = "Var"}), \eqn{\tau_B} depends only on
#' \eqn{(\pi, \mu, \beta)} and \eqn{\tau_W} only on \eqn{(\pi, \sigma, \lambda)},
#' so the off-diagonal split parts (\code{Delta_beta_W, Delta_lambda_B,
#' Delta_mu_W, Delta_sigma_B}) are exactly zero. For CV\eqn{^2}, both
#' \eqn{\tau_B} and \eqn{\tau_W} share the grand-mean denominator, so all five
#' parameters can contribute to both sides; the split components capture this
#' coupling exactly.
#'
#' By construction, the split components sum to the cross-sectional change:
#' \eqn{\sum_p \text{Delta\_<p>\_B} = \tau_B(t) - \tau_B(t_0)} and similarly
#' for the W side, so \code{Delta_total} equals
#' \eqn{(\tau_B(t) + \tau_W(t)) - (\tau_B(t_0) + \tau_W(t_0))}.
#'
#' Aggregate parameter components (\code{Delta_beta, Delta_lambda, Delta_mu,
#' Delta_sigma}) are reported as the sum of their B and W parts; for V they
#' equal the previous (single-side) values exactly. Compositional effects
#' remain split as \code{Delta_pi_B} and \code{Delta_pi_W}.
#'
#' @param params An \code{ineqx_params} object with multiple time periods
#' @param order Character vector of length 3, a permutation of
#'   \code{c("behavioral", "compositional", "pretreatment")} specifying
#'   the order in which parameter groups are switched from baseline to
#'   time-t values. The mapping to the underlying 5-parameter sequence is:
#'   \code{behavioral -> (beta, lambda)}, \code{compositional -> (pi)},
#'   \code{pretreatment -> (mu, sigma)}, with \eqn{\beta} switched before
#'   \eqn{\lambda} and \eqn{\mu} switched before \eqn{\sigma} within their
#'   meta-levels. Use \code{order = "shapley"} to average over all 6
#'   meta-orderings. Default: \code{c("behavioral", "compositional", "pretreatment")}.
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
                                               "pretreatment"),
                                     ref = NULL) {

  stopifnot(inherits(params, "ineqx_params"))
  stopifnot(params$type == "longitudinal")

  # Validate ordering
  valid_components <- c("behavioral", "compositional", "pretreatment")
  order <- match.arg(order, valid_components, several.ok = TRUE)
  if (length(order) != 3 || !setequal(order, valid_components)) {
    stop("'order' must be a permutation of c('behavioral', 'compositional', 'pretreatment')")
  }

  if (is.null(ref)) {
    stop("'ref' is required for longitudinal decomposition.")
  }
  ref_time <- ref
  ystat <- params$ystat
  d <- params$data

  # Reference period data
  d0 <- d[d$time == ref_time, ]

  order_5 <- .expand_order_to_5(order)

  results <- list()
  for (t in setdiff(params$times, ref_time)) {
    dt <- d[d$time == t, ]

    # Unified 5-parameter sequential decomposition: each switch contributes
    # to both tau_B and tau_W (off-diagonal contributions are 0 for V).
    parts <- .decompose_sequential(d0, dt, order_5, ystat)

    # Cross-sectional decomposition at t0 and t for context
    wb_t0 <- compute_delta_WB(d0$pi, d0$mu0, d0$sigma0, d0$beta, d0$lambda, ystat)
    wb_t  <- compute_delta_WB(dt$pi, dt$mu0, dt$sigma0, dt$beta, dt$lambda, ystat)

    # Cross-sectional sub-components at each time (exact, not linearized)
    comps_t0 <- .compute_cross_components(d0$pi, d0$mu0, d0$sigma0, d0$beta, d0$lambda, ystat)
    comps_t  <- .compute_cross_components(dt$pi, dt$mu0, dt$sigma0, dt$beta, dt$lambda, ystat)

    result_t <- list(
      time = t,
      # Split components: B and W parts of each parameter switch.
      # For V, off-diagonals (Delta_beta_W, Delta_lambda_B, Delta_mu_W,
      # Delta_sigma_B) are exactly 0.
      Delta_beta_B   = parts$beta$B,
      Delta_beta_W   = parts$beta$W,
      Delta_lambda_B = parts$lambda$B,
      Delta_lambda_W = parts$lambda$W,
      Delta_pi_B     = parts$pi$B,
      Delta_pi_W     = parts$pi$W,
      Delta_mu_B     = parts$mu$B,
      Delta_mu_W     = parts$mu$W,
      Delta_sigma_B  = parts$sigma$B,
      Delta_sigma_W  = parts$sigma$W,
      # Aggregate parameter components (sums of split parts).
      # For V these equal the previous single-side values exactly.
      Delta_beta   = parts$beta$B   + parts$beta$W,
      Delta_lambda = parts$lambda$B + parts$lambda$W,
      Delta_mu     = parts$mu$B     + parts$mu$W,
      Delta_sigma  = parts$sigma$B  + parts$sigma$W,
      # 3 combined components
      Delta_behavioral    = parts$beta$B   + parts$beta$W +
                            parts$lambda$B + parts$lambda$W,
      Delta_compositional = parts$pi$B + parts$pi$W,
      Delta_pretreatment  = parts$mu$B    + parts$mu$W +
                            parts$sigma$B + parts$sigma$W,
      # 4-component grouping (paper Eq. 28)
      Delta_pi  = parts$pi$B + parts$pi$W,
      Delta_pre = parts$mu$B    + parts$mu$W +
                  parts$sigma$B + parts$sigma$W,
      # Between and within totals (= cross-sectional changes by construction)
      Delta_B = parts$totals$delta_B,
      Delta_W = parts$totals$delta_W,
      # Grand total (= (tau_B(t) + tau_W(t)) - (tau_B(t0) + tau_W(t0)))
      Delta_total = parts$totals$delta_T,
      # Cross-sectional at each time for reference
      tau_B_t0 = wb_t0$tau_B,
      tau_W_t0 = wb_t0$tau_W,
      tau_B_t  = wb_t$tau_B,
      tau_W_t  = wb_t$tau_W,
      # Cross-sectional sub-components at each time
      components_t0 = comps_t0,
      components_t  = comps_t
    )

    results[[as.character(t)]] <- result_t
  }

  # Compute delta method SEs if vcov is available
  se_list <- NULL
  cross_se <- NULL
  lambda_tests <- NULL
  beta_tests <- NULL
  if (!is.null(params$vcov)) {
    se_list <- delta_method_se(params, type = "longit", order = order, ref = ref)
    beta_scale <- if (!is.null(params$transform)) params$transform else "identity"

    # Per-time cross-sectional SEs and Wald tests
    cross_se <- list()
    lambda_tests <- list()
    beta_tests <- list()
    for (t in params$times) {
      t_char <- as.character(t)
      cross_se[[t_char]] <- .delta_method_cross(params, ref = t)

      # Wald test for lambda=0 at this time
      dt <- d[d$time == t, ]
      J <- nrow(dt)
      V <- .get_vcov_for_time(params, t)
      lambda_idx <- (2 * J + 1):(3 * J)
      V_lambda <- V[lambda_idx, lambda_idx, drop = FALSE]
      lambda_tests[[t_char]] <- tryCatch({
        W_stat <- as.numeric(t(dt$lambda) %*% solve(V_lambda) %*% dt$lambda)
        p_value <- pchisq(W_stat, df = J, lower.tail = FALSE)
        list(statistic = W_stat, df = J, p_value = p_value)
      }, error = function(e) NULL)

      # Wald test for beta homogeneous across groups at this time (scale-aware)
      beta_tests[[t_char]] <- tryCatch({
        if (J < 2L) NULL else {
          V_beta <- V[1:J, 1:J, drop = FALSE]
          Cmat   <- cbind(-1, diag(J - 1L))
          Cb     <- as.numeric(Cmat %*% dt$beta)
          W_stat <- as.numeric(t(Cb) %*% solve(Cmat %*% V_beta %*% t(Cmat)) %*% Cb)
          list(statistic = W_stat, df = J - 1L,
               p_value = pchisq(W_stat, df = J - 1L, lower.tail = FALSE),
               scale = beta_scale)
        }
      }, error = function(e) NULL)
    }
  }

  structure(
    list(
      results = results,
      se = se_list,
      cross_se = cross_se,
      lambda_tests = lambda_tests,
      beta_tests = beta_tests,
      order = order,
      ystat = ystat,
      ref = ref,
      params = params
    ),
    class = "ineqx_causal_longit"
  )
}


# ---------------------------------------------------------------------------- #
# Internal: Unified 5-parameter sequential decomposition
# ---------------------------------------------------------------------------- #

#' Unified sequential decomposition switching all 5 parameters and tracking
#' contributions to both \eqn{\tau_B} and \eqn{\tau_W} at each step.
#'
#' Walks through the parameter sequence in \code{order_5}, switching each
#' parameter from its t0 value to its t value and recording the resulting
#' change in both the cross-sectional between-group treatment effect
#' \eqn{\tau_B} and the cross-sectional within-group treatment effect
#' \eqn{\tau_W}. Telescopes exactly: summing all (B, W) contributions
#' recovers \eqn{(\tau_B(t) - \tau_B(t_0))} and \eqn{(\tau_W(t) - \tau_W(t_0))}.
#'
#' @param d0 Reference period parameters (data.frame with pi, mu0, sigma0, beta, lambda)
#' @param dt Target period parameters (same columns)
#' @param order_5 5-permutation of \code{c("beta", "lambda", "pi", "mu", "sigma")}
#' @param ystat \code{"Var"} or \code{"CV2"}
#' @return Named list. \code{parts[[p]]$B} and \code{parts[[p]]$W} for each
#'   parameter \code{p in c("beta","lambda","pi","mu","sigma")}, plus
#'   \code{parts$totals} with \code{delta_B}, \code{delta_W}, \code{delta_T}.
#' @keywords internal
.decompose_sequential <- function(d0, dt, order_5, ystat) {

  state <- list(
    pi    = d0$pi,
    mu    = d0$mu0,
    beta  = d0$beta,
    sigma = d0$sigma0,
    lambda = d0$lambda
  )

  eval_tau <- function(s) {
    wb <- compute_delta_WB(s$pi, s$mu, s$sigma, s$beta, s$lambda, ystat)
    list(tau_B = wb$tau_B, tau_W = wb$tau_W)
  }

  baseline <- eval_tau(state)
  prev_B <- baseline$tau_B
  prev_W <- baseline$tau_W

  to_map <- list(
    beta   = dt$beta,
    lambda = dt$lambda,
    pi     = dt$pi,
    mu     = dt$mu0,
    sigma  = dt$sigma0
  )

  parts <- list()
  for (p in order_5) {
    state[[p]] <- to_map[[p]]
    new_tau <- eval_tau(state)
    parts[[p]] <- list(
      B = new_tau$tau_B - prev_B,
      W = new_tau$tau_W - prev_W
    )
    prev_B <- new_tau$tau_B
    prev_W <- new_tau$tau_W
  }

  parts$totals <- list(
    delta_B = prev_B - baseline$tau_B,
    delta_W = prev_W - baseline$tau_W,
    delta_T = (prev_B + prev_W) - (baseline$tau_B + baseline$tau_W)
  )
  parts
}


# ---------------------------------------------------------------------------- #
# Internal: Expand a 3-level meta-ordering into the 5-parameter sub-ordering
# ---------------------------------------------------------------------------- #

#' Expand a meta-ordering into the underlying 5-parameter sequence.
#'
#' Mapping: \code{behavioral -> (beta, lambda)},
#' \code{compositional -> (pi)}, \code{pretreatment -> (mu, sigma)}.
#' \eqn{\beta} is switched before \eqn{\lambda} within the behavioral
#' meta-level, and \eqn{\mu} before \eqn{\sigma} within the pretreatment
#' meta-level. Users seeking full ordering invariance should use
#' \code{order = "shapley"}.
#'
#' @param order Length-3 character vector, a permutation of
#'   \code{c("behavioral", "compositional", "pretreatment")}.
#' @return Length-5 character vector, a permutation of
#'   \code{c("beta", "lambda", "pi", "mu", "sigma")}.
#' @keywords internal
.expand_order_to_5 <- function(order) {
  meta_map <- list(
    behavioral    = c("beta", "lambda"),
    compositional = "pi",
    pretreatment  = c("mu", "sigma")
  )
  unlist(lapply(order, function(m) meta_map[[m]]), use.names = FALSE)
}


# ---------------------------------------------------------------------------- #
# Internal: Compute exact cross-sectional sub-components of delta_B and delta_W
# ---------------------------------------------------------------------------- #

#' Compute exact sub-components of the cross-sectional treatment effect
#'
#' For Var: delta_B = Var_pi(beta) + 2*Cov_pi(mu0, beta)  (exact identity)
#'          delta_W = mean_pi(sigma0^2)*mean_pi(f) + Cov_pi(sigma0^2, f)  (exact)
#'
#' @param pi Group shares
#' @param mu0 Baseline group means
#' @param sigma0 Baseline group SDs
#' @param beta Treatment effects on mean
#' @param lambda Treatment effects on log-SD
#' @param ystat "Var" or "CV2"
#' @return Named list of sub-components
#' @keywords internal
.compute_cross_components <- function(pi, mu0, sigma0, beta, lambda, ystat) {

  mu1     <- mu0 + beta
  sigma1  <- sigma0 * exp(lambda)
  f       <- exp(2 * lambda) - 1

  # Population means and weighted moments
  beta_bar     <- sum(pi * beta)
  mu0_bar      <- sum(pi * mu0)
  mu1_bar      <- sum(pi * mu1)

  # V-style numerator pieces (and intermediates exposed for prints/tests)
  Var_pi_beta     <- sum(pi * (beta - beta_bar)^2)
  Cov_pi_mu_beta  <- sum(pi * (mu0 - mu0_bar) * (beta - beta_bar))
  sigma2_0_bar    <- sum(pi * sigma0^2)
  f_bar           <- sum(pi * f)
  Cov_pi_sigma2_f <- sum(pi * (sigma0^2 - sigma2_0_bar) * (f - f_bar))

  # Pre-existing inequality (used as the numerator in the rescaling term)
  B_0 <- sum(pi * (mu0 - mu0_bar)^2)   # Var_pi(mu0)
  W_0 <- sigma2_0_bar

  # CV2 denominators (1 for Var, by construction). The asymmetric three-piece
  # decomposition evaluates het/cov at the post-treatment denominator d_1.
  if (ystat == "Var") {
    d_0 <- 1
    d_1 <- 1
  } else {
    d_0 <- 1 / mu0_bar^2
    d_1 <- 1 / mu1_bar^2
  }

  # Three-piece decomposition. Under V (d_0 = d_1 = 1) the rescale terms are
  # zero by construction, and het/cov collapse to today's V definitions.
  het_B     <- d_1 * Var_pi_beta
  cov_B     <- d_1 * 2 * Cov_pi_mu_beta
  rescale_B <- B_0 * (d_1 - d_0)

  het_W     <- d_1 * W_0 * f_bar
  cov_W     <- d_1 * Cov_pi_sigma2_f
  rescale_W <- W_0 * (d_1 - d_0)

  list(
    Var_pi_beta     = Var_pi_beta,
    Cov_pi_mu_beta  = Cov_pi_mu_beta,
    mean_sigma2_0   = sigma2_0_bar,
    mean_f          = f_bar,
    Cov_pi_sigma2_f = Cov_pi_sigma2_f,
    het_B = het_B, cov_B = cov_B, rescale_B = rescale_B,
    het_W = het_W, cov_W = cov_W, rescale_W = rescale_W
  )
}
