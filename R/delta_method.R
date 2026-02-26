# ============================================================================ #
# Delta method standard errors for causal variance decomposition
# ============================================================================ #

# ---------------------------------------------------------------------------- #
# Main SE computation function
# ---------------------------------------------------------------------------- #

#' Compute delta method standard errors for causal decomposition
#'
#' Given an \code{ineqx_params} object with a variance-covariance matrix,
#' computes standard errors for the decomposition quantities using the
#' delta method: \eqn{Var(g(\hat{\theta})) \approx \nabla g(\theta)^T \Sigma_\theta \nabla g(\theta)}.
#'
#' @param params An \code{ineqx_params} object with a non-NULL \code{vcov} field
#' @param type Character, \code{"cross"} for cross-sectional or \code{"longit"}
#'   for longitudinal decomposition
#' @param order Character vector (for longitudinal), the decomposition ordering
#'
#' @return A list of standard errors for the decomposition quantities
#'
#' @keywords internal
delta_method_se <- function(params, type = c("cross", "longit"), order = NULL) {

  type <- match.arg(type)
  stopifnot(inherits(params, "ineqx_params"))

  if (is.null(params$vcov)) {
    stop("params$vcov is NULL. Provide a variance-covariance matrix to compute SEs.")
  }

  if (type == "cross") {
    .delta_method_cross(params)
  } else {
    if (is.null(order)) {
      order <- c("behavioral", "compositional", "pretreatment")
    }
    .delta_method_longit(params, order)
  }
}


# ---------------------------------------------------------------------------- #
# Cross-sectional delta method SEs
# ---------------------------------------------------------------------------- #

.delta_method_cross <- function(params) {

  d <- params$data
  ystat <- params$ystat

  # If longitudinal, use ref period
  if (params$type == "longitudinal") {
    t_use <- if (!is.null(params$ref)) params$ref else params$times[1]
    d <- d[d$time == t_use, ]
  }

  pi_j     <- d$pi
  mu0_j    <- d$mu0
  sigma0_j <- d$sigma0
  beta_j   <- d$beta
  lambda_j <- d$lambda
  J <- length(pi_j)

  # Get vcov for this time period
  V <- .get_vcov_for_time(params, if (params$type == "longitudinal") params$ref else NULL)

  # Gradients for delta_B (w.r.t. beta, mu0) — 2J vector
  grad_B <- .grad_delta_B(pi_j, mu0_j, beta_j, ystat)

  # Gradients for delta_W (w.r.t. lambda, log_sigma0) — 2J vector
  grad_W <- .grad_delta_W(pi_j, mu0_j, sigma0_j, beta_j, lambda_j, ystat)

  # Full gradient for delta_total — 4J vector
  grad_total <- c(grad_B, grad_W)

  # Extract block vcov matrices
  # V is ordered as: (β_1..J, μ0_1..J, λ_1..J, log_σ0_1..J)
  idx_B <- 1:(2 * J)           # beta and mu0 indices
  idx_W <- (2 * J + 1):(4 * J) # lambda and log_sigma0 indices

  V_B <- V[idx_B, idx_B, drop = FALSE]
  V_W <- V[idx_W, idx_W, drop = FALSE]

  se_B <- .compute_se(grad_B, V_B)
  se_W <- .compute_se(grad_W, V_W)
  se_total <- .compute_se(grad_total, V)

  list(
    se_delta_B = se_B,
    se_delta_W = se_W,
    se_delta_total = se_total
  )
}


# ---------------------------------------------------------------------------- #
# Longitudinal delta method SEs
# ---------------------------------------------------------------------------- #

.delta_method_longit <- function(params, order) {

  ref_time <- params$ref
  ystat <- params$ystat
  d <- params$data
  d0 <- d[d$time == ref_time, ]
  J <- nrow(d0)

  results <- list()

  for (t in setdiff(params$times, ref_time)) {
    dt <- d[d$time == t, ]

    # Get vcov matrices for t0 and t
    V_t0 <- .get_vcov_for_time(params, ref_time)
    V_t  <- .get_vcov_for_time(params, t)

    # Compute gradients for 6 components
    grads <- .grad_longit_components(d0, dt, order, ystat)

    # For each component, SE = sqrt(∇_t0' Σ_t0 ∇_t0 + ∇_t' Σ_t ∇_t)
    # (independent time periods for repeated cross-sections)
    ses <- list()
    for (comp_name in names(grads)) {
      g <- grads[[comp_name]]
      var_t0 <- as.numeric(t(g$grad_t0) %*% V_t0 %*% g$grad_t0)
      var_t  <- as.numeric(t(g$grad_t) %*% V_t %*% g$grad_t)
      ses[[paste0("se_", comp_name)]] <- sqrt(max(var_t0 + var_t, 0))
    }

    results[[as.character(t)]] <- ses
  }

  results
}


# ---------------------------------------------------------------------------- #
# Gradient: delta_B w.r.t. (beta, mu0) — 2J vector
# ---------------------------------------------------------------------------- #

.grad_delta_B <- function(pi, mu0, beta, ystat) {
  J <- length(pi)

  if (ystat == "Var") {
    mu1 <- mu0 + beta
    mu1_bar <- sum(pi * mu1)
    beta_bar <- sum(pi * beta)

    # ∂δ_B/∂β_k = 2·π_k·[(μ0_k + β_k) - μ̄₁]
    grad_beta <- 2 * pi * (mu1 - mu1_bar)
    # ∂δ_B/∂μ0_k = 2·π_k·[β_k - β̄]
    grad_mu0 <- 2 * pi * (beta - beta_bar)

    c(grad_beta, grad_mu0)
  } else {
    # Numerical gradient for CV2
    .numerical_grad(function(theta) {
      b <- theta[1:J]
      m <- theta[(J + 1):(2 * J)]
      compute_delta_B(pi, m, b, "CV2")
    }, c(beta, mu0))
  }
}


# ---------------------------------------------------------------------------- #
# Gradient: delta_W w.r.t. (lambda, log_sigma0) — 2J vector
# ---------------------------------------------------------------------------- #

.grad_delta_W <- function(pi, mu0, sigma0, beta, lambda, ystat) {
  J <- length(pi)

  if (ystat == "Var") {
    sigma0_sq <- sigma0^2
    f <- exp(2 * lambda) - 1

    # ∂δ_W/∂λ_k = 2·π_k·σ²_k(0)·exp(2λ_k)
    grad_lambda <- 2 * pi * sigma0_sq * exp(2 * lambda)
    # ∂δ_W/∂log(σ0_k) = 2·π_k·σ²_k(0)·f_k
    # (chain rule: σ0 = exp(log_σ0), so ∂W/∂log_σ0 = σ0·∂W/∂σ0 = 2·π·σ0²·f)
    grad_log_sigma0 <- 2 * pi * sigma0_sq * f

    c(grad_lambda, grad_log_sigma0)
  } else {
    # Numerical gradient for CV2
    # For CV2, delta_W depends on mu0 and beta (via grand means),
    # but we treat those as fixed when differentiating delta_W w.r.t. (lambda, log_sigma0)
    .numerical_grad(function(theta) {
      l <- theta[1:J]
      ls <- theta[(J + 1):(2 * J)]
      s0 <- exp(ls)
      s1 <- s0 * exp(l)
      mu1 <- mu0 + beta
      gmu1 <- sum(pi * mu1)
      gmu0 <- sum(pi * mu0)
      # delta_W = W(d=1) - W(d=0) for CV2
      sum(pi * s1^2) / gmu1^2 - sum(pi * s0^2) / gmu0^2
    }, c(lambda, log(sigma0)))
  }
}


# ---------------------------------------------------------------------------- #
# Gradients for longitudinal components (sequential switching)
# ---------------------------------------------------------------------------- #

.grad_longit_components <- function(d0, dt, order, ystat) {

  J <- nrow(d0)

  # Compute gradients of each component w.r.t. (θ_t0, θ_t)
  # Each θ = (β_1..J, μ0_1..J, λ_1..J, log_σ0_1..J) — 4J vector

  # Between-group gradients
  between_grads <- .grad_between_sequential(d0, dt, order, ystat)

  # Within-group gradients
  within_grads <- .grad_within_sequential(d0, dt, order, ystat)

  # Combine into named list with grad_t0 and grad_t for each component
  result <- list()

  # 6 detailed components
  result$Delta_beta   <- between_grads$behavioral
  result$Delta_lambda <- within_grads$behavioral
  result$Delta_pi_B   <- between_grads$compositional
  result$Delta_pi_W   <- within_grads$compositional
  result$Delta_mu     <- between_grads$pretreatment
  result$Delta_sigma  <- within_grads$pretreatment

  # 3 combined components
  result$Delta_behavioral <- list(
    grad_t0 = result$Delta_beta$grad_t0 + result$Delta_lambda$grad_t0,
    grad_t  = result$Delta_beta$grad_t  + result$Delta_lambda$grad_t
  )
  result$Delta_compositional <- list(
    grad_t0 = result$Delta_pi_B$grad_t0 + result$Delta_pi_W$grad_t0,
    grad_t  = result$Delta_pi_B$grad_t  + result$Delta_pi_W$grad_t
  )
  result$Delta_pretreatment <- list(
    grad_t0 = result$Delta_mu$grad_t0 + result$Delta_sigma$grad_t0,
    grad_t  = result$Delta_mu$grad_t  + result$Delta_sigma$grad_t
  )

  # Between and within totals
  result$Delta_B <- list(
    grad_t0 = between_grads$total$grad_t0,
    grad_t  = between_grads$total$grad_t
  )
  result$Delta_W <- list(
    grad_t0 = within_grads$total$grad_t0,
    grad_t  = within_grads$total$grad_t
  )

  # Grand total
  result$Delta_total <- list(
    grad_t0 = between_grads$total$grad_t0 + within_grads$total$grad_t0,
    grad_t  = between_grads$total$grad_t  + within_grads$total$grad_t
  )

  result
}


# ---------------------------------------------------------------------------- #
# Gradient of between-group sequential decomposition
# Returns gradients w.r.t. full 4J theta for t0 and t
# ---------------------------------------------------------------------------- #

.grad_between_sequential <- function(d0, dt, order, ystat) {

  J <- nrow(d0)

  # Map component names to between-group parameters
  # current state: beta, pi (fixed), mu0
  # For between: we switch beta, pi, mu0
  # Since pi is fixed, its "gradient" is zero

  # We need to track which parameters are currently at t0 vs t values
  # and compute the gradient of delta_B at each switching step

  param_map <- list(
    behavioral    = list(param = "beta", from = d0$beta,  to = dt$beta),
    compositional = list(param = "pi",   from = d0$pi,    to = dt$pi),
    pretreatment  = list(param = "mu",   from = d0$mu0,   to = dt$mu0)
  )

  # Start with all at t0
  current <- list(beta = d0$beta, pi = d0$pi, mu = d0$mu0)

  # Helper: gradient of delta_B(current params) w.r.t. (beta, mu0, pi)
  # Returns gradient w.r.t. (β, μ0) only (pi is fixed in the delta method sense)
  .eval_grad_B <- function(pars) {
    .grad_delta_B(pars$pi, pars$mu, pars$beta, ystat)
  }

  components <- list()
  prev_grad_B <- .eval_grad_B(current)

  for (step in order) {
    p <- param_map[[step]]$param
    current[[p]] <- param_map[[step]]$to
    new_grad_B <- .eval_grad_B(current)

    # The component value = delta_B(after) - delta_B(before)
    # So grad(component) = grad(delta_B(after)) - grad(delta_B(before))
    # But we need to express this w.r.t. the full theta at t0 and t

    # The gradient of this step w.r.t. the actual parameters is complex
    # because which parameters come from t0 vs t changes at each step.
    # Use numerical differentiation for robustness.

    components[[step]] <- .grad_step_numerical(d0, dt, order, step, "between", ystat)
    prev_grad_B <- new_grad_B
  }

  # Total gradient
  total_grad_t0 <- rep(0, 4 * J)
  total_grad_t  <- rep(0, 4 * J)
  for (step in order) {
    total_grad_t0 <- total_grad_t0 + components[[step]]$grad_t0
    total_grad_t  <- total_grad_t  + components[[step]]$grad_t
  }
  components$total <- list(grad_t0 = total_grad_t0, grad_t = total_grad_t)

  components
}


.grad_within_sequential <- function(d0, dt, order, ystat) {

  J <- nrow(d0)

  components <- list()
  for (step in order) {
    components[[step]] <- .grad_step_numerical(d0, dt, order, step, "within", ystat)
  }

  total_grad_t0 <- rep(0, 4 * J)
  total_grad_t  <- rep(0, 4 * J)
  for (step in order) {
    total_grad_t0 <- total_grad_t0 + components[[step]]$grad_t0
    total_grad_t  <- total_grad_t  + components[[step]]$grad_t
  }
  components$total <- list(grad_t0 = total_grad_t0, grad_t = total_grad_t)

  components
}


# ---------------------------------------------------------------------------- #
# Numerical gradient of a single sequential switching step
# w.r.t. full theta vectors at t0 and t
# ---------------------------------------------------------------------------- #

.grad_step_numerical <- function(d0, dt, order, target_step, sub_decomp, ystat) {

  J <- nrow(d0)

  # theta_t0 = (β_1..J, μ0_1..J, λ_1..J, log_σ0_1..J)
  theta_t0 <- c(d0$beta, d0$mu0, d0$lambda, log(d0$sigma0))
  theta_t  <- c(dt$beta, dt$mu0, dt$lambda, log(dt$sigma0))

  # Function that computes the value of the target step
  # given theta_t0 and theta_t
  f_step <- function(th_t0, th_t) {
    # Reconstruct d0 and dt from theta vectors
    d0_mod <- data.frame(
      group = d0$group,
      pi = d0$pi,
      beta   = th_t0[1:J],
      mu0    = th_t0[(J + 1):(2 * J)],
      lambda = th_t0[(2 * J + 1):(3 * J)],
      sigma0 = exp(th_t0[(3 * J + 1):(4 * J)]),
      stringsAsFactors = FALSE
    )
    dt_mod <- data.frame(
      group = dt$group,
      pi = dt$pi,
      beta   = th_t[(1:J)],
      mu0    = th_t[(J + 1):(2 * J)],
      lambda = th_t[(2 * J + 1):(3 * J)],
      sigma0 = exp(th_t[(3 * J + 1):(4 * J)]),
      stringsAsFactors = FALSE
    )

    if (sub_decomp == "between") {
      comps <- .decompose_between_sequential(d0_mod, dt_mod, order, ystat)
    } else {
      comps <- .decompose_within_sequential(d0_mod, dt_mod, order, ystat)
    }
    comps[[target_step]]
  }

  # Numerical gradient w.r.t. theta_t0
  grad_t0 <- .numerical_grad(function(th) f_step(th, theta_t), theta_t0)

  # Numerical gradient w.r.t. theta_t
  grad_t <- .numerical_grad(function(th) f_step(theta_t0, th), theta_t)

  list(grad_t0 = grad_t0, grad_t = grad_t)
}


# ---------------------------------------------------------------------------- #
# Helper: Get vcov matrix for a specific time period
# ---------------------------------------------------------------------------- #

.get_vcov_for_time <- function(params, time) {
  V <- params$vcov

  if (is.list(V) && !is.matrix(V)) {
    # Named list of per-period vcov matrices
    key <- as.character(time)
    if (!key %in% names(V)) {
      stop("No vcov matrix found for time period ", time,
           ". Available: ", paste(names(V), collapse = ", "))
    }
    return(V[[key]])
  }

  # Single matrix — return as-is
  V
}


# ---------------------------------------------------------------------------- #
# Helper: Compute SE from gradient and vcov
# ---------------------------------------------------------------------------- #

.compute_se <- function(grad, vcov) {
  variance <- as.numeric(t(grad) %*% vcov %*% grad)
  sqrt(max(variance, 0))
}


# ---------------------------------------------------------------------------- #
# Helper: Numerical gradient (forward difference)
# ---------------------------------------------------------------------------- #

.numerical_grad <- function(f, x, h = 1e-7) {
  fx <- f(x)
  vapply(seq_along(x), function(i) {
    x_plus <- x
    x_plus[i] <- x[i] + h
    (f(x_plus) - fx) / h
  }, numeric(1))
}


# ---------------------------------------------------------------------------- #
# Helper: Significance stars
# ---------------------------------------------------------------------------- #

.signif_stars <- function(pval) {
  if (is.na(pval)) return("")
  if (pval < 0.001) return("***")
  if (pval < 0.01)  return("**")
  if (pval < 0.05)  return("*")
  if (pval < 0.1)   return(".")
  ""
}

.pval_from_se <- function(estimate, se) {
  if (is.na(se) || se <= 0) return(NA_real_)
  2 * stats::pnorm(-abs(estimate / se))
}
