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
delta_method_se <- function(params, type = c("cross", "longit"), order = NULL,
                            ref = NULL) {

  type <- match.arg(type)
  stopifnot(inherits(params, "ineqx_params"))

  if (is.null(params$vcov)) {
    stop("params$vcov is NULL. Provide a variance-covariance matrix to compute SEs.")
  }

  if (type == "cross") {
    .delta_method_cross(params, ref = ref)
  } else {
    if (is.null(order)) {
      order <- c("behavioral", "compositional", "pretreatment")
    }
    .delta_method_longit(params, order, ref = ref)
  }
}


# ---------------------------------------------------------------------------- #
# Cross-sectional delta method SEs
# ---------------------------------------------------------------------------- #

.delta_method_cross <- function(params, ref = NULL) {

  d <- params$data
  ystat <- params$ystat

  # If longitudinal, use ref period
  if (params$type == "longitudinal") {
    t_use <- if (!is.null(ref)) ref else params$times[1]
    d <- d[d$time == t_use, ]
  }

  pi_j     <- d$pi
  mu0_j    <- d$mu0
  sigma0_j <- d$sigma0
  beta_j   <- d$beta
  lambda_j <- d$lambda
  J <- length(pi_j)

  # Get vcov for this time period
  V <- .get_vcov_for_time(params, if (params$type == "longitudinal") ref else NULL)

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

  # Sub-component SEs
  sub_grads <- .grad_cross_subcomponents(pi_j, mu0_j, sigma0_j, beta_j, lambda_j, ystat)

  se_sub <- list(
    se_het_B     = .compute_se(sub_grads$het_B,     V_B),
    se_cov_B     = .compute_se(sub_grads$cov_B,     V_B),
    se_rescale_B = .compute_se(sub_grads$rescale_B, V_B),
    se_het_W     = .compute_se(sub_grads$het_W,     V_W),
    se_cov_W     = .compute_se(sub_grads$cov_W,     V_W),
    se_rescale_W = .compute_se(sub_grads$rescale_W, V_W)
  )

  # SEs for inequality levels: ineq[Y|T=0] and ineq[Y|T=1]
  grad_ineq0 <- .grad_ineq_level(pi_j, mu0_j, sigma0_j, beta_j, lambda_j, ystat, treated = FALSE)
  grad_ineq1 <- .grad_ineq_level(pi_j, mu0_j, sigma0_j, beta_j, lambda_j, ystat, treated = TRUE)
  se_ineq_control <- .compute_se(grad_ineq0, V)
  se_ineq_treated <- .compute_se(grad_ineq1, V)

  list(
    se_tau_B = se_B,
    se_tau_W = se_W,
    se_tau_total = se_total,
    se_ineq_control = se_ineq_control,
    se_ineq_treated = se_ineq_treated,
    se_sub = se_sub
  )
}


# ---------------------------------------------------------------------------- #
# Longitudinal delta method SEs
# ---------------------------------------------------------------------------- #

.delta_method_longit <- function(params, order, ref = NULL) {

  if (is.null(ref)) {
    stop("'ref' is required for longitudinal delta method SEs.")
  }
  ref_time <- ref
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
# Gradients for cross-sectional sub-components
# Returns gradients for het_B, cov_B (w.r.t. beta, mu0 — 2J)
#                    and het_W, cov_W (w.r.t. lambda, log_sigma0 — 2J)
# ---------------------------------------------------------------------------- #

.grad_cross_subcomponents <- function(pi, mu0, sigma0, beta, lambda, ystat) {
  J <- length(pi)

  if (ystat == "Var") {
    # Under V (d_0 = d_1 = 1), the new het/cov definitions reduce to today's
    # closed forms and rescale_B/W are identically zero. Analytical gradients
    # below match the new het_B = Var_pi(beta), cov_B = 2*Cov_pi(mu0, beta),
    # het_W = W_0*f_bar, cov_W = Cov_pi(sigma0^2, f).

    # --- Between-group sub-components (w.r.t. beta, mu0) ---
    beta_bar <- sum(pi * beta)
    grad_het_B_beta <- 2 * pi * (beta - beta_bar)
    grad_het_B_mu0  <- rep(0, J)
    grad_het_B <- c(grad_het_B_beta, grad_het_B_mu0)

    mu0_bar <- sum(pi * mu0)
    grad_cov_B_beta <- 2 * pi * (mu0 - mu0_bar)
    grad_cov_B_mu0  <- 2 * pi * (beta - beta_bar)
    grad_cov_B <- c(grad_cov_B_beta, grad_cov_B_mu0)

    # rescale_B = B_0 * (d_1 - d_0) = 0 under V
    grad_rescale_B <- rep(0, 2 * J)

    # --- Within-group sub-components (w.r.t. lambda, log_sigma0) ---
    sigma0_sq  <- sigma0^2
    f          <- exp(2 * lambda) - 1
    sigma2_bar <- sum(pi * sigma0_sq)
    f_bar      <- sum(pi * f)

    grad_het_W_lambda <- sigma2_bar * pi * 2 * exp(2 * lambda)
    grad_het_W_lsig   <- f_bar * pi * 2 * sigma0_sq
    grad_het_W <- c(grad_het_W_lambda, grad_het_W_lsig)

    grad_cov_W_lambda <- pi * (sigma0_sq - sigma2_bar) * 2 * exp(2 * lambda)
    grad_cov_W_lsig   <- pi * 2 * sigma0_sq * (f - f_bar)
    grad_cov_W <- c(grad_cov_W_lambda, grad_cov_W_lsig)

    # rescale_W = W_0 * (d_1 - d_0) = 0 under V
    grad_rescale_W <- rep(0, 2 * J)

  } else {
    # CV2: numerical gradients. Each W-component closure treats d_1, d_0
    # (functions of B-block parameters) as constants from the outer scope —
    # the same 2J-against-W approximation used previously for cov_W (now
    # rescale_W). Same goes for B-component closures w.r.t. log_sigma0.
    mu0_bar <- sum(pi * mu0)
    mu1     <- mu0 + beta
    mu1_bar <- sum(pi * mu1)
    d_0     <- 1 / mu0_bar^2
    d_1     <- 1 / mu1_bar^2

    # het_B = d_1 * Var_pi(beta), cov_B = d_1 * 2 * Cov_pi(mu0, beta),
    # rescale_B = Var_pi(mu0) * (d_1 - d_0). All three depend on (beta, mu0).
    grad_het_B <- .numerical_grad(function(theta) {
      b <- theta[1:J]; m <- theta[(J + 1):(2 * J)]
      bbar  <- sum(pi * b)
      mbar1 <- sum(pi * (m + b))
      sum(pi * (b - bbar)^2) / mbar1^2
    }, c(beta, mu0))

    grad_cov_B <- .numerical_grad(function(theta) {
      b <- theta[1:J]; m <- theta[(J + 1):(2 * J)]
      bbar  <- sum(pi * b)
      mbar0 <- sum(pi * m)
      mbar1 <- sum(pi * (m + b))
      2 * sum(pi * (m - mbar0) * (b - bbar)) / mbar1^2
    }, c(beta, mu0))

    grad_rescale_B <- .numerical_grad(function(theta) {
      b <- theta[1:J]; m <- theta[(J + 1):(2 * J)]
      mbar0 <- sum(pi * m); mbar1 <- sum(pi * (m + b))
      Vmu0  <- sum(pi * (m - mbar0)^2)
      Vmu0 * (1 / mbar1^2 - 1 / mbar0^2)
    }, c(beta, mu0))

    # het_W = d_1 * W_0 * f_bar, cov_W = d_1 * Cov_pi(sigma0^2, f),
    # rescale_W = W_0 * (d_1 - d_0). W gradients hold d_0, d_1 fixed.
    grad_het_W <- .numerical_grad(function(theta) {
      l <- theta[1:J]; ls <- theta[(J + 1):(2 * J)]
      s0   <- exp(ls)
      W0   <- sum(pi * s0^2)
      fbar <- sum(pi * (exp(2 * l) - 1))
      d_1 * W0 * fbar
    }, c(lambda, log(sigma0)))

    grad_cov_W <- .numerical_grad(function(theta) {
      l <- theta[1:J]; ls <- theta[(J + 1):(2 * J)]
      s0   <- exp(ls)
      s2   <- s0^2
      fv   <- exp(2 * l) - 1
      s2bar <- sum(pi * s2)
      fbar  <- sum(pi * fv)
      d_1 * sum(pi * (s2 - s2bar) * (fv - fbar))
    }, c(lambda, log(sigma0)))

    grad_rescale_W <- .numerical_grad(function(theta) {
      l <- theta[1:J]; ls <- theta[(J + 1):(2 * J)]
      s0 <- exp(ls)
      W0 <- sum(pi * s0^2)
      W0 * (d_1 - d_0)
    }, c(lambda, log(sigma0)))
  }

  list(
    het_B = grad_het_B, cov_B = grad_cov_B, rescale_B = grad_rescale_B,
    het_W = grad_het_W, cov_W = grad_cov_W, rescale_W = grad_rescale_W
  )
}


# ---------------------------------------------------------------------------- #
# Gradient of total inequality level w.r.t. full 4J parameter vector
# (β_1..J, μ0_1..J, λ_1..J, log_σ0_1..J)
# ---------------------------------------------------------------------------- #

.grad_ineq_level <- function(pi, mu0, sigma0, beta, lambda, ystat, treated) {
  J <- length(pi)

  if (ystat == "Var") {
    if (treated) {
      # Var[Y|T=1] = Σπ(μ₁-μ̄₁)² + Σπσ₁²  where μ₁=μ₀+β, σ₁=σ₀·exp(λ)
      mu1 <- mu0 + beta
      mu1_bar <- sum(pi * mu1)
      sigma1_sq <- sigma0^2 * exp(2 * lambda)

      grad_beta    <- 2 * pi * (mu1 - mu1_bar)
      grad_mu0     <- 2 * pi * (mu1 - mu1_bar)
      grad_lambda  <- 2 * pi * sigma1_sq
      grad_lsig    <- 2 * pi * sigma1_sq
    } else {
      # Var[Y|T=0] = Σπ(μ₀-μ̄₀)² + Σπσ₀²
      mu0_bar <- sum(pi * mu0)

      grad_beta    <- rep(0, J)
      grad_mu0     <- 2 * pi * (mu0 - mu0_bar)
      grad_lambda  <- rep(0, J)
      grad_lsig    <- 2 * pi * sigma0^2
    }
    c(grad_beta, grad_mu0, grad_lambda, grad_lsig)

  } else {
    # CV2: numerical gradient w.r.t. full 4J vector
    .numerical_grad(function(theta) {
      b  <- theta[1:J]
      m  <- theta[(J + 1):(2 * J)]
      l  <- theta[(2 * J + 1):(3 * J)]
      ls <- theta[(3 * J + 1):(4 * J)]
      s0 <- exp(ls)
      if (treated) {
        CV2T_pi(pi, m + b, s0 * exp(l))
      } else {
        CV2T_pi(pi, m, s0)
      }
    }, c(beta, mu0, lambda, log(sigma0)))
  }
}


# ---------------------------------------------------------------------------- #
# Gradients for longitudinal components (sequential switching)
# ---------------------------------------------------------------------------- #

.grad_longit_components <- function(d0, dt, order, ystat) {

  J <- nrow(d0)
  order_5 <- .expand_order_to_5(order)

  # Compute gradients of each split component w.r.t. (theta_t0, theta_t).
  # theta = (beta_1..J, mu0_1..J, lambda_1..J, log_sigma0_1..J) — 4J vector.
  params <- c("beta", "lambda", "pi", "mu", "sigma")
  sides  <- c("B", "W")

  result <- list()

  # 10 split components
  for (p in params) {
    for (side in sides) {
      key <- paste0("Delta_", p, "_", side)
      result[[key]] <- .grad_unified_step(d0, dt, order_5, p, side, ystat)
    }
  }

  # Aggregate parameter components (sums of B + W).
  # pi keeps its split form; other parameters get an aggregate too.
  for (p in c("beta", "lambda", "mu", "sigma")) {
    keyB <- paste0("Delta_", p, "_B")
    keyW <- paste0("Delta_", p, "_W")
    result[[paste0("Delta_", p)]] <- list(
      grad_t0 = result[[keyB]]$grad_t0 + result[[keyW]]$grad_t0,
      grad_t  = result[[keyB]]$grad_t  + result[[keyW]]$grad_t
    )
  }

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

  # 4-component grouping
  result$Delta_pi <- list(
    grad_t0 = result$Delta_pi_B$grad_t0 + result$Delta_pi_W$grad_t0,
    grad_t  = result$Delta_pi_B$grad_t  + result$Delta_pi_W$grad_t
  )
  result$Delta_pre <- list(
    grad_t0 = result$Delta_mu$grad_t0 + result$Delta_sigma$grad_t0,
    grad_t  = result$Delta_mu$grad_t  + result$Delta_sigma$grad_t
  )

  # Between and within totals (sums over the B and W halves of all 5 params)
  total_B_t0 <- rep(0, 4 * J); total_B_t <- rep(0, 4 * J)
  total_W_t0 <- rep(0, 4 * J); total_W_t <- rep(0, 4 * J)
  for (p in params) {
    keyB <- paste0("Delta_", p, "_B")
    keyW <- paste0("Delta_", p, "_W")
    total_B_t0 <- total_B_t0 + result[[keyB]]$grad_t0
    total_B_t  <- total_B_t  + result[[keyB]]$grad_t
    total_W_t0 <- total_W_t0 + result[[keyW]]$grad_t0
    total_W_t  <- total_W_t  + result[[keyW]]$grad_t
  }
  result$Delta_B <- list(grad_t0 = total_B_t0, grad_t = total_B_t)
  result$Delta_W <- list(grad_t0 = total_W_t0, grad_t = total_W_t)

  # Grand total
  result$Delta_total <- list(
    grad_t0 = total_B_t0 + total_W_t0,
    grad_t  = total_B_t  + total_W_t
  )

  result
}


# ---------------------------------------------------------------------------- #
# Numerical gradient of a single split component (param x side)
# w.r.t. full theta vectors at t0 and t
# ---------------------------------------------------------------------------- #

.grad_unified_step <- function(d0, dt, order_5, target_param, target_side, ystat) {

  J <- nrow(d0)

  # theta_t0 = (beta_1..J, mu0_1..J, lambda_1..J, log_sigma0_1..J)
  theta_t0 <- c(d0$beta, d0$mu0, d0$lambda, log(d0$sigma0))
  theta_t  <- c(dt$beta, dt$mu0, dt$lambda, log(dt$sigma0))

  f_step <- function(th_t0, th_t) {
    d0_mod <- data.frame(
      group  = d0$group,
      pi     = d0$pi,
      beta   = th_t0[1:J],
      mu0    = th_t0[(J + 1):(2 * J)],
      lambda = th_t0[(2 * J + 1):(3 * J)],
      sigma0 = exp(th_t0[(3 * J + 1):(4 * J)]),
      stringsAsFactors = FALSE
    )
    dt_mod <- data.frame(
      group  = dt$group,
      pi     = dt$pi,
      beta   = th_t[1:J],
      mu0    = th_t[(J + 1):(2 * J)],
      lambda = th_t[(2 * J + 1):(3 * J)],
      sigma0 = exp(th_t[(3 * J + 1):(4 * J)]),
      stringsAsFactors = FALSE
    )

    parts <- .decompose_sequential(d0_mod, dt_mod, order_5, ystat)
    parts[[target_param]][[target_side]]
  }

  grad_t0 <- .numerical_grad(function(th) f_step(th, theta_t), theta_t0)
  grad_t  <- .numerical_grad(function(th) f_step(theta_t0, th), theta_t)

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


# ============================================================================ #
# Delta method standard errors for descriptive variance decomposition
# ============================================================================ #

# ---------------------------------------------------------------------------- #
# Helper: Gradient of a descriptive delta component w.r.t. theta at ref and t
# theta = (pi_1..G, mu_1..G, sigma2_1..G)
# component: "mu", "sigma", "pi", or "T" (total = sum of all three)
# ---------------------------------------------------------------------------- #

.grad_desc_delta <- function(pi_ref, mu_ref, s2_ref, pi_t, mu_t, s2_t,
                              component, ystat, order) {
  G <- length(pi_ref)
  use_shapley <- identical(order, "shapley")

  # Inequality evaluation: (pi, mu, sigma) -> list(W, B, Total)
  .eval_ineq <- function(pi, mu, sigma) {
    if (ystat == "Var") {
      W <- VarW_pi(pi, sigma)
      B <- VarB_pi(pi, mu)
    } else {
      W <- CV2W_pi(pi, mu, sigma)
      B <- CV2B_pi(pi, mu)
    }
    list(W = W, B = B, Total = W + B)
  }

  # Single-ordering decomposition (mirrors desc_decompose.R logic)
  .decompose_one_order <- function(p0, m0, s0, pt, mt, st, ord) {
    param_map <- list(
      mu    = list(from = m0, to = mt),
      sigma = list(from = s0, to = st),
      pi    = list(from = p0, to = pt)
    )
    current <- list(pi = p0, mu = m0, sigma = s0)
    prev <- .eval_ineq(current$pi, current$mu, current$sigma)

    components <- list()
    for (step in ord) {
      current[[step]] <- param_map[[step]]$to
      new_val <- .eval_ineq(current$pi, current$mu, current$sigma)
      components[[step]] <- new_val$Total - prev$Total
      prev <- new_val
    }
    components
  }

  # Compute the delta_T value for a given component from theta vectors
  # theta_ref = (pi_ref, mu_ref, sigma2_ref), theta_t = (pi_t, mu_t, sigma2_t)
  .eval_component <- function(th_ref, th_t) {
    p0 <- th_ref[1:G]; p0 <- p0 / sum(p0)  # renormalize after perturbation
    m0 <- th_ref[(G+1):(2*G)]
    s0 <- sqrt(th_ref[(2*G+1):(3*G)])       # sigma2 -> sigma
    pt <- th_t[1:G]; pt <- pt / sum(pt)
    mt <- th_t[(G+1):(2*G)]
    st <- sqrt(th_t[(2*G+1):(3*G)])

    if (use_shapley) {
      perms <- .all_permutations(c("mu", "sigma", "pi"))
      total <- 0
      for (perm in perms) {
        comp <- .decompose_one_order(p0, m0, s0, pt, mt, st, perm)
        if (component == "T") {
          total <- total + comp$mu + comp$sigma + comp$pi
        } else {
          total <- total + comp[[component]]
        }
      }
      total / length(perms)
    } else {
      comp <- .decompose_one_order(p0, m0, s0, pt, mt, st, order)
      if (component == "T") {
        comp$mu + comp$sigma + comp$pi
      } else {
        comp[[component]]
      }
    }
  }

  theta_ref <- c(pi_ref, mu_ref, s2_ref)
  theta_t   <- c(pi_t, mu_t, s2_t)

  grad_ref <- .numerical_grad(function(th) .eval_component(th, theta_t), theta_ref)
  grad_t   <- .numerical_grad(function(th) .eval_component(theta_ref, th), theta_t)

  list(grad_ref = grad_ref, grad_t = grad_t)
}


#' Compute delta method SEs for descriptive decomposition
#'
#' Given group-level statistics (pi, mu, sigma, n), constructs the sampling
#' covariance matrix and applies the delta method to compute SEs for VarW,
#' VarB, VarT, and (if applicable) counterfactual deltas.
#'
#' @param wibe Data.frame with columns: time, group, n, pi, mu, sigma, sigma2
#' @param totals Data.frame with columns: time, VarW, VarB, VarT, etc.
#' @param deltas Data.frame with delta columns, or NULL
#' @param ystat Character, "Var" or "CV2"
#' @param ref Reference time period (needed for delta SEs), or NULL
#' @param order Decomposition order: "shapley" or permutation of c("mu","sigma","pi")
#'
#' @return A list with:
#'   \item{totals}{Named list by time, each containing se_VarW, se_VarB, se_VarT
#'     (or se_CV2W, se_CV2B, se_CV2T)}
#'   \item{deltas}{Named list by time, each containing se_delta_mu, se_delta_sigma,
#'     se_delta_pi, se_delta_T (if deltas provided)}
#'
#' @keywords internal
delta_method_desc_se <- function(wibe, totals, deltas, ystat,
                                  ref = NULL, order = "shapley") {

  time_levels <- if ("time" %in% names(totals)) totals$time else 1L

  # --- SEs for totals (VarW, VarB, VarT) per time ---
  se_totals <- list()
  for (t in time_levels) {
    if ("time" %in% names(wibe)) {
      sub <- wibe[wibe$time == t, ]
    } else {
      sub <- wibe
    }

    pi_vec <- sub$pi
    mu_vec <- sub$mu
    sigma2_vec <- sub$sigma2
    n_vec <- sub$n
    G <- length(pi_vec)
    N <- sum(n_vec)

    # Sampling vcov: theta = (pi_1..G, mu_1..G, sigma2_1..G)
    V <- .desc_sampling_vcov(pi_vec, mu_vec, sigma2_vec, n_vec, N)

    if (ystat == "Var") {
      # Analytical gradients for VarW, VarB, VarT w.r.t. (pi, mu, sigma2)
      grad_W <- .grad_desc_VarW(pi_vec, mu_vec, sigma2_vec)
      grad_B <- .grad_desc_VarB(pi_vec, mu_vec, sigma2_vec)
      grad_T <- grad_W + grad_B

      se_totals[[as.character(t)]] <- list(
        se_VarW = .compute_se(grad_W, V),
        se_VarB = .compute_se(grad_B, V),
        se_VarT = .compute_se(grad_T, V)
      )
    } else {
      # Numerical gradients for CV2
      grad_W <- .numerical_grad(function(th) {
        p <- th[1:G]; m <- th[(G+1):(2*G)]; s2 <- th[(2*G+1):(3*G)]
        CV2W_pi(p, m, sqrt(s2))
      }, c(pi_vec, mu_vec, sigma2_vec))
      grad_B <- .numerical_grad(function(th) {
        p <- th[1:G]; m <- th[(G+1):(2*G)]
        CV2B_pi(p, m)
      }, c(pi_vec, mu_vec, sigma2_vec))
      grad_T <- grad_W + grad_B

      se_totals[[as.character(t)]] <- list(
        se_CV2W = .compute_se(grad_W, V),
        se_CV2B = .compute_se(grad_B, V),
        se_CV2T = .compute_se(grad_T, V)
      )
    }
  }

  # --- SEs for deltas (if available) ---
  se_deltas <- NULL
  if (!is.null(deltas) && !is.null(ref)) {

    # Get reference-period statistics
    if ("time" %in% names(wibe)) {
      ref_sub <- wibe[wibe$time == ref, ]
    } else {
      ref_sub <- wibe
    }
    pi_ref  <- ref_sub$pi
    mu_ref  <- ref_sub$mu
    s2_ref  <- ref_sub$sigma2
    n_ref   <- ref_sub$n
    G       <- length(pi_ref)
    N_ref   <- sum(n_ref)
    V_ref   <- .desc_sampling_vcov(pi_ref, mu_ref, s2_ref, n_ref, N_ref)

    se_deltas <- list()
    for (i in seq_len(nrow(deltas))) {
      t <- deltas$time[i]

      # At the reference period, deltas are 0 by definition → SE = 0
      if (t == ref) {
        se_deltas[[as.character(t)]] <- list(
          se_delta_mu    = 0,
          se_delta_sigma = 0,
          se_delta_pi    = 0,
          se_delta_T     = 0
        )
        next
      }

      # Get target-period statistics
      t_sub <- wibe[wibe$time == t, ]
      pi_t  <- t_sub$pi
      mu_t  <- t_sub$mu
      s2_t  <- t_sub$sigma2
      n_t   <- t_sub$n
      N_t   <- sum(n_t)
      V_t   <- .desc_sampling_vcov(pi_t, mu_t, s2_t, n_t, N_t)

      se_t <- list()
      # Compute SE for each component + total
      for (comp in c("mu", "sigma", "pi", "T")) {
        grads <- .grad_desc_delta(
          pi_ref, mu_ref, s2_ref, pi_t, mu_t, s2_t,
          component = comp, ystat = ystat, order = order
        )
        var_ref <- as.numeric(t(grads$grad_ref) %*% V_ref %*% grads$grad_ref)
        var_t   <- as.numeric(t(grads$grad_t)   %*% V_t   %*% grads$grad_t)
        se_t[[paste0("se_delta_", comp)]] <- sqrt(max(var_ref + var_t, 0))
      }
      se_deltas[[as.character(t)]] <- se_t
    }
  } else if (!is.null(deltas)) {
    # No ref available — cannot compute delta SEs
    se_deltas <- list()
    for (i in seq_len(nrow(deltas))) {
      t <- deltas$time[i]
      se_deltas[[as.character(t)]] <- list(
        se_delta_mu    = NA_real_,
        se_delta_sigma = NA_real_,
        se_delta_pi    = NA_real_,
        se_delta_T     = NA_real_
      )
    }
  }

  list(totals = se_totals, deltas = se_deltas)
}


# ---------------------------------------------------------------------------- #
# Sampling vcov for descriptive group statistics
# theta = (pi_1..G, mu_1..G, sigma2_1..G), a 3G-vector
# ---------------------------------------------------------------------------- #

.desc_sampling_vcov <- function(pi, mu, sigma2, n, N) {
  G <- length(pi)
  d <- 3 * G
  V <- matrix(0, nrow = d, ncol = d)

  # Block 1: Var(pi) — multinomial covariance (G x G)
  for (g in seq_len(G)) {
    for (h in seq_len(G)) {
      if (g == h) {
        V[g, h] <- pi[g] * (1 - pi[g]) / N
      } else {
        V[g, h] <- -pi[g] * pi[h] / N
      }
    }
  }

  # Block 2: Var(mu) — independent across groups (G x G diagonal)
  for (g in seq_len(G)) {
    idx <- G + g
    if (n[g] > 0) {
      V[idx, idx] <- sigma2[g] / n[g]
    }
  }

  # Block 3: Var(sigma2) — under normality (G x G diagonal)
  for (g in seq_len(G)) {
    idx <- 2 * G + g
    if (n[g] > 1) {
      V[idx, idx] <- 2 * sigma2[g]^2 / (n[g] - 1)
    }
  }

  V
}


# ---------------------------------------------------------------------------- #
# Analytical gradient of VarW w.r.t. theta = (pi, mu, sigma2)
# VarW = sum(pi_g * sigma2_g)
# ---------------------------------------------------------------------------- #

.grad_desc_VarW <- function(pi, mu, sigma2) {
  G <- length(pi)
  grad <- numeric(3 * G)

  # d VarW / d pi_g = sigma2_g
  grad[1:G] <- sigma2

  # d VarW / d mu_g = 0
  # (already zero)

  # d VarW / d sigma2_g = pi_g
  grad[(2*G+1):(3*G)] <- pi

  grad
}


# ---------------------------------------------------------------------------- #
# Analytical gradient of VarB w.r.t. theta = (pi, mu, sigma2)
# VarB = sum(pi_g * (mu_g - mu_bar)^2), where mu_bar = sum(pi_g * mu_g)
# ---------------------------------------------------------------------------- #

.grad_desc_VarB <- function(pi, mu, sigma2) {
  G <- length(pi)
  grad <- numeric(3 * G)

  mu_bar <- sum(pi * mu)
  dev <- mu - mu_bar

  # d VarB / d pi_g = (mu_g - mu_bar)^2 - 2 * mu_g * sum(pi_h * (mu_h - mu_bar))
  # Simplification: d mu_bar / d pi_g = mu_g
  # d VarB / d pi_g = dev_g^2 - 2 * sum(pi * dev) * mu_g
  # But sum(pi * dev) = 0 always, so:
  # d VarB / d pi_g = dev_g^2 - 2 * mu_g * 0 = dev_g^2
  # Wait, that's not quite right. Let me redo this properly.
  # VarB = sum(pi * dev^2) where dev = mu - mu_bar, mu_bar = sum(pi * mu)
  # d VarB / d pi_g = dev_g^2 + sum(pi * 2*dev * d(dev)/d(pi_g))
  # d(dev_h) / d(pi_g) = -d(mu_bar)/d(pi_g) = -mu_g
  # d VarB / d pi_g = dev_g^2 + 2 * sum(pi * dev * (-mu_g))
  #                 = dev_g^2 - 2 * mu_g * sum(pi * dev)
  #                 = dev_g^2 - 0 = dev_g^2
  grad[1:G] <- dev^2

  # d VarB / d mu_g = 2 * pi_g * (mu_g - mu_bar) * (1 - pi_g)
  #                 - 2 * sum_{h!=g} pi_h * (mu_h - mu_bar) * pi_g
  # = 2 * pi_g * [dev_g * (1 - pi_g) - pi_g * sum_{h!=g} pi_h * dev_h]
  # = 2 * pi_g * [dev_g - pi_g * dev_g - pi_g * (sum(pi*dev) - pi_g*dev_g)]
  # = 2 * pi_g * [dev_g - pi_g * sum(pi*dev)]
  # = 2 * pi_g * dev_g   (since sum(pi*dev) = 0)
  grad[(G+1):(2*G)] <- 2 * pi * dev

  # d VarB / d sigma2_g = 0
  # (already zero)

  grad
}
