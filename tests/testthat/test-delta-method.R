# ============================================================================ #
# Tests for delta method standard errors
# ============================================================================ #

# ---------------------------------------------------------------------------- #
# Gradient correctness: analytical vs numerical
# ---------------------------------------------------------------------------- #

test_that("analytical grad_delta_B matches numerical gradient (Var)", {
  pi <- c(0.5, 0.5)
  mu0 <- c(500, 1000)
  beta <- c(60, 100)
  J <- 2

  # Analytical gradient
  grad_analytic <- ineqx:::.grad_delta_B(pi, mu0, beta, "Var")

  # Numerical gradient
  grad_numeric <- ineqx:::.numerical_grad(function(theta) {
    b <- theta[1:J]
    m <- theta[(J + 1):(2 * J)]
    ineqx:::compute_delta_B(pi, m, b, "Var")
  }, c(beta, mu0))

  expect_equal(grad_analytic, grad_numeric, tolerance = 1e-5)
})

test_that("analytical grad_delta_W matches numerical gradient (Var)", {
  pi <- c(0.5, 0.5)
  mu0 <- c(500, 1000)
  sigma0 <- c(200, 400)
  beta <- c(60, 100)
  lambda <- c(-0.1, -0.2)
  J <- 2

  # Analytical gradient
  grad_analytic <- ineqx:::.grad_delta_W(pi, mu0, sigma0, beta, lambda, "Var")

  # Numerical gradient
  grad_numeric <- ineqx:::.numerical_grad(function(theta) {
    l <- theta[1:J]
    ls <- theta[(J + 1):(2 * J)]
    s0 <- exp(ls)
    s1 <- s0 * exp(l)
    # delta_W = sum(pi * s1^2) - sum(pi * s0^2)
    sum(pi * s1^2) - sum(pi * s0^2)
  }, c(lambda, log(sigma0)))

  expect_equal(grad_analytic, grad_numeric, tolerance = 1e-5)
})

test_that("grad_delta_B works with 3 groups", {
  pi <- c(0.3, 0.5, 0.2)
  mu0 <- c(100, 200, 300)
  beta <- c(10, 20, 15)

  grad <- ineqx:::.grad_delta_B(pi, mu0, beta, "Var")
  expect_length(grad, 6)  # 2 * J = 6

  # Verify against numerical
  grad_num <- ineqx:::.numerical_grad(function(theta) {
    b <- theta[1:3]
    m <- theta[4:6]
    ineqx:::compute_delta_B(pi, m, b, "Var")
  }, c(beta, mu0))

  expect_equal(grad, grad_num, tolerance = 1e-5)
})

# ---------------------------------------------------------------------------- #
# SE computation: known-answer test
# ---------------------------------------------------------------------------- #

test_that("cross-sectional SEs compute correctly with known vcov", {
  params <- ineqx_params(
    data = data.frame(
      group = c("workers", "managers"),
      pi = c(0.5, 0.5),
      mu0 = c(500, 1000),
      sigma0 = c(200, 400),
      beta = c(60, 100),
      lambda = c(-0.1, -0.2)
    ),
    # Construct a simple diagonal vcov (4J = 8)
    vcov = diag(c(
      10, 10,     # se(beta)^2
      25, 25,     # se(mu0)^2
      0.01, 0.01, # se(lambda)^2
      0.04, 0.04  # se(log_sigma0)^2
    ))
  )

  result <- ineqx:::causal_decompose_cross(params)

  # SEs should be present
  expect_true(!is.null(result$se))
  expect_true(result$se$se_delta_B > 0)
  expect_true(result$se$se_delta_W > 0)
  expect_true(result$se$se_delta_total > 0)

  # Verify manually: grad_B = 2*pi*(mu1 - mu1_bar, beta - beta_bar)
  # mu1 = (560, 1100), mu1_bar = 830
  # grad_beta = 2*0.5*(560-830, 1100-830) = (-270, 270)
  # beta_bar = 80
  # grad_mu0 = 2*0.5*(60-80, 100-80) = (-20, 20)
  # grad_B = (-270, 270, -20, 20)
  # Var(delta_B) = 270^2*10 + 270^2*10 + 20^2*25 + 20^2*25
  #              = 2*270^2*10 + 2*20^2*25 = 1458000 + 20000 = 1478000
  expected_se_B <- sqrt(1478000)
  expect_equal(result$se$se_delta_B, expected_se_B, tolerance = 1e-4)
})

test_that("zero vcov gives zero SEs", {
  params <- ineqx_params(
    data = data.frame(
      group = c("A", "B"),
      pi = c(0.5, 0.5),
      mu0 = c(500, 1000),
      sigma0 = c(200, 400),
      beta = c(60, 100),
      lambda = c(-0.1, -0.2)
    ),
    vcov = matrix(0, 8, 8)
  )

  result <- ineqx:::causal_decompose_cross(params)
  expect_equal(result$se$se_delta_B, 0)
  expect_equal(result$se$se_delta_W, 0)
  expect_equal(result$se$se_delta_total, 0)
})

test_that("block independence: off-diagonal in W block doesn't affect B SE", {
  base_vcov <- diag(c(10, 10, 25, 25, 0.01, 0.01, 0.04, 0.04))

  # Add off-diagonal terms only in the W block (lambda, log_sigma0)
  vcov_modified <- base_vcov
  vcov_modified[5, 6] <- 0.005
  vcov_modified[6, 5] <- 0.005
  vcov_modified[7, 8] <- 0.02
  vcov_modified[8, 7] <- 0.02

  params1 <- ineqx_params(
    data = data.frame(
      group = c("A", "B"),
      pi = c(0.5, 0.5),
      mu0 = c(500, 1000),
      sigma0 = c(200, 400),
      beta = c(60, 100),
      lambda = c(-0.1, -0.2)
    ),
    vcov = base_vcov
  )

  params2 <- ineqx_params(
    data = data.frame(
      group = c("A", "B"),
      pi = c(0.5, 0.5),
      mu0 = c(500, 1000),
      sigma0 = c(200, 400),
      beta = c(60, 100),
      lambda = c(-0.1, -0.2)
    ),
    vcov = vcov_modified
  )

  se1 <- ineqx:::delta_method_se(params1, type = "cross")
  se2 <- ineqx:::delta_method_se(params2, type = "cross")

  # delta_B SE should be unchanged
  expect_equal(se1$se_delta_B, se2$se_delta_B)

  # delta_W SE should change
  expect_false(se1$se_delta_W == se2$se_delta_W)
})

# ---------------------------------------------------------------------------- #
# Longitudinal SEs
# ---------------------------------------------------------------------------- #

test_that("longitudinal delta method SEs are computed", {
  params <- ineqx_params(
    data = data.frame(
      group = rep(c("A", "B"), 2),
      time = rep(c(2000, 2010), each = 2),
      pi = c(0.4, 0.6, 0.45, 0.55),
      mu0 = c(500, 1000, 550, 1100),
      sigma0 = c(200, 400, 210, 380),
      beta = c(60, 100, 80, 90),
      lambda = c(-0.1, -0.2, -0.15, -0.15)
    ),
    ref = 2000,
    vcov = list(
      "2000" = diag(rep(1, 8)),
      "2010" = diag(rep(1, 8))
    )
  )

  result <- ineqx:::causal_decompose_longit(params)

  # SEs should be present
  expect_true(!is.null(result$se))
  expect_true("2010" %in% names(result$se))

  se_2010 <- result$se[["2010"]]
  # All SEs should be non-negative
  for (nm in names(se_2010)) {
    expect_true(se_2010[[nm]] >= 0, info = paste("SE for", nm))
  }

  # Total SE should be present
  expect_true("se_Delta_total" %in% names(se_2010))
})

# ---------------------------------------------------------------------------- #
# Zero treatment effects
# ---------------------------------------------------------------------------- #

test_that("zero treatment effects give zero delta_B but nonzero grad_beta", {
  pi <- c(0.5, 0.5)
  mu0 <- c(500, 1000)
  beta <- c(0, 0)

  # delta_B itself is zero
  expect_equal(ineqx:::compute_delta_B(pi, mu0, beta, "Var"), 0)

  grad <- ineqx:::.grad_delta_B(pi, mu0, beta, "Var")

  # grad_mu0 should be zero (beta - beta_bar = 0 - 0 = 0)
  expect_equal(grad[3], 0)
  expect_equal(grad[4], 0)

  # grad_beta should NOT be zero: 2*pi*(mu0 - mu0_bar) != 0
  # because mu0 values differ
  expect_false(grad[1] == 0)
  expect_false(grad[2] == 0)
  # But they should sum to 0 (opposite signs)
  expect_equal(grad[1] + grad[2], 0)
})

test_that("zero lambda gives zero gradient for delta_W w.r.t. log_sigma0", {
  pi <- c(0.5, 0.5)
  mu0 <- c(500, 1000)
  sigma0 <- c(200, 400)
  beta <- c(60, 100)
  lambda <- c(0, 0)

  grad <- ineqx:::.grad_delta_W(pi, mu0, sigma0, beta, lambda, "Var")

  # grad_lambda should be 2*pi*sigma0^2*exp(0) = 2*pi*sigma0^2
  expect_equal(grad[1], 2 * 0.5 * 200^2)
  expect_equal(grad[2], 2 * 0.5 * 400^2)

  # grad_log_sigma0 should be 0 (f = exp(0)-1 = 0)
  expect_equal(grad[3], 0)
  expect_equal(grad[4], 0)
})

# ---------------------------------------------------------------------------- #
# p-value and significance helpers
# ---------------------------------------------------------------------------- #

test_that("significance stars are correct", {
  expect_equal(ineqx:::.signif_stars(0.0001), "***")
  expect_equal(ineqx:::.signif_stars(0.005), "**")
  expect_equal(ineqx:::.signif_stars(0.03), "*")
  expect_equal(ineqx:::.signif_stars(0.08), ".")
  expect_equal(ineqx:::.signif_stars(0.5), "")
  expect_equal(ineqx:::.signif_stars(NA), "")
})

test_that("p-value from SE is correct", {
  # z = 3 -> p ≈ 0.0027
  pval <- ineqx:::.pval_from_se(3, 1)
  expect_equal(pval, 2 * pnorm(-3), tolerance = 1e-10)

  # SE = 0 -> NA
  expect_true(is.na(ineqx:::.pval_from_se(1, 0)))
  expect_true(is.na(ineqx:::.pval_from_se(1, NA)))
})
