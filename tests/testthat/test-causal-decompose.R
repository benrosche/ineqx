# ============================================================================ #
# Tests for causal decomposition
# ============================================================================ #

# Paper example: Section 3.2.4 (workers and managers)
test_that("cross-sectional between-group matches paper example", {
  params <- ineqx_params(
    data = data.frame(
      group = c("workers", "managers"),
      pi = c(0.5, 0.5),
      mu0 = c(500, 1000),
      sigma0 = c(200, 400),
      beta = c(60, 100),
      lambda = c(-0.1, -0.2)
    )
  )

  result <- ineqx(params = params, se = "none")

  # Paper: Var_pi(beta) = 0.5*(60^2 + 100^2) - 80^2 = 400
  expect_equal(result$components$Var_pi_beta, 400)

  # Paper: Cov_pi(mu0, beta) = 0.5*(500*60 + 1000*100) - 750*80 = 5000
  expect_equal(result$components$Cov_pi_mu_beta, 5000)

  # Paper: delta_B = 400 + 2*5000 = 10400
  expect_equal(result$delta_B, 10400)
})

test_that("cross-sectional within-group matches paper example (approximate)", {
  params <- ineqx_params(
    data = data.frame(
      group = c("workers", "managers"),
      pi = c(0.5, 0.5),
      mu0 = c(500, 1000),
      sigma0 = c(200, 400),
      beta = c(60, 100),
      lambda = c(-0.1, -0.2)
    )
  )

  result <- ineqx(params = params, se = "none")

  # Exact: sum(pi * sigma0^2 * (exp(2*lambda) - 1))
  f_w <- exp(2 * (-0.1)) - 1
  f_m <- exp(2 * (-0.2)) - 1
  expected_W <- 0.5 * 200^2 * f_w + 0.5 * 400^2 * f_m
  expect_equal(result$delta_W, expected_W)

  # Should be negative (treatment compresses dispersion)
  expect_lt(result$delta_W, 0)

  # Total = B + W
  expect_equal(result$delta_total, result$delta_B + result$delta_W)
})

test_that("reversed treatment effects flip delta_B sign", {
  params <- ineqx_params(
    data = data.frame(
      group = c("workers", "managers"),
      pi = c(0.5, 0.5),
      mu0 = c(500, 1000),
      sigma0 = c(200, 400),
      beta = c(100, 60),  # reversed: workers benefit more
      lambda = c(-0.1, -0.2)
    )
  )

  result <- ineqx(params = params, se = "none")

  # Paper: Cov becomes -5000, delta_B = 400 + 2*(-5000) = -9600
  expect_equal(result$delta_B, -9600)
})

test_that("longitudinal decomposition sums to total change", {
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
    ref = 2000
  )

  result <- ineqx(params = params, se = "none")
  r <- result$results[["2010"]]

  # 6 components sum to total
  six_sum <- r$Delta_beta + r$Delta_lambda +
             r$Delta_pi_B + r$Delta_pi_W +
             r$Delta_mu + r$Delta_sigma
  expect_equal(six_sum, r$Delta_total, tolerance = 1e-10)

  # 3 combined components sum to total
  three_sum <- r$Delta_behavioral + r$Delta_compositional + r$Delta_pretreatment
  expect_equal(three_sum, r$Delta_total, tolerance = 1e-10)

  # B + W = total
  expect_equal(r$Delta_B + r$Delta_W, r$Delta_total, tolerance = 1e-10)
})

test_that("all 6 orderings give same total", {
  params <- ineqx_params(
    data = data.frame(
      group = rep(c("A", "B", "C"), 2),
      time = rep(c(1, 2), each = 3),
      pi = c(0.3, 0.5, 0.2, 0.35, 0.4, 0.25),
      mu0 = c(100, 200, 300, 120, 210, 280),
      sigma0 = c(30, 50, 40, 35, 45, 42),
      beta = c(10, 20, 15, 12, 18, 16),
      lambda = c(-0.05, 0.1, -0.02, -0.08, 0.05, -0.03)
    ),
    ref = 1
  )

  perms <- list(
    c("behavioral", "compositional", "pretreatment"),
    c("behavioral", "pretreatment", "compositional"),
    c("compositional", "behavioral", "pretreatment"),
    c("compositional", "pretreatment", "behavioral"),
    c("pretreatment", "behavioral", "compositional"),
    c("pretreatment", "compositional", "behavioral")
  )

  totals <- vapply(perms, function(ord) {
    r <- ineqx(params = params, order = ord, se = "none")
    r$results[["2"]]$Delta_total
  }, numeric(1))

  # All should be equal
  expect_equal(max(totals) - min(totals), 0, tolerance = 1e-10)
})

test_that("Shapley values average correctly", {
  params <- ineqx_params(
    data = data.frame(
      group = rep(c("A", "B"), 2),
      time = rep(c(1, 2), each = 2),
      pi = c(0.4, 0.6, 0.5, 0.5),
      mu0 = c(100, 200, 110, 220),
      sigma0 = c(30, 50, 32, 48),
      beta = c(10, 20, 15, 18),
      lambda = c(-0.05, 0.1, -0.08, 0.08)
    ),
    ref = 1
  )

  shap <- ineqx(params = params, order = "shapley", se = "none")

  # Shapley total should equal the total from any single ordering
  single <- ineqx(params = params, se = "none")
  expect_equal(shap$shapley$Delta_total[1],
               single$results[["2"]]$Delta_total,
               tolerance = 1e-10)

  # Shapley components should sum to total
  s <- shap$shapley[1, ]
  comp_sum <- s$Delta_behavioral + s$Delta_compositional + s$Delta_pretreatment
  expect_equal(comp_sum, s$Delta_total, tolerance = 1e-10)
})

test_that("zero treatment effects give zero decomposition", {
  params <- ineqx_params(
    data = data.frame(
      group = c("A", "B"),
      pi = c(0.5, 0.5),
      mu0 = c(500, 1000),
      sigma0 = c(200, 400),
      beta = c(0, 0),
      lambda = c(0, 0)
    )
  )

  result <- ineqx(params = params, se = "none")

  expect_equal(result$delta_B, 0)
  expect_equal(result$delta_W, 0)
  expect_equal(result$delta_total, 0)
})
