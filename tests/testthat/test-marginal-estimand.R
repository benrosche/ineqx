# ============================================================================ #
# Tests for the marginal vs residual estimand (paper eqs 10/37-38, App. B.7)
# ============================================================================ #

# Helper: a small dataset with a within-group control z that shifts the
# conditional mean (so V_Z[m] > 0 and marginal != residual).
.make_marginal_data <- function(seed = 1, n_per = 200) {
  set.seed(seed)
  g <- rep(c("A", "B", "C"), each = n_per)
  x <- rbinom(length(g), 1, 0.5)
  z <- rnorm(length(g))
  base <- c(A = 0, B = 2, C = 4)[g]
  beta <- c(A = 0.5, B = 1.0, C = 1.5)[g]
  inc <- base + beta * x + 0.8 * z + rnorm(length(g), sd = 1)
  data.frame(inc = inc, group = g, x = x, z = z, stringsAsFactors = FALSE)
}

test_that("estimand defaults to marginal and is recorded on the object", {
  d <- .make_marginal_data()
  pm <- suppressMessages(ineqx_params(
    model = fit_ineqx_model(inc ~ x * factor(group) + z, ~ x * factor(group), data = d),
    data = d, treat = "x", group = "group", vcov = FALSE, verbose = FALSE))
  expect_equal(pm$estimand, "marginal")

  pr <- suppressMessages(ineqx_params(
    model = fit_ineqx_model(inc ~ x * factor(group) + z, ~ x * factor(group), data = d),
    data = d, treat = "x", group = "group", estimand = "residual",
    vcov = FALSE, verbose = FALSE))
  expect_equal(pr$estimand, "residual")
})

test_that("marginal group variance exceeds residual when a control varies within group", {
  d <- .make_marginal_data()
  m <- fit_ineqx_model(inc ~ x * factor(group) + z, ~ x * factor(group), data = d)
  pm <- suppressMessages(ineqx_params(model = m, data = d, treat = "x",
                                      group = "group", estimand = "marginal",
                                      vcov = FALSE, verbose = FALSE))
  pr <- suppressMessages(ineqx_params(model = m, data = d, treat = "x",
                                      group = "group", estimand = "residual",
                                      vcov = FALSE, verbose = FALSE))
  # V_Z[m] term makes every marginal baseline SD strictly larger.
  expect_true(all(pm$data$sigma0 > pr$data$sigma0 + 1e-6))
})

test_that("marginal and residual coincide when there is no within-group covariate variation", {
  d <- .make_marginal_data()
  # No control in either equation -> predictions constant within (group, x) ->
  # V_Z[m] = 0 -> marginal moments equal residual moments.
  m <- fit_ineqx_model(inc ~ x * factor(group), ~ x * factor(group), data = d)
  pm <- suppressMessages(ineqx_params(model = m, data = d, treat = "x",
                                      group = "group", estimand = "marginal",
                                      vcov = FALSE, verbose = FALSE))
  pr <- suppressMessages(ineqx_params(model = m, data = d, treat = "x",
                                      group = "group", estimand = "residual",
                                      vcov = FALSE, verbose = FALSE))
  expect_equal(pm$data$sigma0, pr$data$sigma0, tolerance = 1e-8)
  expect_equal(pm$data$lambda, pr$data$lambda, tolerance = 1e-8)
})

test_that("marginal moments satisfy the identity sigma1 == sigma0 * exp(lambda)", {
  d <- .make_marginal_data()
  m <- fit_ineqx_model(inc ~ x * factor(group) + z, ~ x * factor(group), data = d)
  pm <- suppressMessages(ineqx_params(model = m, data = d, treat = "x",
                                      group = "group", estimand = "marginal",
                                      vcov = FALSE, verbose = FALSE))
  expect_equal(pm$data$sigma1, pm$data$sigma0 * exp(pm$data$lambda),
               tolerance = 1e-10)
})

test_that("numerical (marginal) Jacobian reproduces the analytical (residual) one with no control", {
  # With no within-group control the marginal moments equal the residual ones
  # and there is no Jensen gap (predictions are constant within group), so the
  # numerical marginal Jacobian must reproduce the fast analytical residual SEs.
  d <- .make_marginal_data()
  fm <- suppressMessages(ineqx(y = "inc", treat = "x", group = "group",
                               estimand = "marginal",
                               formula_mu = ~ x * factor(group),
                               formula_sigma = ~ x * factor(group),
                               data = d, se = "delta"))
  fr <- suppressMessages(ineqx(y = "inc", treat = "x", group = "group",
                               estimand = "residual",
                               formula_mu = ~ x * factor(group),
                               formula_sigma = ~ x * factor(group),
                               data = d, se = "delta"))
  expect_equal(fm[["se"]]$se_tau_B, fr[["se"]]$se_tau_B, tolerance = 1e-4)
  expect_equal(fm[["se"]]$se_tau_W, fr[["se"]]$se_tau_W, tolerance = 1e-4)
  expect_equal(fm[["se"]]$se_tau_total, fr[["se"]]$se_tau_total, tolerance = 1e-4)
})

test_that("marginal vcov couples mean- and scale-equation uncertainty (incremental FD)", {
  # The marginal S0/Lambda depend on the mean-equation coefficients (via the
  # V_Z[m] term), so the theta-space vcov gets a non-zero between x within
  # cross-block. The incremental finite-difference Jacobian must preserve this;
  # with no within-group control the coupling vanishes and the cross-block is 0.
  # Theta layout per time: [beta_1..J, mu0_1..J, lambda_1..J, log_sigma0_1..J],
  # so for J groups the B block is 1:(2J) and the W block is (2J+1):(4J).
  d <- .make_marginal_data()
  J <- 3L
  bw <- function(fit) {
    V <- ineqx:::.get_vcov_for_time(fit$params, NULL)
    max(abs(V[1:(2 * J), (2 * J + 1):(4 * J)]))
  }

  fit_ctrl <- suppressMessages(ineqx(
    y = "inc", treat = "x", group = "group", estimand = "marginal",
    formula_mu = ~ x * factor(group) + z, formula_sigma = ~ x * factor(group),
    data = d, se = "delta"))
  fit_noctrl <- suppressMessages(ineqx(
    y = "inc", treat = "x", group = "group", estimand = "marginal",
    formula_mu = ~ x * factor(group), formula_sigma = ~ x * factor(group),
    data = d, se = "delta"))

  expect_gt(bw(fit_ctrl), 1e-6)     # controls -> mean coefs feed S0 -> coupling
  expect_lt(bw(fit_noctrl), 1e-8)   # no controls -> block-diagonal Sigma_theta
})

test_that("marginal delta SEs are finite and positive with a within-group control", {
  d <- .make_marginal_data()
  fm <- suppressMessages(ineqx(y = "inc", treat = "x", group = "group",
                               estimand = "marginal",
                               formula_mu = ~ x * factor(group) + z,
                               formula_sigma = ~ x * factor(group),
                               data = d, se = "delta"))
  expect_true(is.finite(fm[["se"]]$se_tau_W) && fm[["se"]]$se_tau_W > 0)
  expect_true(is.finite(fm[["se"]]$se_tau_B) && fm[["se"]]$se_tau_B > 0)
  # The marginal numerical vcov couples mean- and variance-equation
  # coefficients, so it is a full (non-block-diagonal) 4J x 4J matrix.
  V <- ineqx:::.get_vcov_for_time(fm$params, NULL)
  expect_equal(dim(V), c(12L, 12L))
})
