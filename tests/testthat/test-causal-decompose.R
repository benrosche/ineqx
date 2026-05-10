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

  result <- ineqx("y", group = "g", data = data.frame(y = 1, g = "A"),
                  params = params, se = "none")

  # Paper: Var_pi(beta) = 0.5*(60^2 + 100^2) - 80^2 = 400
  expect_equal(result$components$Var_pi_beta, 400)

  # Paper: Cov_pi(mu0, beta) = 0.5*(500*60 + 1000*100) - 750*80 = 5000
  expect_equal(result$components$Cov_pi_mu_beta, 5000)

  # Paper: tau_B = 400 + 2*5000 = 10400
  expect_equal(result$tau_B, 10400)
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

  result <- ineqx("y", group = "g", data = data.frame(y = 1, g = "A"),
                  params = params, se = "none")

  # Exact: sum(pi * sigma0^2 * (exp(2*lambda) - 1))
  f_w <- exp(2 * (-0.1)) - 1
  f_m <- exp(2 * (-0.2)) - 1
  expected_W <- 0.5 * 200^2 * f_w + 0.5 * 400^2 * f_m
  expect_equal(result$tau_W, expected_W)

  # Should be negative (treatment compresses dispersion)
  expect_lt(result$tau_W, 0)

  # Total = B + W
  expect_equal(result$tau_total, result$tau_B + result$tau_W)
})

test_that("reversed treatment effects flip tau_B sign", {
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

  result <- ineqx("y", group = "g", data = data.frame(y = 1, g = "A"),
                  params = params, se = "none")

  # Paper: Cov becomes -5000, tau_B = 400 + 2*(-5000) = -9600
  expect_equal(result$tau_B, -9600)
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
    )
  )

  result <- ineqx("y", group = "g", data = data.frame(y = 1, g = "A"),
                  params = params, ref = 2000,
                  order = c("behavioral", "compositional", "pretreatment"),
                  se = "none")
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

test_that("V longitudinal: off-diagonal split parts are 0; aggregates unchanged", {
  params <- ineqx_params(
    data = data.frame(
      group = rep(c("A", "B"), 2),
      time = rep(c(2000, 2010), each = 2),
      pi = c(0.4, 0.6, 0.45, 0.55),
      mu0 = c(500, 1000, 550, 1100),
      sigma0 = c(200, 400, 210, 380),
      beta = c(60, 100, 80, 90),
      lambda = c(-0.1, -0.2, -0.15, -0.15)
    )
  )

  result <- ineqx("y", group = "g", data = data.frame(y = 1, g = "A"),
                  params = params, ref = 2000,
                  order = c("behavioral", "compositional", "pretreatment"),
                  se = "none")
  r <- result$results[["2010"]]

  # For V, beta and mu only enter tau_B; lambda and sigma only enter tau_W.
  expect_equal(r$Delta_beta_W,  0, tolerance = 1e-10)
  expect_equal(r$Delta_mu_W,    0, tolerance = 1e-10)
  expect_equal(r$Delta_lambda_B, 0, tolerance = 1e-10)
  expect_equal(r$Delta_sigma_B,  0, tolerance = 1e-10)

  # Aggregate fields equal their non-zero side
  expect_equal(r$Delta_beta,   r$Delta_beta_B,   tolerance = 1e-10)
  expect_equal(r$Delta_lambda, r$Delta_lambda_W, tolerance = 1e-10)
  expect_equal(r$Delta_mu,     r$Delta_mu_B,     tolerance = 1e-10)
  expect_equal(r$Delta_sigma,  r$Delta_sigma_W,  tolerance = 1e-10)

  # Delta_total matches cross-sectional change exactly
  expect_equal(r$Delta_total,
               (r$tau_B_t + r$tau_W_t) - (r$tau_B_t0 + r$tau_W_t0),
               tolerance = 1e-10)
})

test_that("CV2 longitudinal decomposition sums to cross-sectional change", {
  params <- ineqx_params(
    data = data.frame(
      group  = rep(c("A", "B"), 2),
      time   = rep(c(2000, 2010), each = 2),
      pi     = c(0.4, 0.6, 0.45, 0.55),
      mu0    = c(500, 1000, 550, 1100),
      sigma0 = c(200, 400, 210, 380),
      beta   = c(60, 100, 80, 90),
      lambda = c(-0.1, -0.2, -0.15, -0.15)
    ),
    ystat = "CV2"
  )

  result <- ineqx("y", group = "g", data = data.frame(y = 1, g = "A"),
                  params = params, ref = 2000, ystat = "CV2",
                  order = c("behavioral", "compositional", "pretreatment"),
                  se = "none")
  r <- result$results[["2010"]]

  # NEW property for CV2: decomp total equals cross-sectional change
  expect_equal(r$Delta_total,
               (r$tau_B_t + r$tau_W_t) - (r$tau_B_t0 + r$tau_W_t0),
               tolerance = 1e-10)

  # 10 split components sum to total
  ten_sum <- with(r,
    Delta_beta_B + Delta_beta_W + Delta_lambda_B + Delta_lambda_W +
    Delta_pi_B   + Delta_pi_W   + Delta_mu_B     + Delta_mu_W +
    Delta_sigma_B + Delta_sigma_W
  )
  expect_equal(ten_sum, r$Delta_total, tolerance = 1e-10)

  # B-side sums to cross-sectional Delta_B
  b_sum <- with(r,
    Delta_beta_B + Delta_lambda_B + Delta_pi_B + Delta_mu_B + Delta_sigma_B
  )
  expect_equal(b_sum, r$tau_B_t - r$tau_B_t0, tolerance = 1e-10)

  # W-side sums to cross-sectional Delta_W
  w_sum <- with(r,
    Delta_beta_W + Delta_lambda_W + Delta_pi_W + Delta_mu_W + Delta_sigma_W
  )
  expect_equal(w_sum, r$tau_W_t - r$tau_W_t0, tolerance = 1e-10)

  # CV2: at least one off-diagonal split part should be nonzero
  # (shared denominator means parameters cross-contribute)
  off_diag <- abs(r$Delta_beta_W) + abs(r$Delta_mu_W) +
              abs(r$Delta_lambda_B) + abs(r$Delta_sigma_B)
  expect_gt(off_diag, 0)
})

test_that("CV2 longitudinal: shapley preserves cross-sectional total", {
  params <- ineqx_params(
    data = data.frame(
      group  = rep(c("A", "B"), 2),
      time   = rep(c(2000, 2010), each = 2),
      pi     = c(0.4, 0.6, 0.45, 0.55),
      mu0    = c(500, 1000, 550, 1100),
      sigma0 = c(200, 400, 210, 380),
      beta   = c(60, 100, 80, 90),
      lambda = c(-0.1, -0.2, -0.15, -0.15)
    ),
    ystat = "CV2"
  )

  shap <- ineqx("y", group = "g", data = data.frame(y = 1, g = "A"),
                params = params, ref = 2000, ystat = "CV2",
                order = "shapley", se = "none")
  r <- shap$results[["2010"]]

  # Shapley total still matches cross-sectional change
  expect_equal(r$Delta_total,
               (r$tau_B_t + r$tau_W_t) - (r$tau_B_t0 + r$tau_W_t0),
               tolerance = 1e-10)

  # 10 split components sum to total under Shapley averaging
  ten_sum <- with(r,
    Delta_beta_B + Delta_beta_W + Delta_lambda_B + Delta_lambda_W +
    Delta_pi_B   + Delta_pi_W   + Delta_mu_B     + Delta_mu_W +
    Delta_sigma_B + Delta_sigma_W
  )
  expect_equal(ten_sum, r$Delta_total, tolerance = 1e-10)
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
    )
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
    r <- ineqx("y", group = "g", data = data.frame(y = 1, g = "A"),
               params = params, ref = 1, order = ord, se = "none")
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
    )
  )

  shap <- ineqx("y", group = "g", data = data.frame(y = 1, g = "A"),
                params = params, ref = 1, order = "shapley", se = "none")

  # Shapley total should equal the total from any single ordering
  single <- ineqx("y", group = "g", data = data.frame(y = 1, g = "A"),
                  params = params, ref = 1,
                  order = c("behavioral", "compositional", "pretreatment"),
                  se = "none")
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

  result <- ineqx("y", group = "g", data = data.frame(y = 1, g = "A"),
                  params = params, se = "none")

  expect_equal(result$tau_B, 0)
  expect_equal(result$tau_W, 0)
  expect_equal(result$tau_total, 0)
})


# ---------------------------------------------------------------------------- #
# sigma0 = 0 accepted
# ---------------------------------------------------------------------------- #

test_that("sigma0 = 0 is accepted in manual params", {
  params <- ineqx_params(
    data = data.frame(
      group = c("A", "B"),
      pi = c(0.5, 0.5),
      mu0 = c(500, 1000),
      sigma0 = c(0, 0),
      beta = c(60, 100),
      lambda = c(-0.1, -0.2)
    )
  )

  expect_s3_class(params, "ineqx_params")

  # Within-group effect should be zero (sigma0=0 means no baseline dispersion)
  result <- ineqx("y", group = "g", data = data.frame(y = 1, g = "A"),
                  params = params, se = "none")
  expect_equal(result$tau_W, 0)
  # Between-group effect should still work
  expect_true(result$tau_B != 0)
})

test_that("4-component decomposition and sub-components are consistent", {
  params <- ineqx_params(
    data = data.frame(
      group = rep(c("A", "B", "C"), 2),
      time = rep(c(1, 2), each = 3),
      pi = c(0.3, 0.5, 0.2, 0.35, 0.4, 0.25),
      mu0 = c(100, 200, 300, 120, 210, 280),
      sigma0 = c(30, 50, 40, 35, 45, 42),
      beta = c(10, 20, 15, 12, 18, 16),
      lambda = c(-0.05, 0.1, -0.02, -0.08, 0.05, -0.03)
    )
  )

  result <- ineqx("y", group = "g", data = data.frame(y = 1, g = "A"),
                  params = params, ref = 1,
                  order = c("behavioral", "compositional", "pretreatment"),
                  se = "none")
  r <- result$results[["2"]]

  # 4-component: Delta_beta + Delta_lambda + Delta_pi + Delta_pre = total
  four_sum <- r$Delta_beta + r$Delta_lambda + r$Delta_pi + r$Delta_pre
  expect_equal(four_sum, r$Delta_total, tolerance = 1e-10)

  # Delta_pi = Delta_pi_B + Delta_pi_W
  expect_equal(r$Delta_pi, r$Delta_pi_B + r$Delta_pi_W, tolerance = 1e-10)

  # Delta_pre = Delta_mu + Delta_sigma
  expect_equal(r$Delta_pre, r$Delta_mu + r$Delta_sigma, tolerance = 1e-10)

  # Sub-components at each time: Var_pi(beta) + 2*Cov_pi(mu0,beta) = tau_B
  ct <- r$components_t
  tau_B_from_comps <- ct$Var_pi_beta + 2 * ct$Cov_pi_mu_beta
  expect_equal(tau_B_from_comps, r$tau_B_t, tolerance = 1e-10)

  # Sub-components: mean(sigma0^2)*mean(f) + Cov(sigma0^2, f) = tau_W
  tau_W_from_comps <- ct$mean_sigma2_0 * ct$mean_f + ct$Cov_pi_sigma2_f
  expect_equal(tau_W_from_comps, r$tau_W_t, tolerance = 1e-10)

  # Same for reference period
  c0 <- r$components_t0
  tau_B0 <- c0$Var_pi_beta + 2 * c0$Cov_pi_mu_beta
  expect_equal(tau_B0, r$tau_B_t0, tolerance = 1e-10)
  tau_W0 <- c0$mean_sigma2_0 * c0$mean_f + c0$Cov_pi_sigma2_f
  expect_equal(tau_W0, r$tau_W_t0, tolerance = 1e-10)
})

test_that("sigma0 < 0 still errors", {
  expect_error(
    ineqx_params(
      data = data.frame(
        group = c("A", "B"),
        pi = c(0.5, 0.5),
        mu0 = c(500, 1000),
        sigma0 = c(-1, 400),
        beta = c(60, 100),
        lambda = c(-0.1, -0.2)
      )
    ),
    "sigma0.*non-negative"
  )
})


# ---------------------------------------------------------------------------- #
# .blend_causal_params() unit tests
# ---------------------------------------------------------------------------- #

test_that(".blend_causal_params merges single CS user params", {
  user <- ineqx_params(
    data = data.frame(
      group = c("A", "B"),
      pi = c(0.5, 0.5),
      mu0 = c(0, 0),
      sigma0 = c(1, 1),
      beta = c(0, 0),
      lambda = c(0, 0)
    )
  )

  model <- ineqx_params(
    data = data.frame(
      group = rep(c("A", "B"), 2),
      time = rep(c(1, 2), each = 2),
      pi = c(0.4, 0.6, 0.45, 0.55),
      mu0 = c(500, 1000, 550, 1100),
      sigma0 = c(200, 400, 210, 380),
      beta = c(60, 100, 80, 90),
      lambda = c(-0.1, -0.2, -0.15, -0.15)
    )
  )

  blended <- ineqx:::.blend_causal_params(user, model, ref = 0L)

  expect_s3_class(blended, "ineqx_params")
  expect_equal(blended$type, "longitudinal")
  expect_equal(sort(blended$times), c(0, 1, 2))
  expect_null(blended$ref)
  expect_true(blended$blended)

  # time=0 should come from user params
  t0 <- blended$data[blended$data$time == 0, ]
  expect_equal(t0$beta, c(0, 0))
  expect_equal(t0$mu0, c(0, 0))

  # time=1 should come from model
  t1 <- blended$data[blended$data$time == 1, ]
  expect_equal(t1$mu0, c(500, 1000))
})

test_that(".blend_causal_params overrides on overlap", {
  user <- ineqx_params(
    data = data.frame(
      group = rep(c("A", "B"), 2),
      time = rep(c(0, 1), each = 2),
      pi = c(0.5, 0.5, 0.5, 0.5),
      mu0 = c(0, 0, 99, 99),
      sigma0 = c(1, 1, 1, 1),
      beta = c(0, 0, 0, 0),
      lambda = c(0, 0, 0, 0)
    )
  )

  model <- ineqx_params(
    data = data.frame(
      group = rep(c("A", "B"), 2),
      time = rep(c(1, 2), each = 2),
      pi = c(0.4, 0.6, 0.45, 0.55),
      mu0 = c(500, 1000, 550, 1100),
      sigma0 = c(200, 400, 210, 380),
      beta = c(60, 100, 80, 90),
      lambda = c(-0.1, -0.2, -0.15, -0.15)
    )
  )

  blended <- ineqx:::.blend_causal_params(user, model, ref = 0)

  # time=1 should come from USER (override), not model
  t1 <- blended$data[blended$data$time == 1, ]
  expect_equal(t1$mu0, c(99, 99))

  # time=2 should come from model (no overlap)
  t2 <- blended$data[blended$data$time == 2, ]
  expect_equal(t2$mu0, c(550, 1100))
})

test_that(".blend_causal_params validates group mismatch", {
  user <- ineqx_params(
    data = data.frame(
      group = c("X", "Y"),
      pi = c(0.5, 0.5),
      mu0 = c(0, 0), sigma0 = c(1, 1),
      beta = c(0, 0), lambda = c(0, 0)
    )
  )

  model <- ineqx_params(
    data = data.frame(
      group = c("A", "B"),
      pi = c(0.4, 0.6),
      mu0 = c(500, 1000), sigma0 = c(200, 400),
      beta = c(60, 100), lambda = c(-0.1, -0.2)
    )
  )

  expect_error(
    ineqx:::.blend_causal_params(user, model, ref = 0),
    "Groups.*do not match"
  )
})


# ---------------------------------------------------------------------------- #
# Causal blending via ineqx() (requires GAMLSS)
# ---------------------------------------------------------------------------- #

test_that("causal blending: single CS counterfactual + model estimation", {
  skip_if_not_installed("gamlss")
  data(cps_sample)
  cps_sample <- na.omit(cps_sample[, c("earnweek", "mother", "edu")])

  groups <- sort(unique(cps_sample[["mother"]]))

  # Counterfactual: zero treatment effect baseline
  ref_params <- ineqx_params(
    data = data.frame(
      group = groups,
      pi = c(0.5, 0.5),
      mu0 = c(500, 500),
      sigma0 = c(100, 100),
      beta = c(0, 0),
      lambda = c(0, 0)
    )
  )

  result <- ineqx("earnweek", treat = "mother", group = "mother",
                   params = ref_params, data = cps_sample,
                   formula_mu = ~ mother * factor(edu),
                   formula_sigma = ~ mother * factor(edu),
                   se = "none")

  # Should be longitudinal (blended creates 2 periods: 0 + observed)
  expect_s3_class(result, "ineqx_causal_longit")
  expect_equal(result$order, "shapley")
  expect_equal(result$ref, 0L)
})

test_that("causal blending: se='delta' falls back to 'none'", {
  skip_if_not_installed("gamlss")
  data(cps_sample)
  cps_sample <- na.omit(cps_sample[, c("earnweek", "mother", "edu")])

  groups <- sort(unique(cps_sample[["mother"]]))

  ref_params <- ineqx_params(
    data = data.frame(
      group = groups,
      pi = c(0.5, 0.5),
      mu0 = c(500, 500),
      sigma0 = c(100, 100),
      beta = c(0, 0),
      lambda = c(0, 0)
    )
  )

  expect_message(
    result <- ineqx("earnweek", treat = "mother", group = "mother",
                     params = ref_params, data = cps_sample,
                     formula_mu = ~ mother * factor(edu),
                     formula_sigma = ~ mother * factor(edu),
                     se = "delta"),
    "Delta method SEs not available"
  )

  expect_equal(result$se_method, "none")
})

test_that("causal blending: multi-period user params require ref", {
  skip_if_not_installed("gamlss")
  data(cps_sample)
  cps_sample <- na.omit(cps_sample[, c("earnweek", "mother", "edu")])

  groups <- sort(unique(cps_sample[["mother"]]))

  ref_params <- ineqx_params(
    data = data.frame(
      group = rep(groups, 2),
      time = rep(c(0, 5), each = 2),
      pi = c(0.5, 0.5, 0.5, 0.5),
      mu0 = c(500, 500, 600, 600),
      sigma0 = c(100, 100, 150, 150),
      beta = c(0, 0, 10, 10),
      lambda = c(0, 0, -0.05, -0.05)
    )
  )

  expect_error(
    ineqx("earnweek", treat = "mother", group = "mother",
           params = ref_params, data = cps_sample,
           formula_mu = ~ mother * factor(edu),
           formula_sigma = ~ mother * factor(edu),
           se = "none"),
    "ref.*required.*multiple time"
  )
})
