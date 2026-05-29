# ============================================================================ #
# Tests for V_L (variance of log earnings) via the split-step workflow.
#
# Covers:
#   (a) fit_ineqx_model(transform = "log") tags the model and survives the
#       round-trip through ineqx_params() onto the params object.
#   (b) Split-step V_L matches integrated V_L bit-for-bit on the same data.
#   (c) ineqx(params, ystat = "VL") errors when the params came from a
#       transform = "identity" (i.e. raw-y) fit.
# ============================================================================ #

skip_if_not_installed("gamlss")

.simulate_pos_cs <- function(n_per_cell = 800, seed = 1) {
  set.seed(seed)
  alpha   <- c(A = log(500), B = log(800))
  beta_D  <- c(A = -0.10,    B = -0.20)
  log_sig <- c(A = log(0.4), B = log(0.5))
  lam_D   <- c(A = 0.05,     B = -0.05)

  rows <- list()
  for (g in c("A", "B")) {
    for (D in c(0, 1)) {
      mu_cell  <- alpha[g] + beta_D[g] * D
      lsd_cell <- log_sig[g] + lam_D[g] * D
      # Lognormal-ish positive y so log(y) is well-defined for VL.
      y <- exp(rnorm(n_per_cell, mean = mu_cell, sd = exp(lsd_cell)))
      rows[[length(rows) + 1]] <- data.frame(
        group = g, treat = D, y = y,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

test_that("fit_ineqx_model(transform = 'log') tags the model and params", {
  d <- .simulate_pos_cs(n_per_cell = 400, seed = 21)
  m <- fit_ineqx_model(
    formula_mu    = y ~ group * treat,
    formula_sigma =   ~ group * treat,
    data = d, transform = "log"
  )
  expect_identical(attr(m, "ineqx_transform"), "log")

  params <- ineqx_params(
    model = m, data = d,
    treat = "treat", group = "group",
    ystat = "Var", vcov = FALSE, verbose = FALSE
  )
  expect_identical(params$transform, "log")
})

test_that("split-step VL matches integrated VL bit-for-bit", {
  d <- .simulate_pos_cs(n_per_cell = 800, seed = 23)

  # Split-step path: fit on log(y), extract params, ineqx with ystat = "VL".
  m_log <- fit_ineqx_model(
    formula_mu    = y ~ group * treat,
    formula_sigma =   ~ group * treat,
    data = d, transform = "log"
  )
  params_log <- ineqx_params(
    model = m_log, data = d,
    treat = "treat", group = "group",
    ystat = "Var", vcov = FALSE, verbose = FALSE
  )
  fit_split <- suppressWarnings(
    ineqx(params = params_log, ystat = "VL", se = "none")
  )

  # Integrated path: ineqx(y, data, ystat = "VL") does its own log transform.
  fit_integ <- suppressWarnings(
    ineqx(
      y = "y", ystat = "VL",
      treat = "treat", group = "group",
      formula_mu    = ~ group * treat,
      formula_sigma = ~ group * treat,
      data = d, se = "none"
    )
  )

  expect_equal(fit_split$tau_B,     fit_integ$tau_B)
  expect_equal(fit_split$tau_W,     fit_integ$tau_W)
  expect_equal(fit_split$tau_total, fit_integ$tau_total)
})

test_that("ineqx(params, ystat = 'VL') errors when params come from a raw-y fit", {
  d <- .simulate_pos_cs(n_per_cell = 400, seed = 27)
  m_raw <- fit_ineqx_model(
    formula_mu    = y ~ group * treat,
    formula_sigma =   ~ group * treat,
    data = d  # transform = "identity" (default)
  )
  params_raw <- ineqx_params(
    model = m_raw, data = d,
    treat = "treat", group = "group",
    ystat = "Var", vcov = FALSE, verbose = FALSE
  )
  expect_identical(params_raw$transform, "identity")
  expect_error(
    ineqx(params = params_raw, ystat = "VL", se = "none"),
    "log\\(y\\) fit"
  )
})
