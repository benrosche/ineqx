# ============================================================================ #
# Tests for the split-step workflow on a stepwise-DiD dataset.
#
# fit_ineqx_model() + ineqx_params(post = ...) + ineqx(params = ...) should
# produce bit-identical point estimates to the integrated ineqx(post = ...,
# data = ..., ...) call on the same data.
# ============================================================================ #

skip_if_not_installed("gamlss")

.simulate_split_did <- function(n_per_cell = 800, seed = 1) {
  set.seed(seed)
  alpha   <- c(A = 500, B = 800)
  beta_D  <- c(A = 30,  B = 60)
  beta_P  <- c(A = 20,  B = 40)
  beta_DP <- c(A = -50, B = -120)
  log_sig <- c(A = log(150), B = log(250))
  lam_D   <- c(A = 0.05, B = 0.05)
  lam_P   <- c(A = 0.02, B = 0.03)
  lam_DP  <- c(A = -0.10, B = -0.20)

  rows <- list()
  for (g in c("A", "B")) {
    for (D in c(0, 1)) {
      for (P in c(0, 1)) {
        mu_cell  <- alpha[g] + beta_D[g] * D + beta_P[g] * P + beta_DP[g] * D * P
        lsd_cell <- log_sig[g] + lam_D[g] * D + lam_P[g] * P + lam_DP[g] * D * P
        y <- rnorm(n_per_cell, mean = mu_cell, sd = exp(lsd_cell))
        rows[[length(rows) + 1]] <- data.frame(
          group = g, treat = D, post = P, y = y,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, rows)
}

test_that("split-step DiD matches integrated DiD bit-for-bit", {
  d <- .simulate_split_did(n_per_cell = 800, seed = 11)

  # Split-step: fit GAMLSS once, derive params, decompose.
  m <- fit_ineqx_model(
    formula_mu    = y ~ group * treat * post,
    formula_sigma =   ~ group * treat * post,
    data = d
  )
  params <- ineqx_params(
    model = m, data = d,
    treat = "treat", group = "group", post = "post",
    ystat = "Var", vcov = FALSE, verbose = FALSE
  )
  fit_split <- ineqx(params = params, ystat = "Var", se = "none")

  # Integrated: one call, same args.
  fit_integ <- ineqx(
    y = "y", ystat = "Var",
    treat = "treat", group = "group", post = "post",
    formula_mu    = ~ group * treat * post,
    formula_sigma = ~ group * treat * post,
    data = d, se = "none"
  )

  # Point estimates must match (GAMLSS is deterministic given the same data
  # and formula; both code paths fit the same model).
  expect_equal(fit_split$tau_B,     fit_integ$tau_B)
  expect_equal(fit_split$tau_W,     fit_integ$tau_W)
  expect_equal(fit_split$tau_total, fit_integ$tau_total)
})

test_that("split-step DiD preserves is_did flag and post column on params", {
  d <- .simulate_split_did(n_per_cell = 400, seed = 13)
  m <- fit_ineqx_model(
    formula_mu    = y ~ group * treat * post,
    formula_sigma =   ~ group * treat * post,
    data = d
  )
  params <- ineqx_params(
    model = m, data = d,
    treat = "treat", group = "group", post = "post",
    ystat = "Var", vcov = FALSE, verbose = FALSE
  )
  expect_true(isTRUE(params$is_did))
  expect_identical(params$post, "post")
  # DiD-only columns present in the params data
  expect_true(all(c("mu0_pre", "mu1_pre", "sigma0_pre", "sigma1_pre") %in%
                  names(params$data)))
})
