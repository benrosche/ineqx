# ============================================================================ #
# Tests for the two-stage bootstrap (bootstrap_params + decompose_boot_params).
#
# The two-stage flow caches the resampled params per replicate and only
# re-decomposes on demand. Given the same seed, the same data, and the same
# (ystat, ref), it must produce numerically identical SEs to the single-stage
# bootstrap_se() that re-fits GAMLSS per replicate.
# ============================================================================ #

skip_if_not_installed("gamlss")

.simulate_cs <- function(n_per_cell = 400, seed = 1) {
  set.seed(seed)
  alpha   <- c(A = 500, B = 800)
  beta_D  <- c(A = -50, B = -120)
  log_sig <- c(A = log(150), B = log(250))
  lam_D   <- c(A = 0.05, B = 0.05)

  rows <- list()
  for (g in c("A", "B")) {
    for (D in c(0, 1)) {
      mu_cell  <- alpha[g] + beta_D[g] * D
      lsd_cell <- log_sig[g] + lam_D[g] * D
      y <- rnorm(n_per_cell, mean = mu_cell, sd = exp(lsd_cell))
      rows[[length(rows) + 1]] <- data.frame(
        group = g, treat = D, y = y,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

test_that("two-stage bootstrap matches single-stage bootstrap_se", {
  d <- .simulate_cs(n_per_cell = 300, seed = 31)

  bp <- bootstrap_params(
    data = d,
    formula_mu    = y ~ group * treat,
    formula_sigma =   ~ group * treat,
    treat = "treat", group = "group",
    ystat = "Var", B = 5, seed = 42, verbose = FALSE
  )
  expect_s3_class(bp, "ineqx_boot_params")
  expect_equal(bp$B, 5L)
  expect_true(bp$B_successful >= 2L)
  expect_length(bp$params_list, bp$B_successful)

  boot_two_stage <- decompose_boot_params(bp, ref = NULL, ystat = "Var")

  boot_one_stage <- ineqx:::bootstrap_se(
    data = d,
    formula_mu    = y ~ group * treat,
    formula_sigma =   ~ group * treat,
    treat = "treat", group = "group",
    ystat = "Var", B = 5, seed = 42, verbose = FALSE
  )

  # Same seed -> same resample indices -> same GAMLSS fits -> same params ->
  # same decomposition -> same SEs.
  expect_equal(boot_two_stage$se, boot_one_stage$se)
  expect_equal(unname(boot_two_stage$replicates),
               unname(boot_one_stage$replicates))
})

test_that("decompose_boot_params reuses cached params for different ystat", {
  d <- .simulate_cs(n_per_cell = 300, seed = 33)

  bp <- bootstrap_params(
    data = d,
    formula_mu    = y ~ group * treat,
    formula_sigma =   ~ group * treat,
    treat = "treat", group = "group",
    ystat = "Var", B = 5, seed = 7, verbose = FALSE
  )

  # Decomposing at CV2 should not error and should yield SEs for tau_B/tau_W.
  se_cv2 <- decompose_boot_params(bp, ystat = "CV2")
  expect_s3_class(se_cv2, "ineqx_boot")
  expect_true(is.numeric(se_cv2$se$se_tau_B))
  expect_true(is.numeric(se_cv2$se$se_tau_W))
})

test_that("bootstrap_params errors gracefully on bad inputs", {
  d <- .simulate_cs(n_per_cell = 100, seed = 35)
  expect_error(
    bootstrap_params(
      data = d,
      formula_mu    = y ~ group * treat,
      formula_sigma =   ~ group * treat,
      treat = "treat", group = "group",
      ystat = "Var", B = 1
    ),
    "B.*at least 2"
  )

  expect_error(
    decompose_boot_params("not a boot_params object"),
    "ineqx_boot_params"
  )
})
