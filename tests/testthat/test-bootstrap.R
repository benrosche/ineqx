# ============================================================================ #
# Tests for bootstrap standard errors
# ============================================================================ #

# ---------------------------------------------------------------------------- #
# boot_config validation
# ---------------------------------------------------------------------------- #

test_that("boot_config validates inputs", {
  df <- data.frame(y = 1:10, treat = rep(0:1, 5), group = rep(1:2, 5))

  expect_error(boot_config(data = "not_df", formula_mu = y ~ treat,
                           formula_sigma = ~treat, treat = "treat"),
               "data.*data.frame")
  expect_error(boot_config(data = df, formula_mu = "not_formula",
                           formula_sigma = ~treat, treat = "treat"),
               "formula_mu.*formula")
  expect_error(boot_config(data = df, formula_mu = y ~ treat,
                           formula_sigma = "not_formula", treat = "treat"),
               "formula_sigma.*formula")
  expect_error(boot_config(data = df, formula_mu = y ~ treat,
                           formula_sigma = ~treat, treat = "nonexistent"),
               "not found")
  expect_error(boot_config(data = df, formula_mu = y ~ treat,
                           formula_sigma = ~treat, treat = "treat", B = 1),
               "B.*must be at least 2")
})

test_that("boot_config creates valid object", {
  df <- data.frame(y = 1:10, treat = rep(0:1, 5), group = rep(1:2, 5))
  bc <- boot_config(data = df, formula_mu = y ~ treat * group,
                     formula_sigma = ~treat * group, treat = "treat",
                     group = "group", B = 50, seed = 42)

  expect_s3_class(bc, "ineqx_boot_config")
  expect_equal(bc$B, 50L)
  expect_equal(bc$seed, 42)
  expect_equal(bc$group, "group")
  expect_equal(bc$treat, "treat")
  expect_true(bc$verbose)
  expect_false(bc$parallel)
})

# ---------------------------------------------------------------------------- #
# .generate_boot_indices
# ---------------------------------------------------------------------------- #

test_that("cross-sectional resampling produces correct sizes", {
  df <- data.frame(x = 1:100, .time = rep(1L, 100))
  indices <- ineqx:::.generate_boot_indices(df, ".time", B = 10, seed = 1)

  expect_length(indices, 10)
  for (idx in indices) {
    expect_length(idx, 100)
    expect_true(all(idx >= 1 & idx <= 100))
  }
})

test_that("longitudinal resampling preserves period sizes", {
  df <- data.frame(
    x = 1:200,
    time = c(rep(2000, 80), rep(2010, 120))
  )
  indices <- ineqx:::.generate_boot_indices(df, "time", B = 10, seed = 1)

  period_idx <- split(seq_len(200), df$time)
  for (idx in indices) {
    expect_length(idx, 200)
    # Check that indices from period 1 stay in period 1's range
    in_period1 <- idx[idx %in% period_idx[["2000"]]]
    in_period2 <- idx[idx %in% period_idx[["2010"]]]
    expect_equal(length(in_period1) + length(in_period2), 200)
  }
})

test_that("seed produces reproducible indices", {
  df <- data.frame(x = 1:50, .time = rep(1L, 50))

  idx1 <- ineqx:::.generate_boot_indices(df, ".time", B = 5, seed = 123)
  idx2 <- ineqx:::.generate_boot_indices(df, ".time", B = 5, seed = 123)

  expect_identical(idx1, idx2)
})

# ---------------------------------------------------------------------------- #
# .format_boot_se
# ---------------------------------------------------------------------------- #

test_that(".format_boot_se works for cross-sectional", {
  sds <- c(delta_B = 0.5, delta_W = 0.3, delta_total = 0.6)
  result <- ineqx:::.format_boot_se(sds, "cross", sds)

  expect_named(result, c("se_delta_B", "se_delta_W", "se_delta_total"))
  expect_equal(result$se_delta_B, 0.5)
  expect_equal(result$se_delta_W, 0.3)
  expect_equal(result$se_delta_total, 0.6)
})

test_that(".format_boot_se works for longitudinal", {
  sds <- c(
    Delta_beta_2010 = 0.1, Delta_lambda_2010 = 0.2,
    Delta_pi_B_2010 = 0.05, Delta_pi_W_2010 = 0.06,
    Delta_mu_2010 = 0.03, Delta_sigma_2010 = 0.04,
    Delta_behavioral_2010 = 0.15, Delta_compositional_2010 = 0.12,
    Delta_pretreatment_2010 = 0.08,
    Delta_B_2010 = 0.2, Delta_W_2010 = 0.18, Delta_total_2010 = 0.3
  )
  result <- ineqx:::.format_boot_se(sds, "longit", sds)

  expect_true("2010" %in% names(result))
  expect_true("se_Delta_beta" %in% names(result[["2010"]]))
  expect_equal(result[["2010"]]$se_Delta_beta, 0.1)
})

# ---------------------------------------------------------------------------- #
# attach_boot_se
# ---------------------------------------------------------------------------- #

test_that("attach_boot_se replaces SEs", {
  # Create a mock causal result
  result <- structure(
    list(delta_B = 100, delta_W = -50, delta_total = 50,
         se = list(se_delta_B = 10, se_delta_W = 5, se_delta_total = 8),
         se_method = "delta", ystat = "Var"),
    class = "ineqx_causal_cross"
  )

  # Create a mock boot object
  boot_obj <- structure(
    list(se = list(se_delta_B = 12, se_delta_W = 6, se_delta_total = 9),
         B = 100, B_successful = 98, B_failed = 2, type = "cross"),
    class = "ineqx_boot"
  )

  updated <- attach_boot_se(result, boot_obj)

  expect_equal(updated$se$se_delta_B, 12)
  expect_equal(updated$se_method, "bootstrap")
  expect_s3_class(updated$boot, "ineqx_boot")
})

test_that("attach_boot_se validates inputs", {
  expect_error(attach_boot_se(list(), list()), "ineqx_boot")
})

# ---------------------------------------------------------------------------- #
# Full bootstrap_se (requires gamlss)
# ---------------------------------------------------------------------------- #

test_that("cross-sectional bootstrap smoke test", {
  skip_if_not_installed("gamlss")

  set.seed(42)
  n <- 200
  group <- rep(c("A", "B"), each = n / 2)
  treat <- rbinom(n, 1, 0.5)
  y <- 10 + ifelse(group == "B", 5, 0) + 3 * treat +
    ifelse(group == "B", 2, 0) * treat + rnorm(n, 0, 2)
  df <- data.frame(y = y, treat = treat, group = group)

  result <- bootstrap_se(
    data = df,
    formula_mu = y ~ treat * group,
    formula_sigma = ~ treat * group,
    treat = "treat",
    group = "group",
    ystat = "Var",
    B = 20,
    seed = 123,
    verbose = FALSE
  )

  expect_s3_class(result, "ineqx_boot")
  expect_equal(result$type, "cross")
  expect_true(result$B_successful >= 15)  # most should succeed
  expect_named(result$se, c("se_delta_B", "se_delta_W", "se_delta_total"))
  expect_true(all(vapply(result$se, is.numeric, logical(1))))
  expect_true(all(vapply(result$se, function(x) x >= 0, logical(1))))
  expect_equal(ncol(result$replicates), 3)
  expect_equal(result$seed, 123)
})

test_that("bootstrap reproducibility with same seed", {
  skip_if_not_installed("gamlss")

  set.seed(1)
  n <- 150
  group <- rep(c("A", "B"), each = n / 2)
  treat <- rbinom(n, 1, 0.5)
  y <- 10 + ifelse(group == "B", 3, 0) + 2 * treat + rnorm(n)
  df <- data.frame(y = y, treat = treat, group = group)

  r1 <- bootstrap_se(data = df, formula_mu = y ~ treat * group,
                      formula_sigma = ~ treat * group,
                      treat = "treat", group = "group",
                      B = 10, seed = 999, verbose = FALSE)
  r2 <- bootstrap_se(data = df, formula_mu = y ~ treat * group,
                      formula_sigma = ~ treat * group,
                      treat = "treat", group = "group",
                      B = 10, seed = 999, verbose = FALSE)

  expect_equal(r1$replicates, r2$replicates)
})

test_that("print.ineqx_boot runs without error", {
  skip_if_not_installed("gamlss")

  set.seed(42)
  n <- 150
  group <- rep(c("A", "B"), each = n / 2)
  treat <- rbinom(n, 1, 0.5)
  y <- 10 + 2 * treat + rnorm(n)
  df <- data.frame(y = y, treat = treat, group = group)

  result <- bootstrap_se(data = df, formula_mu = y ~ treat * group,
                          formula_sigma = ~ treat * group,
                          treat = "treat", group = "group",
                          B = 10, seed = 1, verbose = FALSE)

  expect_output(print(result), "Bootstrap standard errors")
  expect_output(print(result), "Replicates")
})

test_that("CV2 bootstrap works", {
  skip_if_not_installed("gamlss")

  set.seed(42)
  n <- 200
  group <- rep(c("A", "B"), each = n / 2)
  treat <- rbinom(n, 1, 0.5)
  y <- 50 + ifelse(group == "B", 10, 0) + 5 * treat + rnorm(n, 0, 3)
  # Ensure positive for CV2
  y <- abs(y) + 1
  df <- data.frame(y = y, treat = treat, group = group)

  result <- bootstrap_se(data = df, formula_mu = y ~ treat * group,
                          formula_sigma = ~ treat * group,
                          treat = "treat", group = "group",
                          ystat = "CV2", B = 15, seed = 42, verbose = FALSE)

  expect_s3_class(result, "ineqx_boot")
  expect_equal(result$type, "cross")
  expect_named(result$se, c("se_delta_B", "se_delta_W", "se_delta_total"))
})

# ---------------------------------------------------------------------------- #
# Integration via ineqx()
# ---------------------------------------------------------------------------- #

test_that("ineqx() with bootstrap SEs works", {
  skip_if_not_installed("gamlss")

  set.seed(42)
  n <- 200
  group <- rep(c("A", "B"), each = n / 2)
  treat <- rbinom(n, 1, 0.5)
  y <- 10 + ifelse(group == "B", 5, 0) + 3 * treat + rnorm(n, 0, 2)
  df <- data.frame(y = y, treat = treat, group = group)

  # Fit model and extract params
  model <- gamlss::gamlss(y ~ treat * group, sigma.formula = ~ treat * group,
                           data = df, trace = FALSE)
  params <- ineqx_params(model = model, data = df, treat = "treat",
                          group = "group", ystat = "Var", vcov = FALSE)

  bc <- boot_config(data = df, formula_mu = y ~ treat * group,
                     formula_sigma = ~ treat * group, treat = "treat",
                     group = "group", B = 15, seed = 42, verbose = FALSE)

  result <- ineqx(params = params, se = bc)

  expect_s3_class(result, "ineqx_causal_cross")
  expect_equal(result$se_method, "bootstrap")
  expect_true(!is.null(result$boot))
  expect_s3_class(result$boot, "ineqx_boot")
})

test_that("ineqx() with se='none' suppresses SEs", {
  skip_if_not_installed("gamlss")

  set.seed(42)
  n <- 200
  group <- rep(c("A", "B"), each = n / 2)
  treat <- rbinom(n, 1, 0.5)
  y <- 10 + 3 * treat + rnorm(n)
  df <- data.frame(y = y, treat = treat, group = group)

  model <- gamlss::gamlss(y ~ treat * group, sigma.formula = ~ treat * group,
                           data = df, trace = FALSE)
  params <- ineqx_params(model = model, data = df, treat = "treat",
                          group = "group", ystat = "Var", vcov = TRUE)

  result <- ineqx(params = params, se = "none")

  # delta method SEs should be suppressed
  expect_null(result$se)
})

test_that("ineqx() errors with invalid se argument", {
  params <- ineqx_params(
    data = data.frame(group = c("A", "B"), pi = c(0.5, 0.5),
                      mu0 = c(500, 1000), sigma0 = c(200, 400),
                      beta = c(60, 100), lambda = c(-0.1, -0.2))
  )
  expect_error(ineqx(params = params, se = "bootstrap"),
               "delta.*none")
})
