# ============================================================================ #
# Tests for DiD parameter extraction in ineqx_params (model path)
#
# Synthetic data with known beta_DP and lambda_DP coefficients.
# A correct DiD extraction must recover beta_DP as `beta` and lambda_DP as
# `lambda` for each (group, time) cell. Averaging predictions across pre+post
# rows in the cell would instead return a mixture of the cross-sectional
# dummy effect and the DiD interaction, which is the bug this test guards
# against.
# ============================================================================ #

skip_if_not_installed("gamlss")

# ---------------------------------------------------------------------------- #
# Helper: simulate one (group, period) of a 2-group, 2-period DiD GAMLSS
# ---------------------------------------------------------------------------- #
.simulate_did_gamlss <- function(n_per_cell = 2000, seed = 1) {
  set.seed(seed)

  # True parameters per group
  alpha   <- c(A = 500, B = 800)              # baseline mean (P=0, D=0)
  beta_D  <- c(A = 30,  B = 60)               # static D-on-mu (selection)
  beta_P  <- c(A = 20,  B = 40)               # time trend on mu
  beta_DP <- c(A = -50, B = -120)             # DiD ATT on mu  <-- target
  log_sig <- c(A = log(150), B = log(250))    # baseline log SD
  lam_D   <- c(A = 0.05, B = 0.05)            # static D-on-log-sigma
  lam_P   <- c(A = 0.02, B = 0.03)            # time trend on log SD
  lam_DP  <- c(A = -0.10, B = -0.20)          # DiD ATT on log SD <-- target

  rows <- list()
  for (g in c("A", "B")) {
    for (D in c(0, 1)) {
      for (P in c(0, 1)) {
        mu_cell  <- alpha[g] + beta_D[g] * D + beta_P[g] * P + beta_DP[g] * D * P
        lsd_cell <- log_sig[g] + lam_D[g] * D + lam_P[g] * P + lam_DP[g] * D * P
        sd_cell  <- exp(lsd_cell)
        y <- rnorm(n_per_cell, mean = mu_cell, sd = sd_cell)
        rows[[length(rows) + 1]] <- data.frame(
          group = g, treat = D, post = P, y = y,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  out <- do.call(rbind, rows)
  list(
    data    = out,
    beta_DP = beta_DP,
    lam_DP  = lam_DP
  )
}

# ---------------------------------------------------------------------------- #
# Test: extracted beta and lambda equal the true DiD ATT (within group)
# ---------------------------------------------------------------------------- #
test_that("ineqx_params recovers the DiD ATT (beta_DP) when post is supplied", {
  sim <- .simulate_did_gamlss(n_per_cell = 5000, seed = 42)
  d <- sim$data

  # Fit a GAMLSS that exactly matches the data-generating model:
  # mu  ~ group * treat * post
  # log sigma ~ group * treat * post
  m <- gamlss::gamlss(
    y ~ group * treat * post,
    sigma.formula = ~ group * treat * post,
    family = gamlss.dist::NO(),
    data = d,
    trace = FALSE
  )

  params <- ineqx_params(
    data  = d,
    model = m,
    treat = "treat",
    group = "group",
    post  = "post",
    vcov  = FALSE
  )

  # beta in the params data frame should match the true beta_DP per group.
  # Tolerance accounts for sampling noise at n_per_cell = 5000.
  pdata <- params$data
  pdata <- pdata[order(pdata$group), ]

  expect_equal(
    pdata$beta,
    unname(sim$beta_DP[as.character(pdata$group)]),
    tolerance = 0.05  # relative; ~5% of true effect
  )

  expect_equal(
    pdata$lambda,
    unname(sim$lam_DP[as.character(pdata$group)]),
    tolerance = 0.05
  )
})

# ---------------------------------------------------------------------------- #
# Test: simple-difference path is unchanged (no post argument)
# ---------------------------------------------------------------------------- #
test_that("ineqx_params recovers the cross-sectional ATT when post is NULL", {
  set.seed(7)
  alpha  <- c(A = 500, B = 800)
  beta_D <- c(A = -50, B = -120)  # treatment effect on mu (cross-sectional)
  log_sig <- c(A = log(150), B = log(250))
  lam_D  <- c(A = -0.10, B = -0.20)

  rows <- list()
  for (g in c("A", "B")) {
    for (D in c(0, 1)) {
      mu_cell  <- alpha[g] + beta_D[g] * D
      lsd_cell <- log_sig[g] + lam_D[g] * D
      y <- rnorm(5000, mean = mu_cell, sd = exp(lsd_cell))
      rows[[length(rows) + 1]] <- data.frame(
        group = g, treat = D, y = y, stringsAsFactors = FALSE
      )
    }
  }
  d <- do.call(rbind, rows)

  m <- gamlss::gamlss(
    y ~ group * treat,
    sigma.formula = ~ group * treat,
    family = gamlss.dist::NO(),
    data = d,
    trace = FALSE
  )

  params <- ineqx_params(
    data  = d,
    model = m,
    treat = "treat",
    group = "group",
    vcov  = FALSE
  )

  pdata <- params$data
  pdata <- pdata[order(pdata$group), ]

  expect_equal(
    pdata$beta,
    unname(beta_D[as.character(pdata$group)]),
    tolerance = 0.05
  )
  expect_equal(
    pdata$lambda,
    unname(lam_D[as.character(pdata$group)]),
    tolerance = 0.05
  )
})

# ---------------------------------------------------------------------------- #
# Test: pre-period anchor columns are populated for DiD models
# ---------------------------------------------------------------------------- #
test_that("DiD extraction stores mu1_pre/mu0_pre/sigma1_pre/sigma0_pre", {
  sim <- .simulate_did_gamlss(n_per_cell = 2000, seed = 11)
  d <- sim$data

  m <- gamlss::gamlss(
    y ~ group * treat * post,
    sigma.formula = ~ group * treat * post,
    family = gamlss.dist::NO(),
    data = d,
    trace = FALSE
  )

  params <- ineqx_params(
    data  = d,
    model = m,
    treat = "treat",
    group = "group",
    post  = "post",
    vcov  = FALSE
  )

  expect_true("mu1_pre" %in% names(params$data))
  expect_true("mu0_pre" %in% names(params$data))
  expect_true("sigma1_pre" %in% names(params$data))
  expect_true("sigma0_pre" %in% names(params$data))

  # mu0 = mu1_pre + (mu0_pre's pre->post change under control) by construction.
  # Because there is no time dimension here, "pre->post change under control"
  # collapses to mu0(post) - mu0_pre, which equals beta_P. So we should have:
  #   mu0 == mu0_pre + beta_P, where beta_P is the time trend in mu for D=0.
  # Direct check: mu1 - mu0 == beta == beta_DP.
  expect_equal(params$data$mu1 - params$data$mu0, params$data$beta,
               tolerance = 1e-6)
})

# ---------------------------------------------------------------------------- #
# Smoke test: pretrends and DiD outcome.params plots produce ggplot objects
# ---------------------------------------------------------------------------- #
test_that("DiD plot types produce ggplot objects without error", {
  skip_if_not_installed("ggplot2")

  # Build a small longitudinal DiD dataset across two calendar periods
  set.seed(101)
  alpha   <- c(A = 500, B = 800)
  beta_D  <- c(A = 30,  B = 60)
  beta_P  <- c(A = 20,  B = 40)
  beta_DP <- c(A = -50, B = -120)
  log_sig <- c(A = log(150), B = log(250))
  lam_D   <- c(A = 0.05, B = 0.05)
  lam_P   <- c(A = 0.02, B = 0.03)
  lam_DP  <- c(A = -0.10, B = -0.20)

  rows <- list()
  for (period_yr in c(2000, 2010)) {
    for (g in c("A", "B")) {
      for (D in c(0, 1)) {
        for (P in c(0, 1)) {
          mu_cell  <- alpha[g] + beta_D[g] * D + beta_P[g] * P + beta_DP[g] * D * P
          lsd_cell <- log_sig[g] + lam_D[g] * D + lam_P[g] * P + lam_DP[g] * D * P
          y <- rnorm(800, mean = mu_cell, sd = exp(lsd_cell))
          rows[[length(rows) + 1]] <- data.frame(
            group = g, treat = D, post = P, year = period_yr, y = y,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  d <- do.call(rbind, rows)

  res <- ineqx(
    y     = "y",
    ystat = "Var",
    treat = "treat",
    post  = "post",
    group = "group",
    time  = "year",
    ref   = 2000,
    formula_mu    = ~ group * treat * post * factor(year),
    formula_sigma = ~ group * treat * post * factor(year),
    data  = d,
    se    = "none"
  )

  expect_true(isTRUE(res$params$is_did))
  expect_equal(res$params$post, "post")

  # outcome.params (longit) should produce a ggplot
  p1 <- plot(res, type = "outcome.params")
  expect_s3_class(p1, "ggplot")

  # pretrends should produce a ggplot for DiD models
  p2 <- plot(res, type = "pretrends")
  expect_s3_class(p2, "ggplot")
})

# ---------------------------------------------------------------------------- #
# Test: delta-method SEs are produced for DiD models
# ---------------------------------------------------------------------------- #
test_that("DiD delta-method vcov is finite and produces sensible SEs", {
  sim <- .simulate_did_gamlss(n_per_cell = 3000, seed = 23)
  d <- sim$data

  m <- gamlss::gamlss(
    y ~ group * treat * post,
    sigma.formula = ~ group * treat * post,
    family = gamlss.dist::NO(),
    data = d,
    trace = FALSE
  )

  params <- ineqx_params(
    data  = d, model = m,
    treat = "treat", group = "group", post = "post",
    vcov  = TRUE
  )

  expect_false(is.null(params$vcov))
  V <- params$vcov
  expect_true(all(is.finite(diag(V))))
  expect_true(all(diag(V) > 0))

  # SE on beta (the ATT) should be in a sensible range — order of 1-50 dollars
  # for sample size 3000 per cell, given the data-generating process.
  J <- params$n_groups
  se_beta <- sqrt(diag(V)[1:J])
  expect_true(all(se_beta > 0 & se_beta < 50))
})

# ---------------------------------------------------------------------------- #
# Test: is_did flag is set correctly
# ---------------------------------------------------------------------------- #
test_that("ineqx_params sets is_did based on whether post is supplied", {
  sim <- .simulate_did_gamlss(n_per_cell = 1000, seed = 13)
  d <- sim$data

  m <- gamlss::gamlss(
    y ~ group * treat * post,
    sigma.formula = ~ group * treat * post,
    family = gamlss.dist::NO(),
    data = d,
    trace = FALSE
  )

  # With post: is_did TRUE
  params_did <- ineqx_params(
    data  = d, model = m,
    treat = "treat", group = "group", post = "post",
    vcov  = FALSE
  )
  expect_true(isTRUE(params_did$is_did))
  expect_equal(params_did$post, "post")

  # Same model, no post argument: is_did FALSE (treated as simple difference)
  params_sd <- ineqx_params(
    data  = d, model = m,
    treat = "treat", group = "group",
    vcov  = FALSE
  )
  expect_false(isTRUE(params_sd$is_did))
  expect_null(params_sd$post)
})
