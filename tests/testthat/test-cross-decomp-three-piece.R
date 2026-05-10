# ============================================================================ #
# Tests for the three-piece cross-sectional het / cov / rescale decomposition
# ============================================================================ #

# Shared fixture: workers / managers, two groups, treatment lowers means and
# compresses dispersion in both groups, with motherhood penalties differing
# across groups so all three sub-components are non-zero. The explicit ystat=
# in ineqx() now overrides params$ystat, so the fixture itself doesn't need
# to fix a scale.
make_params <- function() {
  ineqx_params(
    data = data.frame(
      group  = c("workers", "managers"),
      pi     = c(0.5, 0.5),
      mu0    = c(500, 1000),
      sigma0 = c(200, 400),
      beta   = c(60,  100),     # workers + managers earn more if treated, but..
      lambda = c(-0.1, -0.2)    # treatment compresses dispersion
    )
  )
}


# ---- Algebraic identity under V -------------------------------------------- #

test_that("three-piece decomposition: rescale is zero under V", {
  params <- make_params()
  res <- ineqx("y", group = "g", data = data.frame(y = 1, g = "A"),
               params = params, ystat = "Var", se = "none")

  expect_equal(res$components$rescale_B, 0)
  expect_equal(res$components$rescale_W, 0)
})

test_that("three-piece decomposition sums to tau under V", {
  params <- make_params()
  res <- ineqx("y", group = "g", data = data.frame(y = 1, g = "A"),
               params = params, ystat = "Var", se = "none")

  cb <- res$components
  expect_equal(cb$het_B + cb$cov_B + cb$rescale_B, res$tau_B)
  expect_equal(cb$het_W + cb$cov_W + cb$rescale_W, res$tau_W)
})

test_that("under V, het and cov match the V-case closed forms", {
  params <- make_params()
  res <- ineqx("y", group = "g", data = data.frame(y = 1, g = "A"),
               params = params, ystat = "Var", se = "none")
  cb <- res$components

  pi    <- c(0.5, 0.5)
  beta  <- c(60, 100)
  bbar  <- sum(pi * beta)
  Vbeta <- sum(pi * (beta - bbar)^2)
  expect_equal(cb$het_B, Vbeta)

  mu0    <- c(500, 1000)
  mbar0  <- sum(pi * mu0)
  Cmubeta <- sum(pi * (mu0 - mbar0) * (beta - bbar))
  expect_equal(cb$cov_B, 2 * Cmubeta)
})


# ---- Algebraic identity under CV² ------------------------------------------ #

test_that("three-piece decomposition sums to tau under CV2", {
  params <- make_params()
  res <- ineqx("y", group = "g", data = data.frame(y = 1, g = "A"),
               params = params, ystat = "CV2", se = "none")

  cb <- res$components
  expect_equal(cb$het_B + cb$cov_B + cb$rescale_B, res$tau_B)
  expect_equal(cb$het_W + cb$cov_W + cb$rescale_W, res$tau_W)
})

test_that("under CV2, het and cov are V-case forms scaled by d_1", {
  params <- make_params()
  res_v   <- ineqx("y", group = "g", data = data.frame(y = 1, g = "A"),
                   params = params, ystat = "Var", se = "none")
  res_cv2 <- ineqx("y", group = "g", data = data.frame(y = 1, g = "A"),
                   params = params, ystat = "CV2", se = "none")

  pi    <- c(0.5, 0.5)
  mu0   <- c(500, 1000)
  beta  <- c(60, 100)
  d_1   <- 1 / sum(pi * (mu0 + beta))^2

  expect_equal(res_cv2$components$het_B, d_1 * res_v$components$het_B)
  expect_equal(res_cv2$components$cov_B, d_1 * res_v$components$cov_B)
  expect_equal(res_cv2$components$het_W, d_1 * res_v$components$het_W)
  expect_equal(res_cv2$components$cov_W, d_1 * res_v$components$cov_W)
})

test_that("under CV2, rescale is non-zero and matches B_0*(d_1-d_0)", {
  params <- make_params()
  res <- ineqx("y", group = "g", data = data.frame(y = 1, g = "A"),
               params = params, ystat = "CV2", se = "none")

  pi     <- c(0.5, 0.5)
  mu0    <- c(500, 1000)
  beta   <- c(60, 100)
  sigma0 <- c(200, 400)
  mbar0  <- sum(pi * mu0)
  mbar1  <- sum(pi * (mu0 + beta))
  d_0    <- 1 / mbar0^2
  d_1    <- 1 / mbar1^2
  B_0    <- sum(pi * (mu0 - mbar0)^2)
  W_0    <- sum(pi * sigma0^2)

  expect_equal(res$components$rescale_B, B_0 * (d_1 - d_0))
  expect_equal(res$components$rescale_W, W_0 * (d_1 - d_0))
  expect_true(abs(res$components$rescale_B) > 0)
  expect_true(abs(res$components$rescale_W) > 0)
})
