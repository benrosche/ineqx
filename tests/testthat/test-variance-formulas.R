# ============================================================================ #
# Tests for variance formula functions
# ============================================================================ #

test_that("VarW_pi computes weighted within-group variance", {
  # Two groups, equal weights
  pi <- c(0.5, 0.5)
  sigma <- c(200, 400)
  expected <- 0.5 * 200^2 + 0.5 * 400^2  # 100,000
  expect_equal(VarW_pi(pi, sigma), expected)
})

test_that("VarB_pi computes weighted between-group variance", {
  pi <- c(0.5, 0.5)
  mu <- c(500, 1000)
  gmu <- 750
  expected <- 0.5 * (500 - 750)^2 + 0.5 * (1000 - 750)^2  # 62,500
  expect_equal(VarB_pi(pi, mu), expected)
})

test_that("VarT_pi = VarW_pi + VarB_pi", {
  pi <- c(0.3, 0.5, 0.2)
  mu <- c(100, 200, 300)
  sigma <- c(10, 20, 30)
  expect_equal(VarT_pi(pi, mu, sigma), VarW_pi(pi, sigma) + VarB_pi(pi, mu))
})

test_that("CV2 decomposition sums correctly", {
  pi <- c(0.4, 0.6)
  mu <- c(500, 1000)
  sigma <- c(100, 200)
  expect_equal(CV2T_pi(pi, mu, sigma),
               CV2W_pi(pi, mu, sigma) + CV2B_pi(pi, mu))
})

test_that("count-based and pi-based formulas agree", {
  n <- c(30, 70)
  pi <- n / sum(n)
  mu <- c(500, 1000)
  sigma <- c(200, 400)

  expect_equal(VarW_n(n, sigma), VarW_pi(pi, sigma))
  expect_equal(VarB_n(n, mu), VarB_pi(pi, mu))
  expect_equal(CV2W_n(n, mu, sigma), CV2W_pi(pi, mu, sigma))
  expect_equal(CV2B_n(n, mu), CV2B_pi(pi, mu))
})
