# Tests for the beta-homogeneity Wald test (beta_test / beta_tests) and the
# proportionality / share plot additions. Uses manual params with an injected
# diagonal vcov so the linear-algebra is exercised without fitting GAMLSS.

# vcov is ordered (beta, mu0, lambda, log sigma0); beta occupies the first J slots.
.inject_vcov <- function(J, var_beta = 100) {
  diag(rep(c(var_beta, 1, 0.01, 0.01), each = J))
}

.cross_with_vcov <- function(beta, transform = NULL) {
  p <- ineqx_params(data = data.frame(
    group  = c("A", "B", "C"), pi = rep(1 / 3, 3),
    mu0    = c(500, 800, 1200), sigma0 = c(200, 300, 400),
    beta   = beta, lambda = c(0, 0, 0)))
  p$vcov <- .inject_vcov(3)
  if (!is.null(transform)) p$transform <- transform
  ineqx("y", group = "g", data = data.frame(y = 1, g = "A"),
        params = p, se = "delta")
}

.longit_with_vcov <- function() {
  grid <- expand.grid(group = c("A", "B", "C"), time = c(2000, 2010, 2020),
                      stringsAsFactors = FALSE)
  bm <- c(A = 500, B = 800, C = 1200); bs <- c(A = 200, B = 300, C = 400)
  grid$pi <- 1 / 3
  grid$mu0 <- bm[grid$group]; grid$sigma0 <- bs[grid$group]
  grid$beta <- -0.08 * grid$mu0; grid$lambda <- 0      # proportional effect
  p <- ineqx_params(data = grid)
  p$vcov <- stats::setNames(
    lapply(unique(grid$time), function(t) .inject_vcov(3)),
    as.character(unique(grid$time)))
  ineqx(params = p, ref = 2000, se = "delta")
}

test_that("beta_test: equal betas are not rejected, df = J-1", {
  r <- .cross_with_vcov(c(-50, -50, -50))
  expect_false(is.null(r$beta_test))
  expect_equal(r$beta_test$df, 2L)
  expect_equal(r$beta_test$statistic, 0, tolerance = 1e-6)
  expect_gt(r$beta_test$p_value, 0.99)
})

test_that("beta_test: unequal betas are rejected", {
  r <- .cross_with_vcov(c(-20, -50, -120))
  expect_equal(r$beta_test$df, 2L)
  expect_lt(r$beta_test$p_value, 0.001)
})

test_that("beta_test: scale label follows params$transform", {
  expect_equal(.cross_with_vcov(c(-50, -50, -50))$beta_test$scale, "identity")
  expect_equal(.cross_with_vcov(c(-50, -50, -50), transform = "log")$beta_test$scale,
               "log")
})

test_that("beta_test is NULL without vcov and with a single group", {
  # No vcov -> no test
  p <- ineqx_params(data = data.frame(
    group = c("A", "B"), pi = c(.5, .5), mu0 = c(500, 800),
    sigma0 = c(200, 300), beta = c(-50, -80), lambda = c(0, 0)))
  r <- ineqx("y", group = "g", data = data.frame(y = 1, g = "A"),
             params = p, se = "none")
  expect_null(r$beta_test)
})

test_that("beta_tests: longitudinal returns a per-time list", {
  r <- .longit_with_vcov()
  expect_false(is.null(r$beta_tests))
  expect_setequal(names(r$beta_tests), c("2000", "2010", "2020"))
  expect_equal(r$beta_tests[["2010"]]$df, 2L)
  # proportional beta on the level scale -> betas differ -> rejected
  expect_lt(r$beta_tests[["2010"]]$p_value, 0.05)
})

test_that("print() shows the beta-homogeneity Wald line", {
  r <- .cross_with_vcov(c(-20, -50, -120))
  out <- paste(capture.output(print(r)), collapse = "\n")
  expect_match(out, "beta homogeneous")
})

# --- plot additions -------------------------------------------------------- #

test_that("type = 'effect.prop' returns a ggplot (longit + cross)", {
  skip_if_not_installed("ggplot2")
  rl <- .longit_with_vcov()
  expect_s3_class(plot(rl, type = "effect.prop", ci = TRUE), "ggplot")
  rc <- .cross_with_vcov(c(-20, -50, -120))
  expect_s3_class(plot(rc, type = "effect.prop", ci = TRUE), "ggplot")
})

test_that("share = TRUE returns a ggplot for causal decomp and cross wibe", {
  skip_if_not_installed("ggplot2")
  rl <- .longit_with_vcov()
  expect_s3_class(plot(rl, type = "decomp", share = TRUE), "ggplot")
  rc <- .cross_with_vcov(c(-20, -50, -120))
  expect_s3_class(plot(rc, type = "wibe", share = TRUE), "ggplot")
})

test_that("descriptive wibe share returns a ggplot", {
  skip_if_not_installed("ggplot2")
  set.seed(1); n <- 200
  bm <- c(A = 500, B = 800, C = 1200)
  raw <- do.call(rbind, lapply(c(2000, 2010, 2020), function(tt)
    do.call(rbind, lapply(c("A", "B", "C"), function(g)
      data.frame(time = tt, grp = g,
                 y = rlnorm(n, meanlog = log(bm[g]), sdlog = 0.4))))))
  desc <- ineqx("y", group = "grp", time = "time", data = raw, se = "none")
  expect_s3_class(plot(desc, type = "wibe", share = TRUE), "ggplot")
})
