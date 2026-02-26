# ============================================================================ #
# Tests for ineqx_params
# ============================================================================ #

test_that("ineqx_params creates valid cross-sectional object", {
  params <- ineqx_params(
    data = data.frame(
      group = c("A", "B"),
      pi = c(0.5, 0.5),
      mu0 = c(500, 1000),
      sigma0 = c(200, 400),
      beta = c(60, 100),
      lambda = c(-0.1, -0.2)
    )
  )

  expect_s3_class(params, "ineqx_params")
  expect_equal(params$type, "cross_sectional")
  expect_equal(params$n_groups, 2)
  expect_null(params$times)
})

test_that("ineqx_params creates valid longitudinal object", {
  params <- ineqx_params(
    data = data.frame(
      group = rep(c("A", "B"), 2),
      time = rep(c(2000, 2010), each = 2),
      pi = rep(c(0.5, 0.5), 2),
      mu0 = c(500, 1000, 550, 1100),
      sigma0 = c(200, 400, 210, 380),
      beta = c(60, 100, 80, 90),
      lambda = c(-0.1, -0.2, -0.15, -0.15)
    ),
    ref = 2000
  )

  expect_s3_class(params, "ineqx_params")
  expect_equal(params$type, "longitudinal")
  expect_equal(params$ref, 2000)
  expect_equal(params$n_times, 2)
})

test_that("ineqx_params validates pi sums to 1", {
  expect_error(
    ineqx_params(data = data.frame(
      group = c("A", "B"), pi = c(0.3, 0.5),
      mu0 = c(1, 2), sigma0 = c(1, 1), beta = c(0, 0), lambda = c(0, 0)
    )),
    "sum to 1"
  )
})

test_that("ineqx_params validates sigma0 > 0", {
  expect_error(
    ineqx_params(data = data.frame(
      group = c("A", "B"), pi = c(0.5, 0.5),
      mu0 = c(1, 2), sigma0 = c(-1, 1), beta = c(0, 0), lambda = c(0, 0)
    )),
    "strictly positive"
  )
})

test_that("ineqx_params requires ref for longitudinal", {
  expect_error(
    ineqx_params(data = data.frame(
      group = rep(c("A", "B"), 2), time = rep(c(1, 2), each = 2),
      pi = rep(c(0.5, 0.5), 2), mu0 = 1:4, sigma0 = rep(1, 4),
      beta = rep(0, 4), lambda = rep(0, 4)
    )),
    "ref.*required"
  )
})

test_that("ineqx_params rejects missing columns", {
  expect_error(
    ineqx_params(data = data.frame(group = c("A", "B"), pi = c(0.5, 0.5))),
    "Missing required"
  )
})
