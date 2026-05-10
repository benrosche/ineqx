# ============================================================================ #
# Tests for descriptive decomposition
# ============================================================================ #

# ---------------------------------------------------------------------------- #
# ineqx_params: auto-detection of descriptive vs causal
# ---------------------------------------------------------------------------- #

test_that("ineqx_params with descriptive columns returns ineqx_desc_params", {
  ref <- ineqx_params(
    data = data.frame(
      group = c("A", "B"),
      pi = c(0.5, 0.5),
      mu = c(100, 200),
      sigma = c(30, 50)
    )
  )
  expect_s3_class(ref, "ineqx_desc_params")
  expect_equal(ref$n_groups, 2)
  expect_equal(ref$groups, c("A", "B"))
  expect_equal(ref$ystat, "Var")
})

test_that("ineqx_params with causal columns still returns ineqx_params", {
  params <- ineqx_params(
    data = data.frame(
      group = c("A", "B"),
      pi = c(0.5, 0.5),
      mu0 = c(100, 200),
      sigma0 = c(30, 50),
      beta = c(10, 20),
      lambda = c(-0.1, -0.2)
    )
  )
  expect_s3_class(params, "ineqx_params")
})

test_that("ineqx_desc_params accepts time column", {
  ref <- ineqx_params(data = data.frame(
    group = c("A", "B"), time = c(1, 1),
    pi = c(0.5, 0.5), mu = c(100, 200), sigma = c(30, 50)
  ))
  expect_s3_class(ref, "ineqx_desc_params")
  expect_equal(ref$times, 1)
})

test_that("ineqx_desc_params validates pi sums to 1", {
  expect_error(
    ineqx_params(data = data.frame(
      group = c("A", "B"), pi = c(0.3, 0.3),
      mu = c(100, 200), sigma = c(30, 50)
    )),
    "pi.*sum"
  )
})

test_that("ineqx_desc_params validates sigma >= 0", {
  expect_error(
    ineqx_params(data = data.frame(
      group = c("A", "B"), pi = c(0.5, 0.5),
      mu = c(100, 200), sigma = c(-1, 50)
    )),
    "sigma.*non-negative"
  )
})

# ---------------------------------------------------------------------------- #
# Descriptive with counterfactual reference
# ---------------------------------------------------------------------------- #

test_that("descriptive with zero counterfactual ref works", {
  set.seed(42)
  n <- 200
  d <- data.frame(
    y = c(rnorm(n, 100, 30), rnorm(n, 200, 50)),
    group = rep(c("A", "B"), each = n)
  )

  ref <- ineqx_params(data = data.frame(
    group = c("A", "B"), pi = c(0.5, 0.5), mu = c(0, 0), sigma = c(0, 0)
  ))

  result <- ineqx("y", group = "group", params = ref, data = d)

  expect_s3_class(result, "ineqx_desc")
  expect_false(is.null(result$deltas))
  expect_false(is.null(result$ref_params))

  # Ref period totals should be zero
  ref_row <- result$totals[result$totals$time == result$ref, ]
  expect_equal(ref_row$VarT, 0)
})

test_that("counterfactual ref deltas sum to observed inequality", {
  set.seed(42)
  n <- 200
  d <- data.frame(
    y = c(rnorm(n, 100, 30), rnorm(n, 200, 50)),
    group = rep(c("A", "B"), each = n)
  )

  ref <- ineqx_params(data = data.frame(
    group = c("A", "B"), pi = c(0.5, 0.5), mu = c(0, 0), sigma = c(0, 0)
  ))

  result <- ineqx("y", group = "group", params = ref, data = d)

  # Observed time period
  obs_time <- setdiff(result$totals$time, result$ref)
  obs_row <- result$totals[result$totals$time == obs_time, ]
  delta_row <- result$deltas[result$deltas$time == obs_time, ]

  # delta_T should equal observed total inequality
  expect_equal(delta_row$delta_T, obs_row$VarT, tolerance = 1e-10)
})

test_that("counterfactual ref with multiple time periods", {
  data(incdat)
  ref <- ineqx_params(data = data.frame(
    group = 1:3, pi = 1/3, mu = 0, sigma = 0
  ))

  result <- ineqx("inc", group = "group", time = "year",
                   params = ref, data = incdat)

  expect_s3_class(result, "ineqx_desc")

  # Should have deltas for each time (ref + observed years)
  obs_years <- sort(unique(incdat$year))

  # Deltas should sum correctly for each year
  for (yr in obs_years) {
    d <- result$deltas[result$deltas$time == yr, ]
    expect_equal(d$delta_mu + d$delta_sigma + d$delta_pi,
                 d$delta_T, tolerance = 1e-10)
    expect_equal(d$delta_mu_W + d$delta_sigma_W + d$delta_pi_W,
                 d$delta_W, tolerance = 1e-10)
    expect_equal(d$delta_mu_B + d$delta_sigma_B + d$delta_pi_B,
                 d$delta_B, tolerance = 1e-10)
  }
})

test_that("counterfactual ref with specific ordering", {
  set.seed(42)
  n <- 200
  d <- data.frame(
    y = c(rnorm(n, 100, 30), rnorm(n, 200, 50), rnorm(n, 300, 80)),
    group = rep(c("A", "B", "C"), each = n)
  )

  ref <- ineqx_params(data = data.frame(
    group = c("A", "B", "C"), pi = 1/3, mu = 0, sigma = 0
  ))

  result <- ineqx("y", group = "group", params = ref,
                   order = c("mu", "sigma", "pi"), data = d)

  obs_time <- setdiff(result$totals$time, result$ref)
  d_row <- result$deltas[result$deltas$time == obs_time, ]
  expect_equal(d_row$delta_mu + d_row$delta_sigma + d_row$delta_pi,
               d_row$delta_T, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------- #
# Validation: counterfactual ref
# ---------------------------------------------------------------------------- #

test_that("desc params + treat errors", {
  ref <- ineqx_params(data = data.frame(
    group = c("A", "B"), pi = c(0.5, 0.5), mu = c(0, 0), sigma = c(0, 0)
  ))
  d <- data.frame(y = 1:10, group = rep(c("A", "B"), 5), treat = rep(0:1, 5))

  expect_error(
    ineqx("y", treat = "treat", group = "group", params = ref, data = d),
    "descriptive params.*treat"
  )
})

test_that("desc params with mismatched groups errors", {
  ref <- ineqx_params(data = data.frame(
    group = c("X", "Y"), pi = c(0.5, 0.5), mu = c(0, 0), sigma = c(0, 0)
  ))
  d <- data.frame(y = 1:10, group = rep(c("A", "B"), 5))

  expect_error(
    ineqx("y", group = "group", params = ref, data = d),
    "Groups.*do not match"
  )
})

test_that("desc params without y errors", {
  ref <- ineqx_params(data = data.frame(
    group = c("A", "B"), pi = c(0.5, 0.5), mu = c(0, 0), sigma = c(0, 0)
  ))
  d <- data.frame(y = 1:10, group = rep(c("A", "B"), 5))

  expect_error(
    ineqx(group = "group", params = ref, data = d),
    "y.*required"
  )
})

# ---------------------------------------------------------------------------- #
# Raw data descriptive (no params, no ref) still works
# ---------------------------------------------------------------------------- #

test_that("raw data descriptive without ref works", {
  set.seed(42)
  n <- 200
  d <- data.frame(
    y = c(rnorm(n, 100, 30), rnorm(n, 200, 50)),
    group = rep(c("A", "B"), each = n)
  )

  result <- ineqx("y", group = "group", data = d)
  expect_s3_class(result, "ineqx_desc")
  expect_null(result$deltas)
  expect_true(result$totals$VarT > 0)
})

test_that("raw data descriptive with ref works", {
  data(incdat)
  result <- ineqx("inc", group = "group", time = "year",
                   ref = 1, data = incdat)

  expect_s3_class(result, "ineqx_desc")
  expect_false(is.null(result$deltas))
  expect_equal(result$ref, 1)
})

test_that("y=NULL without params gives informative error", {
  d <- data.frame(
    group = c("A", "B"), pi = c(0.5, 0.5),
    mu = c(100, 200), sigma = c(30, 50)
  )
  expect_error(
    ineqx(group = "group", data = d),
    "y.*required.*params"
  )
})

# ---------------------------------------------------------------------------- #
# Multi-period descriptive params (blending)
# ---------------------------------------------------------------------------- #

test_that("ineqx_desc_params accepts optional time column", {
  ref <- ineqx_params(data = data.frame(
    group = rep(c("A", "B"), 2),
    time = rep(c(0, 5), each = 2),
    pi = c(0.5, 0.5, 0.6, 0.4),
    mu = c(0, 0, 50, 100),
    sigma = c(0, 0, 10, 20)
  ))

  expect_s3_class(ref, "ineqx_desc_params")
  expect_equal(ref$times, c(0, 5))
  expect_equal(ref$n_times, 2)
})

test_that("multi-period params blend with observed data (no overlap)", {
  data(incdat)
  # Params at time=0 (not in data), data has years 1-5
  ref <- ineqx_params(data = data.frame(
    group = 1:3, pi = 1/3, mu = 0, sigma = 0
  ))

  result <- ineqx("inc", group = "group", time = "year",
                   params = ref, data = incdat)

  # Should have 6 time levels: 0 (counterfactual) + 1,2,3,4,5 (observed)
  expect_equal(sort(unique(result$wibe$time)), c(0, 1, 2, 3, 4, 5))
  expect_equal(result$ref, 0L)
})

test_that("multi-period params override observed periods (overlap)", {
  set.seed(42)
  n <- 100
  d <- data.frame(
    y = c(rnorm(n, 100, 30), rnorm(n, 200, 50),
          rnorm(n, 150, 40), rnorm(n, 250, 60)),
    group = rep(c("A", "B"), 2 * n),
    year = rep(c(1, 2), each = 2 * n)
  )

  # Params provides counterfactual for time=0 AND overrides year=1
  ref <- ineqx_params(data = data.frame(
    group = rep(c("A", "B"), 2),
    time = rep(c(0, 1), each = 2),
    pi = c(0.5, 0.5, 0.5, 0.5),
    mu = c(0, 0, 90, 190),
    sigma = c(0, 0, 25, 45)
  ))

  result <- ineqx("y", group = "group", time = "year",
                   params = ref, ref = 0, data = d)

  # Year=1 should use params values, not estimated from data
  yr1_wibe <- result$wibe[result$wibe$time == 1, ]
  expect_equal(yr1_wibe$mu, c(90, 190))
  expect_equal(yr1_wibe$sigma, c(25, 45))

  # Year=2 should be estimated from data (not in params)
  yr2_wibe <- result$wibe[result$wibe$time == 2, ]
  expect_true(all(yr2_wibe$mu != 0))  # non-zero, estimated
})

test_that("multi-period params with ref pointing to params period", {
  set.seed(42)
  n <- 200
  d <- data.frame(
    y = c(rnorm(n, 100, 30), rnorm(n, 200, 50)),
    group = rep(c("A", "B"), each = n)
  )

  ref <- ineqx_params(data = data.frame(
    group = c("A", "B"), time = c(0, 0),
    pi = c(0.5, 0.5), mu = c(0, 0), sigma = c(0, 0)
  ))

  result <- ineqx("y", group = "group", params = ref, ref = 0, data = d)

  expect_equal(result$ref, 0)
  expect_false(is.null(result$deltas))
})

test_that("multi-period params with ref pointing to observed period", {
  data(incdat)

  # Params provides a hypothetical scenario at time=99
  ref <- ineqx_params(data = data.frame(
    group = 1:3,
    time = rep(99, 3),
    pi = c(0.5, 0.3, 0.2),
    mu = c(500, 800, 1500),
    sigma = c(100, 200, 400)
  ))

  result <- ineqx("inc", group = "group", time = "year",
                   params = ref, ref = 1, data = incdat)

  # ref=1 points to observed year 1
  expect_equal(result$ref, 1)
  # Time=99 should be present
  expect_true(99 %in% result$totals$time)
})

test_that("multi-period params without ref errors", {
  ref <- ineqx_params(data = data.frame(
    group = rep(c("A", "B"), 2),
    time = rep(c(0, 5), each = 2),
    pi = c(0.5, 0.5, 0.6, 0.4),
    mu = c(0, 0, 50, 100),
    sigma = c(0, 0, 10, 20)
  ))

  d <- data.frame(y = rnorm(100), group = rep(c("A", "B"), 50))

  expect_error(
    ineqx("y", group = "group", params = ref, data = d),
    "ref.*required.*multiple time"
  )
})

# ---------------------------------------------------------------------------- #
# Plot decomp CIs
# ---------------------------------------------------------------------------- #

test_that("plot decomp with ci='delta' produces CI ribbons", {
  data(incdat)
  result <- ineqx("inc", group = "group", time = "year",
                   ref = 1, se = "delta", data = incdat)

  p <- plot(result, type = "decomp", ci = "delta")
  expect_s3_class(p, "ggplot")

  # Check that ymin/ymax are in the plot data
  plot_data <- ggplot2::ggplot_build(p)$data
  # At least one layer should have ymin/ymax
  has_ci <- any(vapply(plot_data, function(d) "ymin" %in% names(d), logical(1)))
  expect_true(has_ci)
})

test_that("plot decomp with ci=TRUE uses stored delta SEs", {
  data(incdat)
  result <- ineqx("inc", group = "group", time = "year",
                   ref = 1, se = "delta", data = incdat)

  p <- plot(result, type = "decomp", ci = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that(".resolve_ci accepts 'boot' string", {
  # Minimal mock x with no stored SEs
  x <- list(se = NULL)
  ci_info <- ineqx:::.resolve_ci("boot", x)
  expect_equal(ci_info$method, "boot")
  expect_s3_class(ci_info$boot, "ineqx_boot_config")
})

test_that("deltas sum correctly with blended params", {
  data(incdat)
  ref <- ineqx_params(data = data.frame(
    group = 1:3, pi = 1/3, mu = 0, sigma = 0
  ))

  result <- ineqx("inc", group = "group", time = "year",
                   params = ref, data = incdat)

  for (yr in unique(incdat$year)) {
    d <- result$deltas[result$deltas$time == yr, ]
    expect_equal(d$delta_mu + d$delta_sigma + d$delta_pi,
                 d$delta_T, tolerance = 1e-10)
  }
})
