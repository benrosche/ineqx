# ============================================================================ #
# Tests for the unified `stats` vocabulary in plot(..., type = "wibe")
# ============================================================================ #

# ---- Normalizer unit tests ------------------------------------------------- #

test_that("normalizer returns default when stats is NULL", {
  expect_equal(.normalize_wibe_stats(NULL),
               c("tau", "tau_b", "tau_w"))
  expect_equal(.normalize_wibe_stats(NULL, default = "tau"),
               "tau")
})

test_that("normalizer expands the three het_cov shorthands", {
  expect_equal(.normalize_wibe_stats("het_cov"),
               c("het_b", "cov_b", "het_w", "cov_w"))
  expect_equal(.normalize_wibe_stats("het_cov_b"),
               c("het_b", "cov_b"))
  expect_equal(.normalize_wibe_stats("het_cov_w"),
               c("het_w", "cov_w"))
})

test_that("normalizer expands the three subs shorthands (incl. rescale)", {
  expect_equal(.normalize_wibe_stats("subs"),
               c("het_b", "cov_b", "rescale_b",
                 "het_w", "cov_w", "rescale_w"))
  expect_equal(.normalize_wibe_stats("subs_b"),
               c("het_b", "cov_b", "rescale_b"))
  expect_equal(.normalize_wibe_stats("subs_w"),
               c("het_w", "cov_w", "rescale_w"))
})

test_that("normalizer accepts canonical rescale names", {
  expect_equal(.normalize_wibe_stats("rescale_b"), "rescale_b")
  expect_equal(.normalize_wibe_stats("rescale_w"), "rescale_w")
})

test_that("normalizer mixes shorthands with canonical names and dedupes", {
  out <- .normalize_wibe_stats(c("tau", "het_cov"))
  expect_equal(out, c("tau", "het_b", "cov_b", "het_w", "cov_w"))

  # Duplicate after expansion should be unique-d, original ordering preserved
  out2 <- .normalize_wibe_stats(c("het_b", "het_cov_b"))
  expect_equal(out2, c("het_b", "cov_b"))
})

test_that("normalizer rejects unknown stats", {
  expect_error(.normalize_wibe_stats("foo"),
               "Invalid 'stats' values: foo")
  expect_error(.normalize_wibe_stats(c("tau", "bar")),
               "Invalid 'stats' values: bar")
})


# ---- Fixtures -------------------------------------------------------------- #

cross_fixture <- function() {
  params <- ineqx_params(
    data = data.frame(
      group = c("workers", "managers"),
      pi = c(0.5, 0.5),
      mu0 = c(500, 1000),
      sigma0 = c(200, 400),
      beta = c(60, 100),
      lambda = c(-0.1, -0.2)
    )
  )
  ineqx("y", group = "g", data = data.frame(y = 1, g = "A"),
        params = params, se = "none")
}

longit_fixture <- function() {
  params <- ineqx_params(
    data = data.frame(
      group = rep(c("A", "B"), 2),
      time = rep(c(2000, 2010), each = 2),
      pi = c(0.4, 0.6, 0.45, 0.55),
      mu0 = c(500, 1000, 550, 1100),
      sigma0 = c(200, 400, 210, 380),
      beta = c(60, 100, 80, 90),
      lambda = c(-0.1, -0.2, -0.15, -0.15)
    )
  )
  ineqx("y", group = "g", data = data.frame(y = 1, g = "A"),
        params = params, ref = 2000,
        order = c("behavioral", "compositional", "pretreatment"),
        se = "none")
}


# ---- Cross-sectional smoke tests ------------------------------------------- #

test_that("cross plot accepts every canonical stat individually", {
  cs <- cross_fixture()
  for (s in c("tau", "tau_b", "tau_w",
              "het_b", "cov_b", "rescale_b",
              "het_w", "cov_w", "rescale_w")) {
    p <- plot(cs, type = "wibe", stats = s)
    expect_s3_class(p, "ggplot")
  }
})

test_that("cross plot accepts every shorthand", {
  cs <- cross_fixture()
  for (s in c("het_cov", "het_cov_b", "het_cov_w",
              "subs", "subs_b", "subs_w")) {
    p <- plot(cs, type = "wibe", stats = s)
    expect_s3_class(p, "ggplot")
  }
})

test_that("cross plot accepts mixed canonical + shorthand", {
  cs <- cross_fixture()
  p <- plot(cs, type = "wibe", stats = c("tau_b", "het_b", "cov_b"))
  expect_s3_class(p, "ggplot")
})

test_that("cross plot rejects unknown stats", {
  cs <- cross_fixture()
  expect_error(plot(cs, type = "wibe", stats = "foo"),
               "Invalid 'stats' values: foo")
})


# ---- Cross-sectional backward-compat ----------------------------------------#

test_that("legacy stats = 'tau' renders three τ bars (Total/Between/Within)", {
  cs <- cross_fixture()
  p <- plot(cs, type = "wibe", stats = "tau")
  ld <- ggplot2::layer_data(p, 1L)
  expect_equal(nrow(ld), 3L)
  expect_setequal(as.character(ld$x), c("1", "2", "3"))
  vals_by_x <- setNames(ld$y, as.character(ld$x))
  expect_equal(unname(vals_by_x["1"]), cs$tau_total)
  expect_equal(unname(vals_by_x["2"]), cs$tau_B)
  expect_equal(unname(vals_by_x["3"]), cs$tau_W)
})

test_that("legacy stats = 'het_cov' renders three stacked sub bars", {
  cs <- cross_fixture()
  p <- plot(cs, type = "wibe", stats = "het_cov")
  ld <- ggplot2::layer_data(p, 1L)
  # 2 segments (het + cov) at each of x = 1, 2, 3
  expect_equal(nrow(ld), 6L)
  expect_setequal(as.character(round(ld$x, 6)), c("1", "2", "3"))
})

test_that("legacy c('tau','het_cov') renders dodged τ + stacked sub bars", {
  cs <- cross_fixture()
  p <- plot(cs, type = "wibe", stats = c("tau", "het_cov"))
  # First layer: τ bars (3 rows, dodged left of integer x)
  tau_ld <- ggplot2::layer_data(p, 1L)
  expect_equal(nrow(tau_ld), 3L)
  expect_true(all(tau_ld$x < c(1, 2, 3)))
  # Second layer: stacked sub bars (6 rows, dodged right)
  sub_ld <- ggplot2::layer_data(p, 2L)
  expect_equal(nrow(sub_ld), 6L)
  expect_true(all(sub_ld$x > c(1, 2, 3)[match(round(sub_ld$x), c(1, 2, 3))]))
})

test_that("granular stats = c('tau_b','het_b') dodges at the Between column", {
  cs <- cross_fixture()
  p <- plot(cs, type = "wibe", stats = c("tau_b", "het_b"))
  tau_ld <- ggplot2::layer_data(p, 1L)
  sub_ld <- ggplot2::layer_data(p, 2L)
  expect_equal(nrow(tau_ld), 1L)
  expect_equal(nrow(sub_ld), 1L)
  # τ bar to the left of x = 2, het bar to the right
  expect_lt(tau_ld$x, 2)
  expect_gt(sub_ld$x, 2)
})


# ---- Longitudinal smoke tests ---------------------------------------------- #

test_that("longit plot accepts every shorthand", {
  lt <- longit_fixture()
  for (s in list("het_cov", "het_cov_b", "het_cov_w",
                 "subs", "subs_b", "subs_w",
                 c("tau", "het_cov"),
                 c("tau", "subs"))) {
    p <- plot(lt, type = "wibe", stats = s)
    expect_s3_class(p, "ggplot")
  }
})

test_that("longit plot accepts rescale stats individually", {
  lt <- longit_fixture()
  for (s in c("rescale_b", "rescale_w",
              "tau", "tau_b", "tau_w",
              "het_b", "cov_b", "het_w", "cov_w")) {
    p <- plot(lt, type = "wibe", stats = s)
    expect_s3_class(p, "ggplot")
  }
})

test_that("longit shorthand 'het_cov' equals the four-component long form", {
  lt <- longit_fixture()
  p_short <- plot(lt, type = "wibe", stats = "het_cov")
  p_long  <- plot(lt, type = "wibe",
                  stats = c("het_b", "cov_b", "het_w", "cov_w"))
  ld_short <- ggplot2::layer_data(p_short, 1L)
  ld_long  <- ggplot2::layer_data(p_long,  1L)
  expect_equal(ld_short[order(ld_short$x, ld_short$y), c("x", "y")],
               ld_long[order(ld_long$x, ld_long$y),   c("x", "y")],
               ignore_attr = TRUE)
})

test_that("longit plot rejects unknown stats", {
  lt <- longit_fixture()
  expect_error(plot(lt, type = "wibe", stats = "foo"),
               "Invalid 'stats' values: foo")
})
