test_that(".drop_constant_factor_terms drops single-level factor terms", {
  d <- data.frame(
    earnweekf = 1:6,
    mother    = rep(0:1, 3),
    SES       = factor("A"),          # single level
    year      = rep(c(1982, 2000), 3)
  )

  # Two-sided mu formula: SES and all its interactions are dropped, rest kept
  res <- .drop_constant_factor_terms(earnweekf ~ mother * SES * year, d)
  expect_equal(res$dropped, "SES")
  expect_setequal(attr(stats::terms(res$formula), "term.labels"),
                  c("mother", "year", "mother:year"))
  expect_identical(res$formula[[2L]], quote(earnweekf))  # response preserved

  # One-sided sigma formula behaves the same and stays one-sided
  res_s <- .drop_constant_factor_terms(~ mother * SES * year, d)
  expect_equal(res_s$dropped, "SES")
  expect_length(res_s$formula, 2L)  # one-sided: ~ rhs (no response)
  expect_setequal(attr(stats::terms(res_s$formula), "term.labels"),
                  c("mother", "year", "mother:year"))
})

test_that(".drop_constant_factor_terms leaves multi-level factors untouched", {
  d <- data.frame(
    y      = 1:6,
    mother = rep(0:1, 3),
    SES    = factor(rep(c("A", "B"), 3))   # two levels
  )
  res <- .drop_constant_factor_terms(y ~ mother * SES, d)
  expect_length(res$dropped, 0L)
  expect_setequal(attr(stats::terms(res$formula), "term.labels"),
                  c("mother", "SES", "mother:SES"))
})

test_that(".drop_constant_factor_terms also catches single-level character vars", {
  d <- data.frame(
    y     = 1:4,
    x     = c(1, 2, 3, 4),
    grp   = c("A", "A", "A", "A"),   # character, one level
    stringsAsFactors = FALSE
  )
  res <- .drop_constant_factor_terms(y ~ x + grp, d)
  expect_equal(res$dropped, "grp")
  expect_setequal(attr(stats::terms(res$formula), "term.labels"), "x")
})

test_that(".drop_constant_factor_terms collapses to intercept-only when needed", {
  d <- data.frame(y = 1:4, grp = factor("A"))
  res <- .drop_constant_factor_terms(y ~ grp, d)
  expect_equal(res$dropped, "grp")
  expect_length(attr(stats::terms(res$formula), "term.labels"), 0L)
  expect_equal(attr(stats::terms(res$formula), "intercept"), 1L)
})

test_that(".drop_constant_factor_terms is a no-op when no single-level vars", {
  d <- data.frame(y = 1:4, x = c(1, 2, 3, 4), z = c(5, 6, 7, 8))
  res <- .drop_constant_factor_terms(y ~ x + z, d)
  expect_length(res$dropped, 0L)
  expect_setequal(attr(stats::terms(res$formula), "term.labels"), c("x", "z"))
})
