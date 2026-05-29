# Aggregate bootstrap SEs from cached params

Stage 2 of the two-stage bootstrap. Takes the cached `ineqx_params`
replicates produced by
[`bootstrap_params`](https://benrosche.github.io/ineqx/reference/bootstrap_params.md),
runs the variance decomposition on each at the requested `(ystat, ref)`,
and returns an `ineqx_boot` object with standard errors and percentile
CIs.

## Usage

``` r
decompose_boot_params(
  boot_params,
  ref = NULL,
  order = c("behavioral", "compositional", "pretreatment"),
  ystat = NULL
)
```

## Arguments

- boot_params:

  An object returned by
  [`bootstrap_params`](https://benrosche.github.io/ineqx/reference/bootstrap_params.md).

- ref:

  Numeric, reference time period (longitudinal) or `NULL`
  (cross-section).

- order:

  Character vector of length 3, decomposition ordering for the
  longitudinal decomposition. Default
  `c("behavioral", "compositional", "pretreatment")`. Ignored for
  cross-sectional fits.

- ystat:

  Character, `"Var"` or `"CV2"`. Default: the `ystat` stored on
  `boot_params`.

## Value

An object of class `"ineqx_boot"` matching
[`bootstrap_se`](https://benrosche.github.io/ineqx/reference/bootstrap_se.md)'s
return shape.

## Details

Because the GAMLSS fits are already cached in the input, this call only
pays the cost of B variance-decomposition evaluations — orders of
magnitude cheaper than a fresh
[`bootstrap_se()`](https://benrosche.github.io/ineqx/reference/bootstrap_se.md)
run.

`ystat` can override the value the params were originally extracted at:
the per-replicate `params$ystat` is set to the requested value before
decomposition. This works because the bootstrapped params object's
parametric columns (`mu0`, `sigma0`, `beta`, `lambda`) are
scale-agnostic; `ystat` only selects which decomposition formula to
apply on top.

## See also

[`bootstrap_params`](https://benrosche.github.io/ineqx/reference/bootstrap_params.md),
[`bootstrap_se`](https://benrosche.github.io/ineqx/reference/bootstrap_se.md).
