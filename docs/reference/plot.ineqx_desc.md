# Plot descriptive variance decomposition

Plot descriptive variance decomposition

## Usage

``` r
# S3 method for class 'ineqx_desc'
plot(x, type = "decomp", stats = NULL, ci = FALSE, style = "line", ...)
```

## Arguments

- x:

  An `ineqx_desc` object

- type:

  Character, plot type: `"wibe"` for within/between levels over time,
  `"decomp"` for change contributions over time from delta-mu,
  delta-sigma, and delta-pi relative to the reference period (requires
  ref), `"params"` for group-level means and SDs over time, `"ineq"` for
  overall inequality statistics over time, `"ineq.group"` for per-group
  inequality statistics over time.

- stats:

  For `type = "ineq"` or `"ineq.group"`: a list of inequality statistics
  to compute. Can be character strings (`"V"`, `"VL"`, `"CV2"`,
  `"Gini"`, `"Theil"`) or custom functions with signature
  `f(y, w = NULL)`. Default: `c("V")`.

- ci:

  How to compute confidence intervals. Options: `FALSE` or `"none"`
  (default): no CIs. `TRUE` or `"delta"`: delta method CIs (reads stored
  SEs from
  [`ineqx()`](https://benrosche.github.io/ineqx/reference/ineqx.md)).
  `"boot"`: bootstrap CIs with default settings.
  [`boot_config()`](https://benrosche.github.io/ineqx/reference/boot_config.md):
  bootstrap CIs with custom settings.

- style:

  Character: `"line"` (default) for connected lines with ribbon CIs, or
  `"point"` for points with error bar CIs.

- ...:

  Additional arguments (currently unused)

## Value

A ggplot2 object
