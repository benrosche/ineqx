# Plot cross-sectional causal decomposition

Plot cross-sectional causal decomposition

## Usage

``` r
# S3 method for class 'ineqx_causal_cross'
plot(x, type = "wibe", ci = FALSE, stats = NULL, trim = 0.995, ...)
```

## Arguments

- x:

  An `ineqx_causal_cross` object

- type:

  Character, the plot type:

  `"wibe"`

  :   Within/between treatment effects (bar chart).

  `"wibe.group"`

  :   Group-level contributions to within/between.

  `"treat"`

  :   Predicted treatment effect distributions by group.

  `"treat.params"`

  :   Treatment effect parameters (beta, lambda) by group (bar chart).

  `"outcome"`

  :   Predicted outcome distributions (control vs treated) by group.

  `"outcome.params"`

  :   Predicted means and SDs under control vs treatment. For
      simple-difference models, a bar chart per group. For DiD models, a
      netted-out line chart anchored at the pre-period treated level:
      the post-period gap between observed and counterfactual lines
      equals the DiD ATT.

- ci:

  Whether to show confidence intervals. Accepts `FALSE` or `"none"` (no
  CIs), `TRUE` or `"delta"` (use stored SEs), `"boot"` (bootstrap with
  defaults), or
  [`boot_config()`](https://benrosche.github.io/ineqx/reference/boot_config.md)
  (bootstrap with custom settings). Default `FALSE`.

- stats:

  Character vector. For `type = "wibe"`: components to display. The bar
  plot recognizes the same canonical vocabulary as the longitudinal
  `type = "wibe"` method: `"tau"` (Total τ bar), `"tau_b"` (Between τ
  bar), `"tau_w"` (Within τ bar), `"het_b"`, `"cov_b"`, `"rescale_b"`
  (Between sub-components: heterogeneity, sorting, rescaling), and the
  within mirrors `"het_w"`, `"cov_w"`, `"rescale_w"`. Single-segment
  bars are allowed. The legend labels `"cov_*"` as “Sorting” and
  `"rescale_*"` as “Rescaling”. Shorthands: `"het_cov"` = the four
  heterogeneity / sorting sub-components (no rescaling), `"het_cov_b"` =
  `c("het_b","cov_b")`, `"het_cov_w"` = `c("het_w","cov_w")`, `"subs"` =
  all six sub-components (het + cov + rescale, between and within),
  `"subs_b"` = `c("het_b","cov_b","rescale_b")`, `"subs_w"` =
  `c("het_w","cov_w","rescale_w")`. When `tau` and any sub-component
  appear at the same column, the bars are dodged. Default `"tau"`
  preserves the classic three-bar (Total/Between/Within) layout;
  `c("tau","het_cov")` preserves the dodged classic layout.

- trim:

  Numeric between 0 and 1. For distribution plots (`"treat"`,
  `"outcome"`): quantile at which to trim tails. Default `0.995` (trims
  0.5% from each tail).

- ...:

  Additional arguments (currently unused)

## Value

A ggplot2 object
