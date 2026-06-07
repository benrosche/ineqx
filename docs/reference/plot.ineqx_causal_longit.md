# Plot longitudinal causal decomposition

Plot longitudinal causal decomposition

## Usage

``` r
# S3 method for class 'ineqx_causal_longit'
plot(
  x,
  type = "decomp",
  ci = FALSE,
  style = "line",
  stats = NULL,
  time = NULL,
  trim = 0.995,
  share = FALSE,
  ...
)
```

## Arguments

- x:

  An `ineqx_causal_longit` object (including Shapley-averaged results)

- type:

  Character, the plot type:

  `"decomp"`

  :   Decomposition of changes relative to reference. Use `stats` to
      select components (default: 3 aggregate + total).

  `"wibe"`

  :   Cross-sectional treatment effects on within/between inequality at
      each time (levels, not changes).

  `"treat"`

  :   Predicted treatment effect distributions. When `time` is
      specified, shows per-group distributions. When `time = NULL`,
      shows pi-weighted marginals across all times.

  `"treat.params"`

  :   Treatment effect parameters (beta/lambda) over time (line chart).

  `"effect.prop"`

  :   Proportional treatment effect on the mean by group over time:
      `beta_g / mu0_g` (identity fit) or `exp(beta_g) - 1` (log fit).
      Overlapping lines indicate a proportional (group-invariant)
      effect.

  `"outcome"`

  :   Predicted outcome distributions (control vs treated). When `time`
      is specified, shows per-group distributions. When `time = NULL`,
      shows pi-weighted marginals.

  `"outcome.params"`

  :   Predicted means and SDs under control vs treatment over time (line
      chart). For DiD models, the "Control" line represents the
      DiD-implied counterfactual untreated level for the treated
      post-period subpopulation, so the gap between the lines equals the
      DiD ATT.

  `"pretrends"`

  :   Pre-period predicted Treated and Control levels over time.
      Available only for DiD models. Under parallel trends the two lines
      should evolve with the same slope; their (constant) vertical gap
      reflects pre-existing selection. Diverging slopes indicate a
      violation of parallel trends.

  `"shapley"`

  :   Shapley averages with ordering ranges (only when
      `order = "shapley"`). Use `stats` to select components.

- ci:

  Whether to show confidence intervals. Accepts `FALSE` or `"none"` (no
  CIs), `TRUE` or `"delta"` (use stored SEs), `"boot"` (bootstrap with
  defaults), or
  [`boot_config()`](https://benrosche.github.io/ineqx/reference/boot_config.md)
  (bootstrap with custom settings). Default `FALSE`.

- style:

  Character: `"line"` (default) for connected lines with ribbon CIs, or
  `"point"` for points with error bar CIs.

- stats:

  Character vector of components to display. Meaning depends on `type`:

  For `"decomp"` and `"shapley"`:

  :   Any combination of aggregate components (`"behavioral"`,
      `"compositional"`, `"pretreatment"`, `"total"`) and/or individual
      deltas (`"delta_beta"`, `"delta_lambda"`, `"delta_pi"`,
      `"delta_pi_b"`, `"delta_pi_w"`, `"delta_mu"`, `"delta_sigma"`).
      Default for `"decomp"`:
      `c("behavioral", "compositional", "pretreatment", "total")`.
      Default for `"shapley"`:
      `c("behavioral", "compositional", "pretreatment")`.

  For `"wibe"`:

  :   Components to display, default `c("tau", "tau_b", "tau_w")`.
      Canonical options: `"tau"` (total), `"tau_b"`, `"tau_w"`,
      `"het_b"`, `"cov_b"`, `"rescale_b"`, `"het_w"`, `"cov_w"`,
      `"rescale_w"`. Legend labels `"cov_*"` as “Sorting” and
      `"rescale_*"` as “Rescaling” (rescaling values are zero under
      `ystat = "Var"`). Shorthands: `"het_cov"` = the four het/cov
      sub-components (no rescaling), `"het_cov_b"` =
      `c("het_b","cov_b")`, `"het_cov_w"` = `c("het_w","cov_w")`,
      `"subs"` = all six sub-components, `"subs_b"` =
      `c("het_b","cov_b","rescale_b")`, `"subs_w"` =
      `c("het_w","cov_w","rescale_w")`.

- time:

  Numeric, time point for `type = "treat"` or `type = "outcome"`. If
  `NULL` (default), pi-weighted marginal distributions are shown across
  all time points.

- trim:

  Numeric between 0 and 1. For distribution plots (`"treat"`,
  `"outcome"`): quantile at which to trim tails. Default `0.995`.

- share:

  Logical. For `type = "decomp"`, plot each component as a share of the
  total effect at that time (%) instead of in absolute units. Periods
  where the total is \\\approx 0\\ are blanked, since shares are
  undefined there. Default `FALSE`.

- ...:

  Additional arguments (currently unused)

## Value

A ggplot2 object
