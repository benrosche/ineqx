# Package index

## Main entry points

Unified function for descriptive and causal variance decomposition.

- [`ineqx()`](https://benrosche.github.io/ineqx/reference/ineqx.md) :
  Variance decomposition
- [`ineqx_params()`](https://benrosche.github.io/ineqx/reference/ineqx_params.md)
  : Create an ineqx_params object

## Standard errors

Bootstrap configuration; delta-method SEs are computed automatically
inside ineqx().

- [`boot_config()`](https://benrosche.github.io/ineqx/reference/boot_config.md)
  : Create a bootstrap configuration
- [`bootstrap_params()`](https://benrosche.github.io/ineqx/reference/bootstrap_params.md)
  : Bootstrap resampled parameter objects (without decomposition)
- [`decompose_boot_params()`](https://benrosche.github.io/ineqx/reference/decompose_boot_params.md)
  : Aggregate bootstrap SEs from cached params

## Comparing scenarios

- [`compare()`](https://benrosche.github.io/ineqx/reference/compare.md)
  : Compare multiple ineqx results

## Plotting & theme

- [`plot(`*`<ineqx_causal_cross>`*`)`](https://benrosche.github.io/ineqx/reference/plot.ineqx_causal_cross.md)
  : Plot cross-sectional causal decomposition
- [`plot(`*`<ineqx_causal_longit>`*`)`](https://benrosche.github.io/ineqx/reference/plot.ineqx_causal_longit.md)
  : Plot longitudinal causal decomposition
- [`plot(`*`<ineqx_compare>`*`)`](https://benrosche.github.io/ineqx/reference/plot.ineqx_compare.md)
  : Plot comparison of ineqx results
- [`plot(`*`<ineqx_desc>`*`)`](https://benrosche.github.io/ineqx/reference/plot.ineqx_desc.md)
  : Plot descriptive variance decomposition
- [`theme_ineqx`](https://benrosche.github.io/ineqx/reference/theme_ineqx.md)
  : ineqx ggplot2 theme

## Print & summary

## Datasets

- [`cps_sample`](https://benrosche.github.io/ineqx/reference/cps_sample.md)
  : CPS sample data
- [`incdat`](https://benrosche.github.io/ineqx/reference/incdat.md) :
  Income data

## Shiny demo

- [`runShinyExample()`](https://benrosche.github.io/ineqx/reference/runShinyExample.md)
  : Run Shiny example app
