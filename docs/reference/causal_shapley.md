# Compute Shapley values for the longitudinal causal decomposition

Computes the average decomposition across all 6 possible orderings of
the three component types (behavioral, compositional, pretreatment).
Since orderings are paired (between-group and within-group use the same
structural ordering), there are exactly 6 evaluations.

## Usage

``` r
causal_shapley(params, ref = NULL)
```

## Arguments

- params:

  An `ineqx_params` object with multiple time periods

## Value

An object of class `"ineqx_causal_longit"` with `order = "shapley"` and
additional fields:

- shapley:

  data.frame with Shapley-averaged values for each component at each
  time period

- all_orderings:

  List of 6 `ineqx_causal_longit` results, one per ordering

- ranges:

  data.frame showing min/max for each component across orderings

## Details

Shapley values provide a robustness check against path dependence in the
sequential parameter-switching decomposition. When changes over time are
small relative to levels, the ordering has little practical impact. When
changes are large, Shapley values provide a useful robustness check.
