# Bootstrap CIs for descriptive plot types

Resamples rows of raw_data (stratified by time), applies compute_fn to
each replicate, and returns quantile-based CIs.

## Usage

``` r
.bootstrap_desc_ci(raw_data, B, compute_fn, level = 0.95)
```

## Arguments

- raw_data:

  Data frame with columns y, group, time, w

- B:

  Number of bootstrap replicates

- compute_fn:

  Function taking a resampled data frame and returning a data frame with
  at least a 'value' column plus grouping columns (time, Component,
  stat, group, etc.)

- level:

  Confidence level (default 0.95)

## Value

Data frame with grouping columns + ymin, ymax
