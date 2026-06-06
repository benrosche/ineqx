# Bootstrap standard errors for causal variance decomposition

Computes standard errors by nonparametric bootstrap: resample
individuals, re-estimate the GAMLSS, re-extract parameters, and
re-compute the decomposition B times. Standard errors are the standard
deviation across replicates. Percentile confidence intervals are also
computed.

## Usage

``` r
bootstrap_se(
  data,
  formula_mu,
  formula_sigma,
  treat,
  group,
  time = NULL,
  post = NULL,
  ref = NULL,
  ystat = "Var",
  order = c("behavioral", "compositional", "pretreatment"),
  B = 100L,
  parallel = FALSE,
  ncores = NULL,
  seed = NULL,
  verbose = TRUE,
  cl_type = NULL,
  blend_params = NULL
)
```

## Arguments

- data:

  Data.frame, the original individual-level data

- formula_mu:

  Two-sided formula for the mean (mu) equation

- formula_sigma:

  One-sided formula for the log-SD (sigma) equation

- treat:

  Character, treatment variable name

- group:

  Character, grouping variable name

- time:

  Character, time variable name. NULL for cross-sectional.

- post:

  Character, pre/post indicator for DiD. NULL for simple diff.

- ref:

  Numeric, reference time period for longitudinal decomposition

- ystat:

  Character, `"Var"` or `"CV2"`

- order:

  Character vector of length 3, decomposition ordering

- B:

  Integer, number of bootstrap replicates. Default 100.

- parallel:

  Logical, use parallel computation. Default FALSE.

- ncores:

  Integer, number of cores. Default: all but one.

- seed:

  Integer, random seed. Default NULL.

- verbose:

  Logical, print progress. Default TRUE.

- cl_type:

  Character, parallel cluster type, `"fork"` or `"psock"`. Default
  `NULL` auto-selects (fork off Windows, psock on Windows). See
  [`boot_config`](https://benrosche.github.io/ineqx/reference/boot_config.md)
  for details.

## Value

An object of class `"ineqx_boot"` containing:

- se:

  List of SEs matching
  [`delta_method_se`](https://benrosche.github.io/ineqx/reference/delta_method_se.md)
  output structure

- replicates:

  Matrix (B_successful x K) of replicate estimates

- ci:

  List of 95% percentile confidence intervals

- B:

  Total requested replicates

- B_successful:

  Number of successful replicates

- B_failed:

  Number of failed replicates

- type:

  `"cross"` or `"longit"`

- point_estimates:

  Named vector from original (non-bootstrap) decomposition

- seed:

  The random seed used

## Details

For longitudinal data with repeated cross-sections, resampling is
performed within each time period independently, preserving the sample
size per period.
