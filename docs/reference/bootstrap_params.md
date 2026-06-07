# Bootstrap resampled parameter objects (without decomposition)

Stage 1 of a two-stage bootstrap. Runs B replicates of
`(resample -> fit GAMLSS -> ineqx_params)` and returns the resampled
`ineqx_params` objects together with the original (unresampled) fit and
the resampling indices used. The standard errors come later, from
[`decompose_boot_params`](https://benrosche.github.io/ineqx/reference/decompose_boot_params.md).

## Usage

``` r
bootstrap_params(
  data,
  formula_mu,
  formula_sigma,
  treat,
  group,
  time = NULL,
  post = NULL,
  ystat = "Var",
  estimand = "marginal",
  weights = NULL,
  B = 100L,
  parallel = FALSE,
  ncores = NULL,
  seed = NULL,
  verbose = TRUE,
  cl_type = NULL
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

- ystat:

  Character, `"Var"` or `"CV2"`

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

An object of class `"ineqx_boot_params"` with elements:

- params_list:

  List of B (possibly \< B if some replicates failed) `ineqx_params`
  objects from resampled GAMLSS fits.

- boot_indices:

  List of B integer vectors of row indices used for resampling. Same
  RNG-derived indices as
  [`bootstrap_se()`](https://benrosche.github.io/ineqx/reference/bootstrap_se.md)
  would have produced with the same seed.

- B:

  Requested replicates.

- B_successful:

  Replicates that produced a usable params object.

- B_failed:

  Replicates that errored out.

- type:

  `"cross"` or `"longit"`.

- seed:

  The random seed used.

- ystat:

  The `ystat` the params were extracted at; informational only —
  [`decompose_boot_params()`](https://benrosche.github.io/ineqx/reference/decompose_boot_params.md)
  can override per call.

## Details

Use this when you want to compute bootstrap SEs for several
decomposition views of the same GAMLSS fit (e.g. different `ystat`,
different `ref`) without paying the GAMLSS fit cost B times per view.
The expensive work (GAMLSS fit + counterfactual prediction) runs once
per replicate here; the cheap part (variance decomposition) runs as many
times as you want via
[`decompose_boot_params`](https://benrosche.github.io/ineqx/reference/decompose_boot_params.md).

For the simple "one decomposition, one bootstrap" case, prefer the
integrated
[`bootstrap_se`](https://benrosche.github.io/ineqx/reference/bootstrap_se.md)
which wraps this two-stage flow.

## See also

[`decompose_boot_params`](https://benrosche.github.io/ineqx/reference/decompose_boot_params.md)
to consume this object;
[`bootstrap_se`](https://benrosche.github.io/ineqx/reference/bootstrap_se.md)
for the single-stage shortcut.

## Examples

``` r
if (FALSE) { # \dontrun{
# Two-stage workflow: one boot, many decompositions
bp <- bootstrap_params(
  data = mydata,
  formula_mu    = y ~ treat * group * time + age,
  formula_sigma =   ~ treat * group * time + age,
  treat = "treat", group = "group", time = "time",
  ystat = "Var", B = 500, seed = 42, parallel = TRUE
)

se_var_1980 <- decompose_boot_params(bp, ref = 1980, ystat = "Var")
se_cv2_1980 <- decompose_boot_params(bp, ref = 1980, ystat = "CV2")
se_var_2000 <- decompose_boot_params(bp, ref = 2000, ystat = "Var")
} # }
```
