# Compute delta method SEs for descriptive decomposition

Given group-level statistics (pi, mu, sigma, n), constructs the sampling
covariance matrix and applies the delta method to compute SEs for VarW,
VarB, VarT, and (if applicable) counterfactual deltas.

## Usage

``` r
delta_method_desc_se(
  wibe,
  totals,
  deltas,
  ystat,
  ref = NULL,
  order = "shapley"
)
```

## Arguments

- wibe:

  Data.frame with columns: time, group, n, pi, mu, sigma, sigma2

- totals:

  Data.frame with columns: time, VarW, VarB, VarT, etc.

- deltas:

  Data.frame with delta columns, or NULL

- ystat:

  Character, "Var" or "CV2"

- ref:

  Reference time period (needed for delta SEs), or NULL

- order:

  Decomposition order: "shapley" or permutation of c("mu","sigma","pi")

## Value

A list with:

- totals:

  Named list by time, each containing se_VarW, se_VarB, se_VarT (or
  se_CV2W, se_CV2B, se_CV2T)

- deltas:

  Named list by time, each containing se_delta_mu, se_delta_sigma,
  se_delta_pi, se_delta_T (if deltas provided)
