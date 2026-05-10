# Compute delta method standard errors for causal decomposition

Given an `ineqx_params` object with a variance-covariance matrix,
computes standard errors for the decomposition quantities using the
delta method: \\Var(g(\hat{\theta})) \approx \nabla g(\theta)^T
\Sigma\_\theta \nabla g(\theta)\\.

## Usage

``` r
delta_method_se(params, type = c("cross", "longit"), order = NULL, ref = NULL)
```

## Arguments

- params:

  An `ineqx_params` object with a non-NULL `vcov` field

- type:

  Character, `"cross"` for cross-sectional or `"longit"` for
  longitudinal decomposition

- order:

  Character vector (for longitudinal), the decomposition ordering

## Value

A list of standard errors for the decomposition quantities
