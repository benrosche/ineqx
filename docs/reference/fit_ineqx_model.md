# Fit a GAMLSS model for variance decomposition

A convenience wrapper around
[`gamlss::gamlss()`](https://rdrr.io/pkg/gamlss/man/gamlss.html) for
users who want to fit a GAMLSS model without worrying about the syntax.
The model simultaneously estimates the conditional mean and conditional
(log) standard deviation, which are the inputs to the causal variance
decomposition.

## Usage

``` r
fit_ineqx_model(
  formula_mu,
  formula_sigma,
  data,
  weights = NULL,
  family = NULL,
  ...
)
```

## Arguments

- formula_mu:

  Formula for the mean equation (e.g., `y ~ treat * group + controls`)

- formula_sigma:

  Formula for the log-SD equation (e.g., `~ treat * group + controls`).
  Note: one-sided formula.

- data:

  A data.frame

- weights:

  Optional character string naming the weight variable in `data`, or a
  numeric vector of weights

- family:

  A `gamlss.family` object specifying the distribution. Default:
  [`gamlss.dist::NO()`](https://rdrr.io/pkg/gamlss.dist/man/NO.html)
  (normal distribution).

- ...:

  Additional arguments passed to
  [`gamlss::gamlss()`](https://rdrr.io/pkg/gamlss/man/gamlss.html)

## Value

A fitted `gamlss` object

## Details

Users who want more control over the model specification should fit
their own `gamlss` model directly and use
[`ineqx_params`](https://benrosche.github.io/ineqx/reference/ineqx_params.md)
to extract the decomposition parameters.

The GAMLSS framework models both the conditional mean and conditional
variance of the outcome distribution. The mean equation captures how
treatment affects average outcomes within each group. The sigma equation
captures how treatment affects the dispersion of outcomes within each
group.

For the simple difference estimator, the model is: \$\$\mu\_{ij} =
\alpha_j + \beta\_{D,j} D_i + Z_i \gamma\$\$ \$\$\log(\sigma\_{ij}) =
\alpha_j^{(\sigma)} + \lambda\_{D,j} D_i + Z_i \gamma^{(\sigma)}\$\$

For the DiD estimator, add pre/post interactions: \$\$\mu\_{ij} =
\alpha_j + \beta\_{D,j} D_i + \beta_P P_i + \beta\_{DP,j} D_i P_i + Z_i
\gamma\$\$

## See also

[`ineqx_params`](https://benrosche.github.io/ineqx/reference/ineqx_params.md)
for extracting decomposition parameters from the fitted model
