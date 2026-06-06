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
  transform = c("identity", "log"),
  na.action = stats::na.omit,
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

- transform:

  Either `"identity"` (default) or `"log"`. When `"log"`, the response
  column on the LHS of `formula_mu` is log-transformed in a local copy
  of `data` before fitting, and the returned model is tagged with
  `attr(model, "ineqx_transform") = "log"`. Downstream
  [`ineqx_params`](https://benrosche.github.io/ineqx/reference/ineqx_params.md)
  and [`ineqx`](https://benrosche.github.io/ineqx/reference/ineqx.md)
  read this tag, which is what lets the split-step workflow compute
  \\V_L\\ (variance of log earnings) via `ystat = "VL"`. Requires the
  LHS of `formula_mu` to be a simple variable name (e.g., `y`, not
  `log(y)` or `I(y + 1)`) and that column to be strictly positive.

- na.action:

  A function that handles NAs in `data`. Applied to the subset of `data`
  containing the variables referenced by `formula_mu`, `formula_sigma`,
  and `weights` before fitting. Defaults to
  [`na.omit`](https://rdrr.io/r/stats/na.fail.html), matching
  [`ineqx_params`](https://benrosche.github.io/ineqx/reference/ineqx_params.md)
  and the integrated
  [`ineqx()`](https://benrosche.github.io/ineqx/reference/ineqx.md) path
  so the split-step workflow keeps both steps on identical rows. Pass
  `na.fail` to error on NAs.

- ...:

  Additional arguments passed to
  [`gamlss::gamlss()`](https://rdrr.io/pkg/gamlss/man/gamlss.html)

## Value

A fitted `gamlss` object. When `transform = "log"`, the object carries
`attr(., "ineqx_transform") = "log"` so downstream functions know the
fit lives on the log scale.

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

## Split-step workflow

`fit_ineqx_model()` is the entry point to the split-step workflow, where
the GAMLSS fit is cached once and many decompositions (varying `ystat`
and `ref`) are derived cheaply on top:


    model  <- fit_ineqx_model(y ~ treat * group * time + controls,
                                   ~ treat * group * time + controls,
                               data = d)
    params <- ineqx_params(model = model, data = d,
                           treat = "treat", group = "group", time = "time",
                           ystat = "Var", vcov = TRUE)
    fit_var <- ineqx(params = params, ystat = "Var", ref = 1980)
    fit_cv2 <- ineqx(params = params, ystat = "CV2", ref = 1980)

For the stepwise DiD, pass `post = "post01"` (or whatever your pre/post
column is) to
[`ineqx_params()`](https://benrosche.github.io/ineqx/reference/ineqx_params.md).
For \\V_L\\, fit with `transform = "log"` and pass `ystat = "VL"` to
[`ineqx()`](https://benrosche.github.io/ineqx/reference/ineqx.md).

## Single-level grouping variables

A factor or character predictor with only one observed level cannot
enter a model formula
([`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) errors
with "contrasts can be applied only to factors with 2 or more levels").
Any such term is automatically dropped before fitting, with a message,
since the between-group component it would identify is zero by
construction. The rest of the model is estimated normally and the
decomposition returns `VarB = 0` (within equals total).

## See also

[`ineqx_params`](https://benrosche.github.io/ineqx/reference/ineqx_params.md)
for extracting decomposition parameters from the fitted model;
[`ineqx`](https://benrosche.github.io/ineqx/reference/ineqx.md) for
running the decomposition.
