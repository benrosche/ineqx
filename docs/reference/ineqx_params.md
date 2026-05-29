# Create an ineqx_params object

Constructs a standardized parameter object that serves as the interface
between model estimation and variance decomposition. This decouples the
two stages: users can construct `ineqx_params` from any estimation
method (GAMLSS, Bayesian models, simulation, or manual specification).

## Usage

``` r
ineqx_params(
  data,
  model = NULL,
  treat = NULL,
  group = NULL,
  time = NULL,
  post = NULL,
  ystat = "Var",
  vcov = NULL,
  na.action = stats::na.omit,
  verbose = TRUE
)
```

## Arguments

- data:

  For descriptive mode: a data.frame with columns `group`, `pi`, `mu`,
  `sigma`. For causal manual mode: a data.frame with columns `group`,
  `pi`, `mu0`, `sigma0`, `beta`, `lambda` (and optionally `time`). For
  model mode: the data.frame used to fit the model.

- model:

  A fitted `gamlss` object. If provided, parameters are extracted from
  the model via counterfactual predictions. Default: NULL (manual
  specification).

- treat:

  Character, name of the treatment variable in `data`. Required when
  `model` is provided. Must be coded as 0 (untreated) and 1 (treated),
  or with multiple non-zero treatment levels.

- group:

  Character, name of the grouping variable in `data`. Required when
  `model` is provided.

- time:

  Character, name of the time variable. NULL for cross-sectional. Used
  in both manual mode (if `data` contains a `time` column) and model
  mode.

- post:

  Character, name of the pre/post indicator for DiD designs. Only used
  in model mode. NULL for simple difference estimator.

- ystat:

  Character, either `"Var"` (variance) or `"CV2"` (squared coefficient
  of variation). Default: `"Var"`.

- vcov:

  For manual mode: an optional variance-covariance matrix (or named list
  of matrices for longitudinal data). For model mode: logical, whether
  to extract the vcov via numerical Jacobian (default: TRUE). Set to
  FALSE for faster computation without SEs.

- na.action:

  A function that handles NAs in `data`. Only used in model mode.
  Applied to the subset of `data` containing the variables referenced by
  the model's mu/sigma formulas plus `treat`, `group`, `time`, and
  `post`. Defaults to [`na.omit`](https://rdrr.io/r/stats/na.fail.html),
  matching the integrated
  [`ineqx()`](https://benrosche.github.io/ineqx/reference/ineqx.md)
  path. Pass `na.fail` to error on NAs, or `na.pass` to keep them (and
  let [`predict()`](https://rdrr.io/r/stats/predict.html) propagate
  them).

## Value

An object of class `"ineqx_desc_params"` (descriptive mode) or
`"ineqx_params"` (causal mode)

## Details

There are three usage modes:

- Descriptive (manual):

  Pass a data.frame with columns `group`, `pi`, `mu`, `sigma`. Returns
  an `ineqx_desc_params` object representing a counterfactual reference
  for descriptive decomposition.

- Causal (manual):

  Pass a data.frame with columns `group`, `pi`, `mu0`, `sigma0`, `beta`,
  `lambda` (and optionally `time`). Returns an `ineqx_params` object.

- Causal (from fitted gamlss):

  Pass a fitted `gamlss` object via `model`, along with the raw data and
  variable names. The function extracts treatment effects via
  counterfactual predictions.

The mode is auto-detected from the columns in `data`: if `mu` and
`sigma` are present (without `mu0`, `sigma0`), the descriptive path is
used; otherwise the causal path is used.

When `model` is provided, the function generates predictions under
counterfactual treatment assignments. For each observation, it predicts
the outcome under treat=0 (baseline) and treat=1 (treated), for both the
mean (mu) and standard deviation (sigma) equations. The average marginal
effects are computed within each group-time cell.

For the simple difference estimator, the ATT is: \$\$\beta\_{D,j} =
E\[\hat{\mu}(D=1) - \hat{\mu}(D=0) \| G=j\]\$\$

## Examples

``` r
# Descriptive counterfactual reference (equal groups, no inequality)
desc_ref <- ineqx_params(
  data = data.frame(
    group = c("workers", "managers"),
    pi = c(0.5, 0.5),
    mu = c(0, 0),
    sigma = c(0, 0)
  )
)

# Causal manual specification for a two-group example
params <- ineqx_params(
  data = data.frame(
    group = c("workers", "managers"),
    pi = c(0.5, 0.5),
    mu0 = c(500, 1000),
    sigma0 = c(200, 400),
    beta = c(60, 100),
    lambda = c(-0.1, -0.2)
  )
)

if (FALSE) { # \dontrun{
# From a fitted gamlss model
params <- ineqx_params(
  model = my_gamlss, data = mydata,
  treat = "x", group = "SES",
  ystat = "Var", vcov = TRUE
)

# Stepwise DiD: pass `post`, the indicator for the post-treatment
# observation in a paired pre/post panel. ineqx_params() then runs the
# 4-prediction DiD extraction (parallel-trends adjusted) instead of the
# 2-prediction simple-difference extraction.
did_params <- ineqx_params(
  model = my_did_gamlss, data = panel_data,
  treat = "treat", group = "edu", time = "year5",
  post = "post01", ystat = "Var", vcov = TRUE
)

# V_L (variance of log earnings) via the split-step workflow:
log_model  <- fit_ineqx_model(
  formula_mu    = y ~ treat * group * time,
  formula_sigma =   ~ treat * group * time,
  data = mydata, transform = "log"
)
log_params <- ineqx_params(model = log_model, data = mydata,
                           treat = "treat", group = "group", time = "time",
                           ystat = "Var", vcov = TRUE)
fit_vl <- ineqx(params = log_params, ystat = "VL", ref = 1980)
} # }
```
