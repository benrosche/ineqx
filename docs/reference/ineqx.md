# Variance decomposition

A unified function for both descriptive and causal variance
decomposition. If `treat` is not specified (and no `params` are
provided), a descriptive within/between decomposition is performed. If
`treat` is specified or `params` are provided, a causal decomposition of
the treatment effect on inequality is performed.

## Usage

``` r
ineqx(
  y = NULL,
  ystat = "Var",
  treat = NULL,
  post = NULL,
  group = NULL,
  time = NULL,
  ref = NULL,
  order = "shapley",
  formula_mu = NULL,
  formula_sigma = NULL,
  estimand = c("marginal", "residual"),
  params = NULL,
  weights = NULL,
  se = "delta",
  data = NULL,
  ...
)
```

## Arguments

- y:

  Character, name of the outcome variable in `data`.

- ystat:

  Character, one of `"Var"` (default), `"CV2"`, or `"VL"` (variance of
  log). `"VL"` is supported for descriptive decomposition and for
  integrated causal estimation (when `treat`, `formula_mu`, and
  `formula_sigma` are supplied with raw `y`/`data`). It is implemented
  by log-transforming `y` and running the standard `"Var"`
  decomposition; output is labelled as `"VL"`. `"VL"` is not supported
  when `params` is supplied; in that case, fit your model on `log(y)`
  yourself and pass `ystat = "Var"`.

- treat:

  Character, name of the treatment variable in `data`. Coded 0/1. If
  NULL, a descriptive decomposition is performed.

- post:

  Character, pre/post indicator for DiD designs. NULL for simple
  difference estimator. Only used in integrated estimation mode. For
  `ystat = "CV2"`, keep `y` on the outcome's level scale and use `post`
  for the DiD contrast; applying `"CV2"` to first-differenced outcomes
  is not recommended because it targets relative dispersion in changes
  and can be unstable when mean changes are near zero.

- group:

  Character, name of the grouping variable in `data`.

- time:

  Character, name of the time variable in `data`. If NULL, a single
  cross-section is assumed.

- ref:

  Numeric, reference time period. Required for longitudinal
  decomposition (both descriptive and causal).

- order:

  Decomposition ordering. Either:

  - `"shapley"` (default): averages across all 6 possible orderings.

  - For descriptive: a permutation of `c("mu", "sigma", "pi")`.

  - For causal: a permutation of
    `c("behavioral", "compositional", "pretreatment")`.

  Only relevant for longitudinal data with a reference period.

- formula_mu:

  One-sided formula for the mean equation (integrated estimation mode).
  E.g., `~ treat * group + controls`. The outcome `y` is prepended
  internally.

- formula_sigma:

  One-sided formula for the log-SD equation (integrated estimation
  mode). E.g., `~ treat * group + controls`.

- estimand:

  Character, either `"marginal"` (default) or `"residual"`. Selects
  whether within-group dispersion is the marginal counterfactual
  variance (law of total variance over the covariate distribution;
  controls contribute to within-group inequality) or the
  residual/conditional scale (paper Appendix B.7). See
  [`ineqx_params`](https://benrosche.github.io/ineqx/reference/ineqx_params.md).
  Only used in the integrated- and blending-estimation paths.

- params:

  A parameter object created by
  [`ineqx_params`](https://benrosche.github.io/ineqx/reference/ineqx_params.md).
  Either an `ineqx_desc_params` (descriptive counterfactual reference)
  or an `ineqx_params` (causal decomposition). When an `ineqx_params` is
  provided together with `treat`, `formula_mu`, and `formula_sigma`,
  blending mode is used: the params define a counterfactual baseline,
  and a GAMLSS model is fitted to estimate parameters for observed
  periods. See Details.

- weights:

  Character, name of the weight variable in `data`. If NULL, equal
  weights are used.

- se:

  Standard error method. One of:

  - `"delta"` or `TRUE` (default): delta method SEs. For causal,
    requires vcov in params. For descriptive, uses sampling-based
    covariance.

  - `"none"` or `FALSE`: skip SE computation

  - `"boot"`: bootstrap SEs with default settings (equivalent to
    `se = boot_config()`)

  - A
    [`boot_config`](https://benrosche.github.io/ineqx/reference/boot_config.md)
    object: bootstrap SEs with custom settings (causal only)

- data:

  A data.frame containing the variables. For manual descriptive mode
  (when `y = NULL`), must contain columns `group`, `pi`, `mu`, `sigma`,
  and optionally a time column.

- ...:

  Additional arguments passed to
  [`gamlss::gamlss()`](https://rdrr.io/pkg/gamlss/man/gamlss.html) for
  the integrated- and blending-estimation paths. Useful for bumping the
  iteration limit on saturated models (`n.cyc = 100`) or selecting a
  non-default family.

## Value

For descriptive decomposition: an `ineqx_desc` object. For
cross-sectional causal: an `ineqx_causal_cross` object. For longitudinal
causal (both specific orderings and Shapley): an `ineqx_causal_longit`
object. Shapley results include additional fields `$shapley`,
`$all_orderings`, and `$ranges`.

## Details

There are five usage modes:

- Descriptive (raw data):

  Provide `y` but omit `treat` and `params`. Decomposes total inequality
  into within- and between-group components. For longitudinal data with
  a reference period, further decomposes the change over time into
  contributions from changing means, dispersions, and group composition.

- Descriptive (counterfactual reference):

  Provide `y`, `data`, and an `ineqx_desc_params` object via `params`.
  The params define a counterfactual baseline; observed group-level
  statistics are computed from the data. The decomposition shows how
  each parameter contributes to the change from counterfactual to
  observed inequality.

- Causal (integrated estimation):

  Provide `treat`, `formula_mu`, `formula_sigma`, and `data`. The
  function fits a GAMLSS model, extracts parameters, and performs the
  causal decomposition.

- Causal (counterfactual reference):

  Provide an `ineqx_params` object via `params` together with `treat`,
  `formula_mu`, and `formula_sigma`. The params define a counterfactual
  baseline (e.g., zero treatment effects), and a GAMLSS model is fitted
  to estimate parameters for observed periods. The two are blended for a
  longitudinal decomposition showing how treatment effects changed from
  the baseline. Delta method SEs are not available; use
  [`boot_config()`](https://benrosche.github.io/ineqx/reference/boot_config.md)
  for bootstrap SEs.

- Causal (externally estimated model):

  Provide an `ineqx_params` object via `params` (without
  `treat`/formulas). Use
  [`ineqx_params`](https://benrosche.github.io/ineqx/reference/ineqx_params.md)
  to create one manually or from a fitted gamlss model.

## Examples

``` r
data(incdat)

# Descriptive decomposition from raw data
ineqx("inc", group = "group", time = "year", data = incdat, ref = 1)
#> Computing descriptive decomposition...
#> Finished.
#> Descriptive variance decomposition
#> Inequality measure: Var 
#> Reference period: 1 
#> Ordering: shapley 
#> 
#> Totals by time:
#>  time   VarW     VarB   VarT
#>     1  40.28  1.23480  41.52
#>     2  59.28  0.03878  59.32
#>     3  90.50  1.83958  92.34
#>     4 136.19  6.49804 142.69
#>     5 192.58 12.92182 205.50
#> 
#> Decomposition of changes in Var relative to ref = 1:
#> 
#>   time 1: (reference)
#> 
#>   time 2:
#>     Between-group (delta_mu):                   -1.1960
#>     Within-group (delta_sigma):                 18.9966
#>     Compositional (delta_pi):                    0.0000
#>       Between-group:                             0.0000
#>       Within-group:                              0.0000
#>     Total:                                      17.8005
#> 
#>   time 3:
#>     Between-group (delta_mu):                    0.6048
#>     Within-group (delta_sigma):                 50.2203
#>     Compositional (delta_pi):                    0.0000
#>       Between-group:                             0.0000
#>       Within-group:                              0.0000
#>     Total:                                      50.8251
#> 
#>   time 4:
#>     Between-group (delta_mu):                    5.2632
#>     Within-group (delta_sigma):                 95.9118
#>     Compositional (delta_pi):                    0.0000
#>       Between-group:                             0.0000
#>       Within-group:                              0.0000
#>     Total:                                     101.1750
#> 
#>   time 5:
#>     Between-group (delta_mu):                   11.6870
#>     Within-group (delta_sigma):                152.2994
#>     Compositional (delta_pi):                    0.0000
#>       Between-group:                             0.0000
#>       Within-group:                              0.0000
#>     Total:                                     163.9864

# Descriptive with counterfactual reference
ref <- ineqx_params(data = data.frame(
  group = 1:3, pi = 1/3, mu = 0, sigma = 0
))
ineqx("inc", group = "group", time = "year",
      params = ref, data = incdat)
#> Computing descriptive decomposition...
#> Finished.
#> Descriptive variance decomposition
#> Inequality measure: Var 
#> Reference period: 0 
#> Ordering: shapley 
#> 
#> Totals by time:
#>  time   VarW     VarB   VarT
#>     0   0.00  0.00000   0.00
#>     1  40.28  1.23480  41.52
#>     2  59.28  0.03878  59.32
#>     3  90.50  1.83958  92.34
#>     4 136.19  6.49804 142.69
#>     5 192.58 12.92182 205.50
#> 
#> Decomposition of changes in Var relative to counterfactual reference:
#> 
#>   time 0: (reference)
#> 
#>   time 1:
#>     Between-group (delta_mu):                    1.2348
#>     Within-group (delta_sigma):                 40.2821
#>     Compositional (delta_pi):                    0.0000
#>       Between-group:                             0.0000
#>       Within-group:                              0.0000
#>     Total:                                      41.5169
#> 
#>   time 2:
#>     Between-group (delta_mu):                    0.0388
#>     Within-group (delta_sigma):                 59.2786
#>     Compositional (delta_pi):                    0.0000
#>       Between-group:                             0.0000
#>       Within-group:                              0.0000
#>     Total:                                      59.3174
#> 
#>   time 3:
#>     Between-group (delta_mu):                    1.8396
#>     Within-group (delta_sigma):                 90.5024
#>     Compositional (delta_pi):                    0.0000
#>       Between-group:                             0.0000
#>       Within-group:                              0.0000
#>     Total:                                      92.3419
#> 
#>   time 4:
#>     Between-group (delta_mu):                    6.4980
#>     Within-group (delta_sigma):                136.1939
#>     Compositional (delta_pi):                    0.0000
#>       Between-group:                             0.0000
#>       Within-group:                              0.0000
#>     Total:                                     142.6919
#> 
#>   time 5:
#>     Between-group (delta_mu):                   12.9218
#>     Within-group (delta_sigma):                192.5815
#>     Compositional (delta_pi):                    0.0000
#>       Between-group:                             0.0000
#>       Within-group:                              0.0000
#>     Total:                                     205.5033

# Causal: from manual params (externally estimated model)
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
ineqx("inc", group = "group", data = incdat, params = params, se = "none")
#> Computing decomposition...
#> Finished.
#> Cross-sectional causal variance decomposition
#> Inequality measure: Var 
#> Estimand: marginal 
#> 
#> Treatment effect on outcome variance:
#>     Var[Y | T = 0]:                      162500.0000
#>     Var[Y | T = 1]:                      142900.2187
#>     Total effect (tau_T):                -19599.7813
#>       Between-group (tau_B):             10400.0000
#>       Within-group (tau_W):              -29999.7813
#> 
#> Between-group sub-components:
#>     Var_pi(beta):                          400.0000
#>     2*Cov_pi(mu0, beta):                 10000.0000
#> 
#> Within-group sub-components:
#>     mean(sigma0^2) * mean(f):            -25547.4600
#>     Cov_pi(sigma0^2, f):                 -4452.3212

if (FALSE) { # \dontrun{
# Causal: integrated estimation
ineqx("inc", treat = "x", group = "group", data = incdat,
      formula_mu = ~ x * factor(group),
      formula_sigma = ~ x * factor(group),
      se = "delta")

# Causal: counterfactual reference (blending)
cf_ref <- ineqx_params(data = data.frame(
  group = 1:3, pi = 1/3,
  mu0 = 500, sigma0 = 100, beta = 0, lambda = 0
))
ineqx("inc", treat = "x", group = "group",
      params = cf_ref,
      formula_mu = ~ x * factor(group),
      formula_sigma = ~ x * factor(group),
      se = "none", data = incdat)
} # }
```
