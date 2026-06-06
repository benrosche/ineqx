# Split-step workflow: fit once, decompose many

[`ineqx()`](https://benrosche.github.io/ineqx/reference/ineqx.md) has
two entry points that produce the same numbers in different ways. The
**integrated** entry takes raw data and a formula, fits a GAMLSS model
internally, and returns the decomposition in one call. The
**split-step** entry takes a pre-built `ineqx_params` object and returns
the decomposition without re-fitting GAMLSS. When you only need one
`(ystat, ref)` view you can use either; when you need several views of
the same underlying GAMLSS fit, the split-step path saves you N GAMLSS
refits.

This vignette shows the split-step pattern for the three settings where
it matters most:

- [Cross-section](#cross-section)
- [Stepwise DiD](#did)
- [V_L (variance of log earnings)](#vl)
- [When to prefer the integrated path](#integrated)

``` r

library(ineqx)
data("cps_sample")
```

------------------------------------------------------------------------

### Cross-section: one GAMLSS fit, many decompositions

The three core ingredients are
[`fit_ineqx_model()`](https://benrosche.github.io/ineqx/reference/fit_ineqx_model.md),
`ineqx_params(model = ...)`, and `ineqx(params = ...)`. The first two
are expensive (they fit GAMLSS and extract counterfactual predictions);
the third is cheap (it just runs the decomposition formula on the
already-extracted parameters).

``` r

formula_mu    <- earnweekf ~ mother * SES * year
formula_sigma <-           ~ mother * SES * year

# Step 1 — fit GAMLSS once and cache.
model <- fit_ineqx_model(
  formula_mu    = formula_mu,
  formula_sigma = formula_sigma,
  data    = cps_sample,
  weights = "earnwtf"
)
#> GAMLSS-RS iteration 1: Global Deviance = 7950375215 
#> GAMLSS-RS iteration 2: Global Deviance = 7949095852 
#> GAMLSS-RS iteration 3: Global Deviance = 7949094587 
#> GAMLSS-RS iteration 4: Global Deviance = 7949094586 
#> GAMLSS-RS iteration 5: Global Deviance = 7949094586

# Step 2 — extract counterfactual params, including the vcov needed for
# delta-method standard errors downstream. Also cache.
params_var <- ineqx_params(
  model = model, data = cps_sample,
  treat = "mother", group = "SES", time = "year",
  ystat = "Var", vcov = TRUE
)

# Step 3 — cheap. Run as many decompositions as you want.
fit_var_1982 <- ineqx(params = params_var, ystat = "Var", ref = 1982)
fit_var_2000 <- ineqx(params = params_var, ystat = "Var", ref = 2000)
```

To get a CV² decomposition from the same model, extract a CV² params
object once (this is also cheap relative to the GAMLSS step), then
re-use it the same way:

``` r

params_cv2 <- ineqx_params(
  model = model, data = cps_sample,
  treat = "mother", group = "SES", time = "year",
  ystat = "CV2", vcov = TRUE
)
fit_cv2_1982 <- ineqx(params = params_cv2, ystat = "CV2", ref = 1982)
```

------------------------------------------------------------------------

### Stepwise DiD via the split path

For a panel dataset where each subject is observed once before and once
after treatment, pass `post = "<your pre/post indicator>"` to
[`ineqx_params()`](https://benrosche.github.io/ineqx/reference/ineqx_params.md).
The function detects DiD mode (`is_did` flag) and extracts four
counterfactual predictions per row instead of two, applying the
parallel-trends identity to recover . The downstream call to
`ineqx(params = ...)` then runs the same DiD decomposition the
integrated `ineqx(post = ..., ...)` call would produce.

If you plan to decompose `CV2`, keep the modeled outcome in levels.
Applying `CV2` to first-differenced outcomes is not recommended: it
targets relative dispersion in changes, not relative inequality in the
outcome level, and can be unstable when mean changes are near zero.

``` r

did_model <- fit_ineqx_model(
  formula_mu    = y ~ mother * post * SES * year + age + race + married,
  formula_sigma =   ~ mother * post * SES * year + age + race + married,
  data    = panel_data,
  weights = "wtfinl"
)

did_params <- ineqx_params(
  model = did_model, data = panel_data,
  treat = "mother", group = "SES", time = "year",
  post  = "post01",        # <- DiD switch
  ystat = "Var", vcov = TRUE
)

fit_did_var <- ineqx(params = did_params, ystat = "Var", ref = 1982)
fit_did_cv2 <- ineqx(params = did_params, ystat = "CV2", ref = 1982)
```

The DiD params object carries `is_did = TRUE` and the name of the `post`
column, so downstream plotting and printing know to surface pre-period
anchors.

------------------------------------------------------------------------

### V_L (variance of log earnings)

`V_L` requires a GAMLSS fit on `log(y)`, not on `y`. To make this work
in the split-step path without losing the connection to `ystat = "VL"`,
pass `transform = "log"` to
[`fit_ineqx_model()`](https://benrosche.github.io/ineqx/reference/fit_ineqx_model.md).
The function log-transforms the response in a local copy of `data` and
tags the returned model so downstream
[`ineqx_params()`](https://benrosche.github.io/ineqx/reference/ineqx_params.md)
carries the tag onto the params object, and
`ineqx(params, ystat = "VL")` knows the parameters already live on the
log scale.

``` r

log_model <- fit_ineqx_model(
  formula_mu    = earnweekf ~ mother * SES * year,
  formula_sigma =           ~ mother * SES * year,
  data      = cps_sample,
  weights   = "earnwtf",
  transform = "log"           # <- log(earnweekf) before fitting
)

log_params <- ineqx_params(
  model = log_model, data = cps_sample,
  treat = "mother", group = "SES", time = "year",
  ystat = "Var", vcov = TRUE  # ystat = "Var" because params live on log scale
)

fit_vl <- ineqx(params = log_params, ystat = "VL", ref = 1982)
```

The numbers from this call are identical to those from the integrated
`ineqx(y = "earnweekf", data = cps_sample, ystat = "VL", ...)` call. The
split-step advantage is that you only pay one log-scale GAMLSS fit even
if you later want VL decompositions at multiple `ref` values.

`transform = "log"` requires the LHS of `formula_mu` to be a simple
variable name (e.g. `earnweekf`, not `log(earnweekf)`) and that column
to be strictly positive. Construct the log column yourself before
fitting if either condition fails.

------------------------------------------------------------------------

### When to prefer the integrated path

If you only need one `(ystat, ref)` view and don’t want to cache
anything, the integrated call is shorter:

``` r

ineqx(
  y = "earnweekf", ystat = "CV2",
  treat = "mother", group = "SES", time = "year", ref = 1982,
  formula_mu    = ~ mother * SES * year,
  formula_sigma = ~ mother * SES * year,
  weights = "earnwtf",
  data    = cps_sample
)
```

The integrated path also handles `ystat = "VL"` transparently — it
log-transforms `data[[y]]` internally, runs the Var decomposition, and
relabels to VL on output. Use it when the analysis is a one-shot; reach
for the split-step pattern when you need to cache an expensive GAMLSS
fit and reuse it across `ystat`, `ref`, or other downstream choices.
