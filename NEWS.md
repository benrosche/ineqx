# ineqx 0.4.0

## Breaking changes

* Merged `ineq()` and `ineqx()` into a single `ineqx()` function
  - If no `treat` is specified: descriptive decomposition (was `ineq()`)
  - If `treat` is specified or `params` provided: causal decomposition
* `ineq()` removed — use `ineqx()` without `treat` instead
* `y` is now a separate argument (character, column name in `data`); `formula_mu` is now one-sided (like `formula_sigma`)
* `boot_config()` now requires `y` argument and one-sided `formula_mu`
* Default `order` changed from `c("behavioral", "compositional", "pretreatment")` to `"shapley"`
* Argument order: `y, ystat, treat, post, group, time, ref, order, formula_mu, formula_sigma, params, weights, se, data`
* Renamed modes: "Mode A" → "integrated estimation", "Mode B" → "externally estimated model"

## Improvements

* Fully parallel argument structure: descriptive and causal paths share `y`, `group`, `time`, `ref`, `ystat`, `order`, `data`
* Public API reduced to 2 main functions: `ineqx()` and `ineqx_params()`

# ineqx 0.3.0

## Breaking changes

* Streamlined to 3 user-facing functions: `ineq()`, `ineqx()`, `ineqx_params()`
* `desc_decompose()` replaced by `ineq()` with enhanced counterfactual decomposition
* `extract_params()` merged into `ineqx_params(model = ...)`
* `causal_decompose_cross()`, `causal_decompose_longit()`, `causal_shapley()`, `fit_ineqx_model()` are now internal
* `ineqx()` `se_method`/`boot` arguments replaced by unified `se` argument: `"delta"` (default), `"none"`, or `boot_config()` object
* `ineqx()` no longer returns `ineqx_result` wrapper; returns decomposition object directly

## New features

* `ineq()`: descriptive W/B decomposition with counterfactual parameter-switching decomposition (mu, sigma, pi) and Shapley averaging
* `ineqx()` integrated estimation: GAMLSS fitting via `formula_mu`/`formula_sigma`/`data` arguments
* `order = "shapley"` supported for ordering-robust decomposition
* `ineqx_params()` accepts either manual data.frame or fitted gamlss model (via `model` arg)

## Improvements

* Cleaner print/plot methods for all result classes

# ineqx 0.2.0

## Changes

* Complete package redesign with new modular architecture
* Old functions removed: `wibe()`, `calcAME()`, `plot.ineqx()`
* Replaced by: `desc_decompose()`, `extract_params()`, new S3 plot methods
* The `c.`/`i.` variable prefix convention has been removed; use plain column names
* `zeallot` dependency removed; `gamlss` moved from `Depends` to `Suggests`

## New features

* `ineqx_params()`: Standardized input format that decouples estimation from decomposition. Users can construct parameters from any estimation method.
* `causal_decompose_cross()`: Cross-sectional causal decomposition (treatment effect on variance = within + between)
* `causal_decompose_longit()`: Longitudinal causal decomposition with configurable decomposition ordering (6 possible orderings via the `order` parameter)
* `causal_shapley()`: Shapley values averaged across all 6 orderings as robustness check
* `desc_decompose()`: Clean descriptive within/between decomposition
* `extract_params()`: Extract decomposition parameters from fitted GAMLSS models
* `fit_ineqx_model()`: Convenience wrapper for fitting GAMLSS models
* Both `Var` (variance) and `CV2` (squared coefficient of variation) supported throughout
* Print, summary, and plot methods for all result classes

## Bug fixes

* CV2 decomposition corrected (was marked as incorrect in v0.1.0)

# ineqx 0.1.0

* Initial release
