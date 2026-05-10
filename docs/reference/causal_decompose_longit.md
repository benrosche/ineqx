# Longitudinal causal variance decomposition

Decomposes the change in treatment effect on inequality from a reference
period to each subsequent period into behavioral, compositional, and
pre-treatment components, using a unified sequential parameter switching
scheme that tracks each switch's contribution to both \\\tau_B\\ and
\\\tau_W\\ simultaneously.

## Usage

``` r
causal_decompose_longit(
  params,
  order = c("behavioral", "compositional", "pretreatment"),
  ref = NULL
)
```

## Arguments

- params:

  An `ineqx_params` object with multiple time periods

- order:

  Character vector of length 3, a permutation of
  `c("behavioral", "compositional", "pretreatment")` specifying the
  order in which parameter groups are switched from baseline to time-t
  values. The mapping to the underlying 5-parameter sequence is:
  `behavioral -> (beta, lambda)`, `compositional -> (pi)`,
  `pretreatment -> (mu, sigma)`, with \\\beta\\ switched before
  \\\lambda\\ and \\\mu\\ switched before \\\sigma\\ within their
  meta-levels. Use `order = "shapley"` to average over all 6
  meta-orderings. Default:
  `c("behavioral", "compositional", "pretreatment")`.

## Value

An object of class `"ineqx_causal_longit"` containing:

- results:

  List keyed by time period, each containing the 6 components plus 3
  combined components and the total

- order:

  The ordering used

- ystat:

  The inequality measure

- params:

  The input ineqx_params object

## Details

For each parameter \\p \in \\\beta, \lambda, \pi, \mu, \sigma\\\\, the
result reports two split components:

- Delta\_\<p\>\_B:

  Contribution of switching p (from t0 to t) to the change in \\\tau_B\\

- Delta\_\<p\>\_W:

  Contribution of switching p (from t0 to t) to the change in \\\tau_W\\

For the variance (`ystat = "Var"`), \\\tau_B\\ depends only on \\(\pi,
\mu, \beta)\\ and \\\tau_W\\ only on \\(\pi, \sigma, \lambda)\\, so the
off-diagonal split parts
(`Delta_beta_W, Delta_lambda_B, Delta_mu_W, Delta_sigma_B`) are exactly
zero. For CV\\^2\\, both \\\tau_B\\ and \\\tau_W\\ share the grand-mean
denominator, so all five parameters can contribute to both sides; the
split components capture this coupling exactly.

By construction, the split components sum to the cross-sectional change:
\\\sum_p \text{Delta\\\<p\>\\B} = \tau_B(t) - \tau_B(t_0)\\ and
similarly for the W side, so `Delta_total` equals \\(\tau_B(t) +
\tau_W(t)) - (\tau_B(t_0) + \tau_W(t_0))\\.

Aggregate parameter components
(`Delta_beta, Delta_lambda, Delta_mu, Delta_sigma`) are reported as the
sum of their B and W parts; for V they equal the previous (single-side)
values exactly. Compositional effects remain split as `Delta_pi_B` and
`Delta_pi_W`.
