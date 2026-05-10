# Unified sequential decomposition switching all 5 parameters and tracking contributions to both \\\tau_B\\ and \\\tau_W\\ at each step.

Walks through the parameter sequence in `order_5`, switching each
parameter from its t0 value to its t value and recording the resulting
change in both the cross-sectional between-group treatment effect
\\\tau_B\\ and the cross-sectional within-group treatment effect
\\\tau_W\\. Telescopes exactly: summing all (B, W) contributions
recovers \\(\tau_B(t) - \tau_B(t_0))\\ and \\(\tau_W(t) -
\tau_W(t_0))\\.

## Usage

``` r
.decompose_sequential(d0, dt, order_5, ystat)
```

## Arguments

- d0:

  Reference period parameters (data.frame with pi, mu0, sigma0, beta,
  lambda)

- dt:

  Target period parameters (same columns)

- order_5:

  5-permutation of `c("beta", "lambda", "pi", "mu", "sigma")`

- ystat:

  `"Var"` or `"CV2"`

## Value

Named list. `parts[[p]]$B` and `parts[[p]]$W` for each parameter
`p in c("beta","lambda","pi","mu","sigma")`, plus `parts$totals` with
`delta_B`, `delta_W`, `delta_T`.
