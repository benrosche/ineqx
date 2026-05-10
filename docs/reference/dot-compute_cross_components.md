# Compute exact sub-components of the cross-sectional treatment effect

For Var: delta_B = Var_pi(beta) + 2\*Cov_pi(mu0, beta) (exact identity)
delta_W = mean_pi(sigma0^2)\*mean_pi(f) + Cov_pi(sigma0^2, f) (exact)

## Usage

``` r
.compute_cross_components(pi, mu0, sigma0, beta, lambda, ystat)
```

## Arguments

- pi:

  Group shares

- mu0:

  Baseline group means

- sigma0:

  Baseline group SDs

- beta:

  Treatment effects on mean

- lambda:

  Treatment effects on log-SD

- ystat:

  "Var" or "CV2"

## Value

Named list of sub-components
