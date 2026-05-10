# Compute treatment effect on within-group inequality

Compute treatment effect on within-group inequality

## Usage

``` r
compute_delta_W(pi, mu0, sigma0, lambda, ystat)
```

## Arguments

- pi:

  Group shares among treated

- mu0:

  Baseline group means (needed for CV2)

- sigma0:

  Baseline group SDs sigma_j(0)

- lambda:

  Treatment effects on log-SD

- ystat:

  "Var" or "CV2"

## Value

Scalar delta_W^D
