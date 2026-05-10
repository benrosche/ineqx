# Compute treatment effect on both W and B components

Compute treatment effect on both W and B components

## Usage

``` r
compute_delta_WB(pi, mu0, sigma0, beta, lambda, ystat)
```

## Arguments

- pi:

  Group shares among treated

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

List with delta_B, delta_W, delta_total
