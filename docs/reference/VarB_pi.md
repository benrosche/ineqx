# Between-group variance (pi-based)

Between-group variance (pi-based)

## Usage

``` r
VarB_pi(pi, mu)
```

## Arguments

- pi:

  Numeric vector of group shares (must sum to 1)

- mu:

  Numeric vector of group means

## Value

Scalar between-group variance: sum(pi \* (mu - grand_mean)^2)
