# Within-group variance (pi-based)

Within-group variance (pi-based)

## Usage

``` r
VarW_pi(pi, sigma)
```

## Arguments

- pi:

  Numeric vector of group shares (must sum to 1)

- sigma:

  Numeric vector of group standard deviations

## Value

Scalar within-group variance: sum(pi \* sigma^2)
