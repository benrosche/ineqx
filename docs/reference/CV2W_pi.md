# Within-group CV^2 (pi-based)

Within-group CV^2 (pi-based)

## Usage

``` r
CV2W_pi(pi, mu, sigma)
```

## Arguments

- pi:

  Numeric vector of group shares

- mu:

  Numeric vector of group means

- sigma:

  Numeric vector of group standard deviations

## Value

Scalar within-group CV^2: sum(pi \* sigma^2) / (sum(pi \* mu))^2
