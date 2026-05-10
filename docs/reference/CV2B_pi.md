# Between-group CV^2 (pi-based)

Between-group CV^2 (pi-based)

## Usage

``` r
CV2B_pi(pi, mu)
```

## Arguments

- pi:

  Numeric vector of group shares

- mu:

  Numeric vector of group means

## Value

Scalar between-group CV^2: sum(pi \* (mu - grand_mean)^2) / grand_mean^2
