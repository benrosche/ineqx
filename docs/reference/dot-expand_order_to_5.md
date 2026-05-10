# Expand a meta-ordering into the underlying 5-parameter sequence.

Mapping: `behavioral -> (beta, lambda)`, `compositional -> (pi)`,
`pretreatment -> (mu, sigma)`. \\\beta\\ is switched before \\\lambda\\
within the behavioral meta-level, and \\\mu\\ before \\\sigma\\ within
the pretreatment meta-level. Users seeking full ordering invariance
should use `order = "shapley"`.

## Usage

``` r
.expand_order_to_5(order)
```

## Arguments

- order:

  Length-3 character vector, a permutation of
  `c("behavioral", "compositional", "pretreatment")`.

## Value

Length-5 character vector, a permutation of
`c("beta", "lambda", "pi", "mu", "sigma")`.
