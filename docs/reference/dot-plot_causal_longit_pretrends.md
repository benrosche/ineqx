# Longitudinal pretrends diagnostic plot

For DiD models, plots the pre-period predicted Treated and Control
levels across calendar time. Under parallel trends, the two lines should
evolve with the same slope; their (constant) vertical gap reflects
pre-existing selection (\\\beta_D\\). Diverging slopes indicate a
violation of parallel trends.

## Usage

``` r
.plot_causal_longit_pretrends(x, ci = FALSE, style = "line")
```

## Details

Both panels (mu and sigma) are shown.
