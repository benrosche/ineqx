# Plot comparison of ineqx results

Visualizes differences between scenarios. Color encodes the
decomposition component, linetype encodes the scenario.

## Usage

``` r
# S3 method for class 'ineqx_compare'
plot(x, type = NULL, style = "line", stats = NULL, ...)
```

## Arguments

- x:

  An `ineqx_compare` object from
  [`compare`](https://benrosche.github.io/ineqx/reference/compare.md).

- type:

  Character: plot type. Defaults depend on class:

  - `ineqx_causal_longit`: `"decomp"`, `"wibe"`

  - `ineqx_causal_cross`: `"wibe"`

  - `ineqx_desc`: `"wibe"`, `"decomp"`

- style:

  Character: `"line"` or `"point"` (longitudinal only).

- stats:

  Character vector of decomposition components to display (for
  `type = "decomp"` with longitudinal data). See
  [`plot.ineqx_causal_longit`](https://benrosche.github.io/ineqx/reference/plot.ineqx_causal_longit.md)
  for valid values. Default:
  `c("behavioral", "compositional", "pretreatment", "total")`.

- ...:

  Additional arguments (currently unused).

## Value

A ggplot2 object
