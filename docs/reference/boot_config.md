# Create a bootstrap configuration

Groups bootstrap settings into a single object that can be passed to
[`ineqx`](https://benrosche.github.io/ineqx/reference/ineqx.md) via the
`se` argument, or to plot methods via the `ci` argument. All model
arguments (data, formulas, etc.) are inherited from the
[`ineqx()`](https://benrosche.github.io/ineqx/reference/ineqx.md) call
automatically.

## Usage

``` r
boot_config(
  B = 100L,
  parallel = FALSE,
  ncores = NULL,
  seed = NULL,
  verbose = TRUE
)
```

## Arguments

- B:

  Integer, number of bootstrap replicates. Default 100.

- parallel:

  Logical, use parallel computation. Default FALSE.

- ncores:

  Integer, number of cores for parallel. Default: all but one.

- seed:

  Integer, random seed for reproducibility. Default NULL.

- verbose:

  Logical, print progress messages. Default TRUE.

## Value

An object of class `"ineqx_boot_config"`

## Examples

``` r
if (FALSE) { # \dontrun{
# Bootstrap SEs with 500 replicates
ineqx("income", treat = "treat", group = "group", data = mydata,
      formula_mu = ~ treat * group,
      formula_sigma = ~ treat * group,
      se = boot_config(B = 500, seed = 42))

# Bootstrap CIs in plot (explicit config or "boot" shorthand)
plot(desc_result, ci = boot_config(B = 100))
plot(desc_result, ci = "boot")  # equivalent to boot_config()
} # }
```
