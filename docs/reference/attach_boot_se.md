# Attach bootstrap SEs to a causal decomposition result

Replaces the standard errors in a decomposition result object with
bootstrap SEs. This is useful when you've computed the decomposition and
bootstrap SEs separately.

## Usage

``` r
attach_boot_se(result, boot)
```

## Arguments

- result:

  An `ineqx_causal_cross` or `ineqx_causal_longit` object

- boot:

  An `ineqx_boot` object

## Value

The `result` object with `$se` replaced by bootstrap SEs and additional
fields `$se_method` and `$boot`.
