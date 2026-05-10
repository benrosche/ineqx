# Extract block-diagonal coefficient vcov from a gamlss model

Bypasses `gamlss::vcov.gamlss()` / `gen.likelihood()` which use
[`get()`](https://rdrr.io/r/base/get.html) to resolve `data` from the
model's call — this fails when the model was created inside
[`fit_ineqx_model()`](https://benrosche.github.io/ineqx/reference/fit_ineqx_model.md)
because the call references local variable names that no longer exist.

## Usage

``` r
.extract_raw_vcov(model)
```

## Arguments

- model:

  A fitted gamlss object

## Value

A (p_mu + p_sigma) x (p_mu + p_sigma) block-diagonal vcov matrix

## Details

Instead, we extract the per-equation vcov directly from the QR
decompositions stored on the model (`model$mu.qr`, `model$sigma.qr`) and
combine them into a block-diagonal matrix.
