# Cross-sectional causal variance decomposition

Decomposes the treatment effect on total variance (or CV^2) into
between-group and within-group components at a single timepoint.

## Usage

``` r
causal_decompose_cross(params, ref = NULL)
```

## Arguments

- params:

  An `ineqx_params` object. For cross-sectional decomposition, this
  should contain a single time period (or no time column).

## Value

An object of class `"ineqx_causal_cross"` containing:

- tau_B:

  Scalar, between-group treatment effect on inequality

- tau_W:

  Scalar, within-group treatment effect on inequality

- tau_total:

  Scalar, total treatment effect (tau_B + tau_W)

- components:

  List of sub-components for interpretive analysis

- by_group:

  data.frame of group-level contributions

- params:

  The input ineqx_params object

## Details

The between-group treatment effect is: \$\$\delta_B^D =
\text{Var}\_\pi(\beta_D) + 2\\\text{Cov}\_\pi(\mu(0), \beta_D)\$\$

The within-group treatment effect is: \$\$\delta_W^D = \sum_j \pi_j
\sigma_j^2(0) \[\exp(2\lambda\_{D,j}) - 1\]\$\$
