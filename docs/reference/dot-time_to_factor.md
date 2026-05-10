# Convert time column to ordered factor for discrete x-axis

Ensures time points are displayed as evenly-spaced categories, which
avoids large gaps when counterfactual time values (e.g. 0) are mixed
with calendar years (e.g. 1980, 1985).

## Usage

``` r
.time_to_factor(time_vals)
```
