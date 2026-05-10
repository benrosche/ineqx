# Cross-sectional DiD outcome.params: netted-out classic DiD line chart

Two lines per group, both anchored at the pre-period treated level
(`mu1_pre` / `sigma1_pre`). The "Treated" line ends at the observed
treated post-period level (`mu1`); the "Counterfactual" line ends at the
DiD-implied counterfactual untreated post-period level (`mu0`). Lines
coincide at "pre" by construction; the post-period gap equals the ATT
(beta) on the mu panel and lambda on the log-SD panel.

## Usage

``` r
.plot_causal_cross_outcome_did(x, ci = FALSE)
```
