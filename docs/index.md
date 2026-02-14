# ineqx: Descriptive and causal variance decompositions

![](reference/figures/ineqx-hexagon.jpg)

The ineqx package allows to analyze how inequality in an outcome (e.g.,
income) splits into inequality within and between groups (e.g., gender).
It is possible to decompose inequality at a single point in time and to
decompose changes in inequality over time. In addition to this
descriptive decomposition, the ineqx packages allows to analyze how
treatment effects (i.e., binary predictors) impact within- and
between-group inequality and how this effect changes over time.

Existing approaches to analyzing inequality often ignore within-group
inequality by solely analyzing mean differences between groups.
Approaches that do allow examining both changes in within- and
between-group inequality (e.g., [Western & Bloome
2009](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1467-9531.2009.01222.x)),
in turn, are limited in addressing causal questions about why inequality
is changing.

[Rosche (2022)](https://osf.io/preprints/socarxiv/f53kz/) introduces a
novel approach to analyzing how a treatment variable affects both
changes in within- and between-group inequality and decomposing these
changes into compositional and behavioral effects. The procedure
combines a classic variance decomposition with the
Kitagawa-Blinder-Oaxaca (KBO) decomposition approach. Compared to KBO,
however, the method allows analyzing treatment effects not only on the
mean but on the whole conditional distribution.

The ineqx packages implements both the descriptive ([Western & Bloome
2009](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1467-9531.2009.01222.x))
and causal variance decomposition ([Rosche
2022](https://osf.io/preprints/socarxiv/f53kz/)). The package allows
decomposing both the variance and the squared coefficient of variation
(CV²).

*With the ineqx package you can analyze*

- how overall inequality at a single point in time splits into a within-
  and between-group component (e.g., Does income inequality differ more
  within or between gender categories?)
- whether the overall change in inequality over time stems from changes
  in within-group inequality, between-group inequality, or changes in
  the composition of the groups (e.g., )
- the degree to which changes in inequality are due to changes in the
  effect of a treatment (i.e., binary predictor) on within- and
  between-group inequality, due to changes in the composition of the
  groups, and due to changes in pre-treatment inequality
- analyze the effect of a treatment on inequality (i.e. variability in
  an outcome) rather than just the mean

This is how the
[`ineqx()`](https://benrosche.github.io/ineqx/reference/ineqx.md)
function looks like:

``` r

data(incdat)
ineqx(treat="X", post="post_X", y="income", ystat="Var", group="female", time="year", ref=1990, dat)
```

### Try the ineqx package in this Shiny app

Settings

[Toggle navigation](_w_1c04c0467a4d4613b92d9722db50420f/#)

Causal variance decomposition

Number of groups

2 3 4

Number of timepoints

1 2 3 4

### Try the ineqx package!

Run ineqx

  
  

### Developers

I welcome contributions to the package! Feel free to submit changes for
review or contact me if you have any questions.

### Issues or Feature Requests

If you would like to log an issue or submit a feature request, please
create a new issue or comment on an existing issue on [GitHub
Issues](https://github.com/benrosche/ineqx/issues) on this repo.

### Changelog

See [NEWS.md](https://github.com/benrosche/ineqx/news/index.html) for
the package changelog.

More information can be found
[here](http://benrosche.com/projects/ineqx/)
