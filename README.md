
<!-- README.md is generated from README.Rmd. Please edit that file -->

# <img src="man/figures/ineqx-hexagon.jpg" width="200" align="right" />

### ineqx: Descriptive and causal variance decompositions

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
(CV<sup>2</sup>).

<u>With the ineqx package you can analyze</u>

-   how overall inequality at a single point in time splits into a
    within- and between-group component (e.g., Does income inequality
    differ more within or between gender categories?)
-   whether the overall change in inequality over time stems from
    changes in within-group inequality, between-group inequality, or
    changes in the composition of the groups (e.g., )
-   the degree to which changes in inequality are due to changes in the
    effect of a treatment (i.e., binary predictor) on within- and
    between-group inequality, due to changes in the composition of the
    groups, and due to changes in pre-treatment inequality
-   analyze the effect of a treatment on inequality (i.e. variability in
    an outcome) rather than just the mean

This is how the `ineqx()` function looks like:

``` r
data(incdat)
ineqx(treat="X", post="post_X", y="income", ystat="Var", group="female", time="year", ref=1990, dat)
```

### Try the ineqx package in this Shiny app

<style>
iframe {
  width: 100%;
  height: 1000px;
  border: none;
}
</style>
<iframe src="https://benrosche.shinyapps.io/ineqx-app/" />

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
