# Synthetic panel of women's weekly earnings

A simulated two-wave panel of women's weekly earnings, used to
demonstrate the descriptive and causal variance decompositions. Its
structure and variable names follow the Current Population Survey (CPS):
each woman is observed in two consecutive years, and some women
(`mother = 1`) have their first child between the two waves. All values
are simulated; no real survey records are included.

## Usage

``` r
data(earnings_panel)
```

## Format

A data frame with 15,200 rows (7,600 women x 2 waves) and 12 variables:

- id:

  Person identifier.

- year:

  Survey year (1982–2020).

- decade:

  Decade of the first wave (1980, 1990, 2000, 2010); the same for both
  waves of a person.

- post:

  Wave indicator: 0 = first wave, 1 = second wave (after the birth, for
  mothers).

- mother:

  1 if the woman has her first child between the two waves, 0 otherwise.

- earnweek:

  Weekly earnings in dollars (0 if not employed).

- earnwt:

  Earnings weight.

- educ:

  Education: factor with levels `"HS or less"`, `"Some college"`,
  `"BA+"`.

- ses:

  Socioeconomic status: factor with levels `"Low"`, `"Middle"`,
  `"High"`.

- age:

  Age in years.

- race:

  Race/ethnicity: factor with levels `"White"`, `"Black"`, `"Hispanic"`,
  `"Other"`.

- married:

  1 if married, 0 otherwise.

## Details

The data-generating process is stylized: motherhood lowers employment
and mean earnings and raises earnings dispersion, more so for
less-educated women, and these penalties shrink over time; the share of
mothers falls, the share with a BA rises, and within-group dispersion
grows. Draws are stratified so that results are smooth despite the
modest sample size. See `create/create_earnings_panel.R` in the package
source for the full data-generating process.

## Examples

``` r
data(earnings_panel)
head(earnings_panel)
#>   id year decade post mother earnweek earnwt       educ ses age  race married
#> 1  1 1982   1980    0      0   263.81   2846 HS or less Low  40 White       1
#> 2  1 1983   1980    1      0     0.00   2846 HS or less Low  41 White       1
#> 3  2 1982   1980    0      0     0.00   3108 HS or less Low  24 White       1
#> 4  2 1983   1980    1      0     0.00   3108 HS or less Low  25 White       1
#> 5  3 1982   1980    0      0     0.00   2395 HS or less Low  24 White       0
#> 6  3 1983   1980    1      0     0.00   2395 HS or less Low  25 White       0
```
