# Calculate AME

\[...\]

## Usage

``` r
calcAME(treat, group, time, what, vfr, dat)
```

## Arguments

- treat:

  Character string. Treatment variable.

- group:

  Grouping variable

- time:

  Time variable

- what:

  "Mu" or "Sigma"

- vfr:

  gamlss output

- dat:

  Dataframe

## Value

List of length 2. Element 1 returns AME by group and time. Elements 2
returns the AME by time.

## Author

Benjamin Rosche \<benjamin.rosche@gmail.com\>

## Examples

``` r
data(incdat)
vfr <- gamlss()
#> Error in terms(formula, specials = .gamlss.sm.list): promise already under evaluation: recursive default argument reference or earlier problems?
AME_mu <- calcAME("treat", group="SES", time="year", what="mu", vfr, incdat)
#> Error in calcAME("treat", group = "SES", time = "year", what = "mu", vfr,     incdat): treat not in dataset
```
