# Descriptive within/between decomposition

\[...\]

## Usage

``` r
wibe(
  y = NULL,
  group = NULL,
  time = NULL,
  weights = NULL,
  ref = F,
  long = F,
  dat
)
```

## Arguments

- y:

  Dependent variable

- group:

  Grouping variable to decompose variance into within- and between-group
  components

- time:

  Time variable to analyze change over time

- weights:

  Probability weights

- ref:

  Number or FALSE. Should values be reported in reference to a specific
  time?

- long:

  Logical. Should output be in long format?

- dat:

  Dataframe

- smoothDat:

  Logical. Should data be smoothed?

## Value

List of length 2. Element 1 returns the decomposition by group and time.
Elements 2 returns the decomposition by time.

## Details

...

## Author

Benjamin Rosche \<benjamin.rosche@gmail.com\>

## Examples

``` r
data(incdat)
wibe1 <- wibe(y="inc", group="SES", time="year", dat=dat1)
#> Error: object 'dat1' not found
```
