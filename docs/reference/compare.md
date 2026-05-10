# Compare multiple ineqx results

Stacks two or more ineqx result objects into a single comparison object
with a scenario column. Useful for comparing observed vs counterfactual
decompositions, or multiple counterfactual scenarios.

## Usage

``` r
compare(..., names = NULL)
```

## Arguments

- ...:

  Named ineqx result objects, or a single named list of them. All
  objects must be the same class (e.g., all `ineqx_causal_longit`).

- names:

  Optional character vector of scenario labels. If provided, overrides
  names from `...`.

## Value

An object of class `ineqx_compare` containing:

- data:

  Long-format data.frame with a `scenario` column

- class_type:

  Character: the shared class of the input objects

- scenarios:

  Character vector of scenario names

- ystat:

  Character: shared inequality measure

- ref:

  Numeric: shared reference period (for longitudinal)

- objects:

  List of original input objects

## Examples

``` r
if (FALSE) { # \dontrun{
comp <- compare(Observed = res1, Counterfactual = res2)
print(comp)
plot(comp)
} # }
```
