# ineqx ggplot2 theme

A custom ggplot2 theme for consistent styling of ineqx plots. Use as
`+ theme_ineqx` in ggplot2 pipelines.

## Usage

``` r
theme_ineqx
```

## Format

An object of class `theme` (inherits from
[`ggplot2::theme`](https://ggplot2.tidyverse.org/reference/theme.html),
`gg`, `S7_object`) of length 144.

## Examples

``` r
library(ggplot2)
ggplot(mtcars, aes(wt, mpg)) + geom_point() + theme_ineqx

```
