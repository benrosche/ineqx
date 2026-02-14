# plot function

\[...\]

## Usage

``` r
# S3 method for class 'ineqx'
plot(ineqx.out, type, yscale = 1)
```

## Arguments

- ineqx.out:

  ineqx.out object from ineqx()

- type:

  Character string. Plot type. Choose from: dMuP, dMuT, dSigmaP,
  dSigmaT, dWP, dWT, dBP, dBT, dPP, dPT, dT, dPA

- yscale:

  Either 1 or 2. Choose value on y-axis. 1: Effect in units of ystat. 2:
  Effect as % of value at reference time (ref)

## Value

Returns a ggplot2 object

## Details

The y-axis can be scaled in two ways:

`yscale=1` (default): absolute points of Var or CV2.  
`yscale=2`: % of the value of Var or CV2 at reference time.  
(Total Var or CV2 with `type="dT"` and within-group & between-group Var
or CV2 with `type="dW"` & `type="dB"`, respectively.)

## Author

Benjamin Rosche \<benjamin.rosche@gmail.com\>

## Examples

``` r
data(incdat)
i1 <- ineqx(y=inc, group=SES, time=year, ref=1, dat=incdat)
#> Error: object 'inc' not found
plot(i1, type="dMuSigma2")
#> Error: object 'i1' not found
```
