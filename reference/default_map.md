# Default parameter map

Construct the default `map` list used to control which parameters are
estimated, fixed, or coupled in TMB. The returned list has the same
overall structure as the model parameter list and is intended to be
passed to
[`TMB::MakeADFun()`](https://rdrr.io/pkg/TMB/man/MakeADFun.html) via its
`map` argument.

By default:

- taxis spline coefficients (`alpha`) are mapped using
  `.make_alpha_map()`, with the first row fixed and remaining
  coefficients estimated independently;

- the taxis scaling parameter (`logKappa`) is fixed;

- diffusion spline coefficients (`beta`) are mapped using
  `.make_beta_map()`;

- advection coefficients (`gamma`) are either fixed or, if advection is
  enabled, coupled between the \\x\\- and \\y\\-directions within each
  covariate;

- observation-error parameters (`logSdO`) are fixed unless estimation is
  enabled for the corresponding tag type, in which case \\x\\- and
  \\y\\-direction standard deviations are coupled by default.

## Usage

``` r
default_map(dat, conf, par)
```

## Arguments

- dat:

  A data list as produced by
  [`setup_data()`](https://tokami.github.io/admove/reference/setup_data.md).

- conf:

  A configuration list as produced by
  [`default_conf()`](https://tokami.github.io/admove/reference/default_conf.md).

- par:

  A parameter list with initial values as produced by
  [`default_par()`](https://tokami.github.io/admove/reference/default_par.md).

## Value

A named list of factors with elements `alpha`, `logKappa`, `beta`,
`gamma`, and `logSdO`. Entries with `NA` are fixed, while equal factor
levels are estimated as the same parameter.

## Examples

``` r
map <- with(skjepo$sim, default_map(dat, conf, par))
```
