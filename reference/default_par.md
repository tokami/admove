# Default initial parameter values for admove

Creates a named list of model parameters with default initial values for
estimation in `admove`.

## Usage

``` r
default_par(dat, conf = NULL, verbose = TRUE)
```

## Arguments

- dat:

  A data list containing model input data, as produced by
  [`setup_data()`](https://tokami.github.io/admove/reference/setup_data.md).

- conf:

  An optional configuration list, typically created by
  [`default_conf()`](https://tokami.github.io/admove/reference/default_conf.md).
  If `NULL`, a default configuration is generated from `dat`.

- verbose:

  Logical; if `TRUE`, informative messages are printed.

## Value

A named list of initial parameter values.

## Details

This function generates a named list of initial parameter values based
on the supplied data and model configuration. These values can be used
as starting values for model fitting and are intended to provide
reasonable defaults for the selected model setup.

The taxis scaling parameter `logKappa` is fixed during estimation (see
[`default_map()`](https://tokami.github.io/admove/reference/default_map.md))
and therefore its initial value is its final value. Because `kappa` has
units of \\\[\text{distance}^2 / \text{time}\]\\, a value of 1 is only
appropriate when coordinates are already on a unit scale. For projected
coordinates such as UTM (metres), `kappa = 1` makes the taxis
contribution negligible and renders the taxis spline coefficients
unidentifiable. The default is therefore set to
`kappa = cellsize^2 / median_dt`, which ensures that a unit covariate
gradient over one grid-cell width produces movement of one cell width
per median time step. Override via `par$logKappa <- log(<value>)` after
calling this function.

## Examples

``` r
par <- with(skjepo$sim, default_par(dat, conf))
```
