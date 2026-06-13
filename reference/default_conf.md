# Default model configuration for admove

Creates a configuration list with default settings for the `admove`
model based on the supplied input data.

## Usage

``` r
default_conf(dat, verbose = TRUE)
```

## Arguments

- dat:

  A data list containing model input data, as produced by
  [`setup_data()`](https://tokami.github.io/admove/reference/setup_data.md).

- verbose:

  Logical; if `TRUE`, informative messages may be printed.

## Value

A named list of default model configuration settings.

## Details

The default configuration is determined from the available data. In
particular, the function detects which tag types are present and sets
corresponding model flags.

The returned list contains logical flags controlling which data sources
and movement components are used, settings for observation-variance
estimation, the estimation engine, the CTMC approximation method, and
default seasonal settings for covariates and spline effects.

Observation uncertainty is **off by default**
(`obs_var_type = c(0L, 0L, 0L)`), meaning tag locations are treated as
exact. To estimate observation variance for a tag type, set the
corresponding element to `1L` (all but last observation) or `2L` (all
observations). For example, to estimate observation variance for
data-storage tags: `conf$obs_var_type[1] <- 1L`.

## Examples

``` r
conf <- default_conf(skjepo$dat)
```
