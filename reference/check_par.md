# Check parameter dimensions for admove

Checks whether the parameters supplied in `par` have the expected names
and dimensions for the provided data and model configuration.

## Usage

``` r
check_par(par, dat, conf = NULL, verbose = TRUE)
```

## Arguments

- par:

  A named list of model parameters to be checked.

- dat:

  A data list containing model input data, as produced by
  [`setup_data()`](https://tokami.github.io/admove/reference/setup_data.md).

- conf:

  An optional configuration list, typically created by
  [`default_conf()`](https://tokami.github.io/admove/reference/default_conf.md).
  If `NULL`, a default configuration is generated from `dat`.

- verbose:

  Logical; if `TRUE`, informative messages are printed when `conf` is
  generated internally.

## Value

Invisibly returns `TRUE` if all parameter names and dimensions are
valid.

## Details

The function constructs the expected parameter structure using
[`default_par()`](https://tokami.github.io/admove/reference/default_par.md)
and compares it with the user-supplied `par` list. It checks that all
required parameters are present, flags unexpected parameters, and
verifies that each parameter has the correct dimensions.

If any mismatch is found, the function stops with an informative error
message describing the problem.

## Examples

``` r
## If all checks passed, the function returns an invisible TRUE:
with(skjepo$sim, check_par(par, dat, conf))

## If there is a problem, the function returns an error:
if (FALSE) { # \dontrun{
par <- with(skjepo$sim, default_par(dat, conf))
par$alpha <- matrix(0, 3, 1)
check_par(par, skjepo$sim$dat, skjepo$sim$conf)
} # }
```
