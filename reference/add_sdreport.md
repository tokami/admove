# Add an RTMB sdreport to a fitted admove model

Runs `RTMB::sdreport()` for a fitted `admove` model and adds the result
to the fitted object. Estimated values and standard deviations are also
extracted as named lists and stored in the returned object.

## Usage

``` r
add_sdreport(fit, save_covariance = FALSE)
```

## Arguments

- fit:

  A fitted model object of class `"admove"`, as returned by
  [`admove()`](https://tokami.github.io/admove/reference/admove.md).

- save_covariance:

  Logical; if `TRUE`, the full covariance matrix from the sdreport is
  retained. If `FALSE` (default), the covariance matrix is removed to
  reduce memory use.

## Value

An updated object of class `"admove"` with sdreport results added.

## Details

This function adds three components to the fitted object:

- `sdrep`: the full sdreport object returned by `RTMB::sdreport()`

- `pl`: a named list of estimates extracted with `as.list(sdrep, "Est")`

- `plsd`: a named list of standard deviations extracted with
  `as.list(sdrep, "Std")`

If these components already exist in `fit`, they are overwritten.
