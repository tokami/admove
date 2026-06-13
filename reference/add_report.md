# Add reported RTMB quantities to a fitted admove model

Runs `obj$report()` for a fitted `admove` model and adds the reported
quantities to the fitted object.

## Usage

``` r
add_report(fit)
```

## Arguments

- fit:

  A fitted model object of class `"admove"`, as returned by
  [`admove()`](https://tokami.github.io/admove/reference/admove.md).

## Value

An updated object of class `"admove"` with reported quantities added.

## Details

The reported quantities are stored in the `rep` component of the
returned object. If this component already exists, it is overwritten.
