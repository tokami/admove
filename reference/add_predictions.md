# Add model predictions to a fitted admove model

Computes and adds model predictions for a fitted `admove` model.

## Usage

``` r
add_predictions(fit)
```

## Arguments

- fit:

  A fitted model object of class `"admove"`, as returned by
  [`admove()`](https://tokami.github.io/admove/reference/admove.md).

## Value

An updated object of class `"admove"` with model predictions added.

## Details

This function evaluates the model prediction step and stores the
resulting predicted quantities in the fitted object.
