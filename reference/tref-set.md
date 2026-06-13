# Set time reference on an object

Replacement function for
[`tref()`](https://tokami.github.io/admove/reference/tref.md).

## Usage

``` r
# S3 method for class 'admove_cov_list'
tref(x) <- value

tref(x) <- value

# Default S3 method
tref(x) <- value
```

## Arguments

- x:

  An object to modify.

- value:

  An object of class `admove_tref`.

## Value

`x` with an updated `"tref"` attribute.
