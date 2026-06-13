# Set spatial reference on an object

Replacement function for
[`sref()`](https://tokami.github.io/admove/reference/sref.md).

## Usage

``` r
# S3 method for class 'admove_cov_list'
sref(x) <- value

sref(x) <- value

# Default S3 method
sref(x) <- value
```

## Arguments

- x:

  An object to modify.

- value:

  An object of class `admove_sref`.

## Value

`x` with an updated `"sref"` attribute.
