# Get time reference from an object

Generic function to extract time reference information from an object.
admove objects typically store this information in `attr(x, "tref")`.

## Usage

``` r
# S3 method for class 'admove_cov_list'
tref(x, ...)

tref(x, ...)

# Default S3 method
tref(x, ...)
```

## Arguments

- x:

  An object.

- ...:

  Further arguments passed to methods (currently unused).

## Value

An object of class `admove_tref`.
