# Get spatial reference from an object

Generic function to extract spatial reference information from an
object. admove objects typically store this information in
`attr(x, "sref")`.

## Usage

``` r
# S3 method for class 'admove_cov_list'
sref(x, ...)

sref(x, ...)

# Default S3 method
sref(x, ...)
```

## Arguments

- x:

  An object.

- ...:

  Further arguments passed to methods (currently unused).

## Value

An object of class `admove_sref`.
