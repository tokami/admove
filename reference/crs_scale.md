# Get CRS scale

Generic function to extract the scaling factor relating CRS units to the
working spatial units used by an object.

## Usage

``` r
crs_scale(x, ...)

# Default S3 method
crs_scale(x, ...)

# S3 method for class 'admove_sref'
crs_scale(x, ...)

# S3 method for class 'admove_grid'
crs_scale(x, ...)

# S3 method for class 'admove_cov'
crs_scale(x, ...)

# S3 method for class 'admove_tags'
crs_scale(x, ...)

# S3 method for class 'admove_data'
crs_scale(x, ...)
```

## Arguments

- x:

  An object.

- ...:

  Further arguments passed to methods.

## Value

A numeric scaling factor.
