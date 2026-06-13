# Get spatial units

Generic function to extract spatial units from a supported object.

## Usage

``` r
units_space(x, ...)

# Default S3 method
units_space(x, ...)

# S3 method for class 'admove_sref'
units_space(x, ...)

# S3 method for class 'admove_grid'
units_space(x, ...)

# S3 method for class 'admove_cov'
units_space(x, ...)

# S3 method for class 'admove_tags'
units_space(x, ...)

# S3 method for class 'admove_data'
units_space(x, ...)
```

## Arguments

- x:

  An object from which to extract spatial units.

- ...:

  Further arguments passed to methods.

## Value

A character string describing the spatial units, such as `"m"`, `"km"`,
or `"degree"`.
