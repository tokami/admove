# Get temporal origin

Generic function to extract the temporal origin from a supported object.

## Usage

``` r
origin(x, ...)

# Default S3 method
origin(x, ...)

# S3 method for class 'admove_tref'
origin(x, ...)

# S3 method for class 'admove_cov'
origin(x, ...)

# S3 method for class 'admove_tags'
origin(x, ...)

# S3 method for class 'admove_data'
origin(x, ...)
```

## Arguments

- x:

  An object from which to extract the temporal origin.

- ...:

  Further arguments passed to methods.

## Value

A `POSIXct` time origin, or another stored origin representation.
