# Get temporal period

Generic function to extract the temporal period from a supported object.

## Usage

``` r
period(x, ...)

# Default S3 method
period(x, ...)

# S3 method for class 'admove_tref'
period(x, ...)

# S3 method for class 'admove_cov'
period(x, ...)

# S3 method for class 'admove_tags'
period(x, ...)

# S3 method for class 'admove_data'
period(x, ...)
```

## Arguments

- x:

  An object from which to extract the period.

- ...:

  Further arguments passed to methods.

## Value

A numeric value giving the number of time steps per year, or another
stored period representation.
