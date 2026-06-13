# Get temporal units

Generic function to extract temporal units from a supported object.

## Usage

``` r
units_time(x, ...)

# Default S3 method
units_time(x, ...)

# S3 method for class 'admove_tref'
units_time(x, ...)

# S3 method for class 'admove_cov'
units_time(x, ...)

# S3 method for class 'admove_tags'
units_time(x, ...)

# S3 method for class 'admove_data'
units_time(x, ...)
```

## Arguments

- x:

  An object from which to extract temporal units.

- ...:

  Further arguments passed to methods.

## Value

A character string describing the temporal units, such as `"day"`,
`"month"`, `"year"`, or `"week"`.
