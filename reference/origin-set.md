# Set temporal origin

Replacement function for
[`origin()`](https://tokami.github.io/admove/reference/origin.md).

## Usage

``` r
origin(x) <- value

# Default S3 method
origin(x) <- value

# S3 method for class 'admove_tref'
origin(x) <- value

# S3 method for class 'admove_cov'
origin(x) <- value

# S3 method for class 'admove_tags'
origin(x) <- value

# S3 method for class 'admove_data'
origin(x) <- value
```

## Arguments

- x:

  An object to modify.

- value:

  A temporal origin, typically of class `POSIXct`.

## Value

`x` with updated temporal origin information.
