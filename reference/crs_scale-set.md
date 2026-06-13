# Set CRS scale

Replacement function for
[`crs_scale()`](https://tokami.github.io/admove/reference/crs_scale.md).

## Usage

``` r
crs_scale(x) <- value

# Default S3 method
crs_scale(x) <- value

# S3 method for class 'admove_sref'
crs_scale(x) <- value

# S3 method for class 'admove_grid'
crs_scale(x) <- value

# S3 method for class 'admove_cov'
crs_scale(x) <- value

# S3 method for class 'admove_tags'
crs_scale(x) <- value

# S3 method for class 'admove_data'
crs_scale(x) <- value
```

## Arguments

- x:

  An object to modify.

- value:

  A numeric scaling factor.

## Value

`x` with updated CRS scaling information.
