# Set coordinate reference system (CRS)

Replacement function for
[`crs()`](https://tokami.github.io/admove/reference/crs.md).

## Usage

``` r
crs(x) <- value

# Default S3 method
crs(x) <- value

# S3 method for class 'admove_sref'
crs(x) <- value

# S3 method for class 'admove_grid'
crs(x) <- value

# S3 method for class 'admove_cov'
crs(x) <- value

# S3 method for class 'admove_tags'
crs(x) <- value

# S3 method for class 'admove_data'
crs(x) <- value
```

## Arguments

- x:

  An object to modify.

- value:

  A CRS specification.

## Value

`x` with updated CRS information.
