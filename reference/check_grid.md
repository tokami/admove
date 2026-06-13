# Check or coerce an admove grid

Checks whether an object is a valid `admove_grid`. If the object is not
already of class `"admove_grid"`, the function attempts to convert it
using
[`create_grid()`](https://tokami.github.io/admove/reference/create_grid.md).

## Usage

``` r
check_grid(x)
```

## Arguments

- x:

  An object to be checked or converted.

## Value

An object of class `"admove_grid"`, or `NULL` if `x` is `NULL`.

## Details

If `x` is `NULL`, it is returned unchanged. If `x` is not an
`admove_grid`, the function tries to construct one from `x` using
[`create_grid()`](https://tokami.github.io/admove/reference/create_grid.md).
If this fails, an error is thrown.
