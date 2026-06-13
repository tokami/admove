# sf-compatible CRS accessor for admove objects

S3 method for
[`st_crs`](https://r-spatial.github.io/sf/reference/st_crs.html) that
returns the CRS stored in `sref(x)$crs` as an
[`sf::crs`](https://r-spatial.github.io/sf/reference/coerce-methods.html)
object.

## Usage

``` r
st_crs.admove_grid(x, ...)

st_crs.admove_cov(x, ...)

st_crs.admove_tags(x, ...)
```

## Arguments

- x:

  An object with spatial reference (`sref`).

- ...:

  Further arguments (unused).

## Value

An
[`sf::crs`](https://r-spatial.github.io/sf/reference/coerce-methods.html)
object.

## Details

Requires the sf package to be installed.
