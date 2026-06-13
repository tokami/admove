# Change coordinate units by scaling an admove object

Multiply the stored coordinates of an `admove_*` object by `scale` and
update its spatial reference accordingly. The spatial reference is
updated by:

- multiplying `sref(x)$crs_scale` by `scale`

- updating `sref(x)$units` where possible (e.g. m \<-\> km)

## Usage

``` r
scale_sref(x, scale = 1, units = NULL, verbose = TRUE)
```

## Arguments

- x:

  An `admove_grid`, `admove_*tags`, `admove_cov`, or `admove_data`
  object.

- scale:

  Numeric scalar. Coordinates are multiplied by this factor. For
  example, to convert meters to kilometers (numerically), use
  `scale = 0.001`.

- units:

  units

- verbose:

  Logical; if `TRUE`, prints informational messages.

## Value

`x` with rescaled coordinates and updated `sref(x)`.

## Examples

``` r
grid <- create_grid(cellsize = 5e3, xrange = c(0, 5e4), yrange = c(0, 5e4))

grid <- add_sref(grid, list(crs = 32631, units = "m", crs_scale = 1))

grid <- scale_sref(grid, scale = 0.001)

units_space(grid)
#> [1] "km"

crs_scale(grid)
#> [1] 0.001

```
