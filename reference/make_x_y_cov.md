# Create x- and y-coordinate covariate fields

`make_x_y_cov()` creates two simple covariate fields representing the
spatial x and y coordinates of a grid. The function returns these as a
list of covariate arrays that can be used as input to *admove*, for
example as simple spatial trend covariates or for testing and
demonstration purposes.

The first covariate varies along the x dimension and the second varies
along the y dimension. Both covariates are created on the spatial domain
and cell centres of the supplied object.

## Usage

``` r
make_x_y_cov(x, tref = NULL)
```

## Arguments

- x:

  An object providing spatial dimensions and cell centres, typically an
  `admove_grid` or another object for which
  [`x_centers()`](https://tokami.github.io/admove/reference/x_centers.md),
  [`y_centers()`](https://tokami.github.io/admove/reference/y_centers.md),
  and [`sref()`](https://tokami.github.io/admove/reference/sref.md) are
  defined.

- tref:

  Optional time reference information to attach to the returned
  covariates. Default: `NULL`.

## Value

A list-like object of class `admove_cov` containing two covariate
fields: one for the x coordinate and one for the y coordinate.

## Details

Both covariates are created as single-time-slice fields with
`times = 0`. Spatial reference information is copied from `x`, and
optional temporal reference information can be attached via `tref`.

## Examples

``` r
## xy_cov <- make_x_y_cov(grid)
```
