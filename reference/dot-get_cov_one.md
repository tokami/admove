# Simulate a single covariate field

Simulates a single spatial covariate field on an `admove_grid`.

## Usage

``` r
.get_cov_one(grid, sd, h, nu, rho, delta, matern, simple)
```

## Arguments

- grid:

  A grid object of class `"admove_grid"`, as returned by
  [`create_grid()`](https://tokami.github.io/admove/reference/create_grid.md).

- sd:

  Standard deviation used in the simulation of the random field.

- h:

  Parameter controlling the precision matrix or covariance structure.

- nu:

  Smoothness parameter of the Matérn covariance structure.

- rho:

  Spatial range parameter of the Matérn covariance structure.

- delta:

  Small positive value added to the precision matrix diagonal to improve
  numerical stability.

- matern:

  Logical; if `TRUE`, a Matérn-based covariance structure is used.
  Otherwise, a simpler neighbour-based precision matrix is used.

- simple:

  Logical; if `TRUE`, a simple deterministic spatial field is generated
  instead of a random field.

## Value

A matrix containing one simulated covariate field on the grid.

## Details

If `simple = TRUE`, the simulated field is a smooth radial surface
centred on the middle of the grid. Otherwise, a random field is
generated from a precision matrix defined by
[`.get_precision_matrix()`](https://tokami.github.io/admove/reference/dot-get_precision_matrix.md).
