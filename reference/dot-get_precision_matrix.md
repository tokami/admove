# Construct a precision matrix for a grid

Constructs a spatial precision matrix for an `admove_grid`, either from
a Matérn covariance structure or from a simple neighbour-based
structure.

## Usage

``` r
.get_precision_matrix(grid, h, nu, rho, delta, matern = TRUE)
```

## Arguments

- grid:

  A grid object of class `"admove_grid"`, as returned by
  [`create_grid()`](https://tokami.github.io/admove/reference/create_grid.md).

- h:

  Distance argument used internally in the Matérn covariance structure.

- nu:

  Smoothness parameter of the Matérn covariance structure.

- rho:

  Spatial range parameter of the Matérn covariance structure.

- delta:

  Small positive value added to the diagonal to ensure numerical
  stability and positive definiteness.

- matern:

  Logical; if `TRUE`, a Matérn-based covariance structure is used.
  Otherwise, a neighbour-based graph Laplacian is constructed.

## Value

A precision matrix for the active cells of the grid.

## Details

When `matern = TRUE`, pairwise distances between grid-cell centres are
used to build a Matérn covariance matrix, which is then inverted to
obtain a precision matrix. When `matern = FALSE`, a sparse
neighbour-based precision matrix is constructed from the grid adjacency
structure.

In both cases, `delta * I` is added to the diagonal.
