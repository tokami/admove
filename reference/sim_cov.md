# Simulate covariate fields

Simulates one or more spatial covariate fields on an `admove_grid`,
optionally with temporal autocorrelation.

## Usage

``` r
sim_cov(
  grid = NULL,
  nt = 1,
  simple = FALSE,
  rho_t = 0.85,
  sd = 2,
  h = 0.2,
  nu = 2,
  rho_s = NULL,
  delta = 0.1,
  zrange = c(20, 28),
  trange = NULL,
  matern = TRUE,
  sim_buffer = FALSE,
  tref = NULL,
  verbose = TRUE
)
```

## Arguments

- grid:

  Optional spatial grid. If not already an `admove_grid`, it is
  converted using
  [`create_grid()`](https://tokami.github.io/admove/reference/create_grid.md).

- nt:

  Number of covariate fields (time steps) to simulate.

- simple:

  Logical; if `TRUE`, a simple deterministic spatial field is simulated
  instead of a random field.

- rho_t:

  Temporal autocorrelation coefficient between successive covariate
  fields.

- sd:

  Standard deviation used in the simulation of the spatial random field.

- h:

  Parameter controlling the precision matrix or covariance structure.

- nu:

  Smoothness parameter of the Matérn covariance structure.

- rho_s:

  Spatial range parameter of the Matérn covariance structure. If `NULL`,
  a default is chosen from the grid resolution.

- delta:

  Small positive value added to the precision matrix diagonal to improve
  numerical stability.

- zrange:

  Numeric vector of length 2 giving the target range to which the
  simulated fields are rescaled.

- trange:

  Optional numeric vector of length 2 giving the time range covered by
  the simulated covariate series. If `NULL`, the range `c(0, nt - 1)` is
  used.

- matern:

  Logical; if `TRUE`, a Matérn-based covariance structure is used.
  Otherwise, a simpler neighbour-based precision matrix is used.

- sim_buffer:

  Logical; if `TRUE`, the grid is extended by a one-cell buffer before
  simulation.

- tref:

  Optional temporal reference attached to the returned covariate object.

- verbose:

  Logical; if `TRUE`, informative messages may be printed.

## Value

An `admove_cov` object containing the simulated covariate fields.

## Details

The first covariate field is simulated independently. If `nt > 1`,
subsequent fields are generated using an AR(1)-type temporal dependence
structure with coefficient `rho_t`. All fields are then rescaled to the
interval given by `zrange`.

The returned object is processed with
[`prep_cov()`](https://tokami.github.io/admove/reference/prep_cov.md)
and includes spatial and temporal metadata.

## Examples

``` r
cov <- sim_cov()
#> Consider providing time 'units' (tref$units) and/or 'period' (tref$period).
```
