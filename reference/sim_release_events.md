# Simulate release events

Simulates release locations and release times for tagging experiments on
an `admove_grid`.

## Usage

``` r
sim_release_events(
  grid,
  trange_rel = NULL,
  xrange_rel = NULL,
  yrange_rel = NULL,
  n_release_events = 10,
  use_reject = TRUE,
  n_reject = 100
)
```

## Arguments

- grid:

  A grid object of class `"admove_grid"`, as returned by
  [`create_grid()`](https://tokami.github.io/admove/reference/create_grid.md).

- trange_rel:

  Optional numeric vector of length 2 giving the time range within which
  release times are generated. If `NULL`, the default is `c(0, 1)`.

- xrange_rel:

  Optional numeric vector of length 2 giving the x-range within which
  release locations are generated. If `NULL`, the full x-range of `grid`
  is used.

- yrange_rel:

  Optional numeric vector of length 2 giving the y-range within which
  release locations are generated. If `NULL`, the full y-range of `grid`
  is used.

- n_release_events:

  Number of release events to simulate.

- use_reject:

  Logical; if `TRUE`, candidate release locations that fall into invalid
  grid cells are rejected and redrawn.

- n_reject:

  Maximum number of rejection attempts for each release event.

## Value

A numeric matrix with columns `x0`, `y0`, and `t0`, giving the simulated
release positions and release times.

## Details

Release positions are drawn uniformly from the specified x- and
y-ranges, and release times are drawn uniformly from `trange_rel`. If
`use_reject = TRUE`, locations falling into `NA` cells of the grid are
rejected and resampled.

## Examples

``` r
rel_events <- sim_release_events(create_grid())
```
