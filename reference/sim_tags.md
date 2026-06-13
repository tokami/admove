# Simulate tagging data

Simulates tagging data under the `admove` movement model for one of the
supported tag types: data-storage tags, mark-resight tags, or
conventional mark-recapture tags.

## Usage

``` r
sim_tags(
  tag_type,
  grid = NULL,
  cov = NULL,
  par = NULL,
  dat = NULL,
  conf = NULL,
  n_tags = 1,
  n_resightings = c(1, 5),
  trange = NULL,
  dt_tags = NULL,
  trange_rel = NULL,
  xrange_rel = NULL,
  yrange_rel = NULL,
  n_release_events = 1,
  release_events = NULL,
  trange_rec = NULL,
  knots_tax = NULL,
  knots_dif = NULL,
  funcs = NULL,
  use_reject = FALSE,
  n_reject = 20,
  sim_engine = 1,
  ctmc_method = 2,
  sref = NULL,
  tref = NULL,
  target_dif_frac = 1/500,
  target_tax_frac = 1/10,
  target_sdO_frac = 1/30,
  plot = FALSE,
  plot_land = FALSE,
  verbose = TRUE
)
```

## Arguments

- tag_type:

  Character string specifying the tag type to simulate. Supported values
  are `"d"` or `"dtags"` for data-storage tags, `"s"` or `"stags"` for
  mark-resight tags, and `"c"` or `"ctags"` for conventional
  mark-recapture tags.

- grid:

  Optional spatial grid used for simulation. If not already of class
  `"admove_grid"`, it is converted using
  [`create_grid()`](https://tokami.github.io/admove/reference/create_grid.md).

- cov:

  Optional covariate fields used to define spatially and temporally
  varying movement rates.

- par:

  Optional named list of simulation parameters. Missing parameters are
  filled using
  [`default_sim_par()`](https://tokami.github.io/admove/reference/default_sim_par.md).

- dat:

  Optional `admove_data` object. If `NULL`, a default data object is
  constructed internally from the supplied inputs.

- conf:

  Optional configuration list controlling the movement model. If `NULL`,
  [`default_conf()`](https://tokami.github.io/admove/reference/default_conf.md)
  is used.

- n_tags:

  Number of tags to simulate.

- n_resightings:

  Integer vector giving the minimum and maximum number of resightings
  for mark-resight tags.

- trange:

  Numeric vector of length 2 giving the simulation time range. If
  `NULL`, it is derived from the covariate time dimension when possible,
  or defaults to `c(0, 1)`.

- dt_tags:

  Optional simulation time step for tag trajectories. If `NULL`, a
  default value of `0.1` is used.

- trange_rel:

  Optional time range within which release times are generated.

- xrange_rel:

  Optional x-range within which release positions are generated.

- yrange_rel:

  Optional y-range within which release positions are generated.

- n_release_events:

  Number of release events to simulate if `release_events` is not
  supplied.

- release_events:

  Optional data frame specifying release events. If `NULL`, release
  events are generated internally.

- trange_rec:

  Optional time range within which recapture or final observation times
  are generated.

- knots_tax:

  Optional knot locations for the taxis component.

- knots_dif:

  Optional knot locations for the diffusion component.

- funcs:

  Optional named list of simulation functions. If `NULL`, defaults are
  created with
  [`default_sim_funcs()`](https://tokami.github.io/admove/reference/default_sim_funcs.md).

- use_reject:

  Logical; if `TRUE`, invalid movement proposals are rejected and
  resampled.

- n_reject:

  Maximum number of rejection attempts per step if `use_reject = TRUE`.

- sim_engine:

  Integer specifying the simulation engine: `1` for continuous-space
  simulation and `2` for CTMC-based grid simulation.

- ctmc_method:

  Integer controlling the matrix-exponential method used for CTMC
  simulation.

- sref:

  Optional spatial reference to attach to the simulated data.

- tref:

  Optional temporal reference to attach to the simulated data.

- target_dif_frac:

  Target diffusion strength as a fraction of the characteristic spatial
  scale squared per unit time, used by
  [`default_sim_par()`](https://tokami.github.io/admove/reference/default_sim_par.md).

- target_tax_frac:

  Target taxis strength as a fraction of the characteristic spatial
  scale per unit time, used by
  [`default_sim_par()`](https://tokami.github.io/admove/reference/default_sim_par.md).

- target_sdO_frac:

  Target observation error as a fraction of the characteristic spatial
  scale, used by
  [`default_sim_par()`](https://tokami.github.io/admove/reference/default_sim_par.md).

- plot:

  Logical; if `TRUE`, the simulated tags are plotted.

- plot_land:

  Logical; if `TRUE`, land masses are added to the plot.

- verbose:

  Logical; if `TRUE`, informative messages are printed.

## Value

An object of class `"admove_sim"` containing the simulated tags, grid,
covariates, simulation parameters, data object, and default fitting
components.

## Details

The function first sets up the spatial grid, time range, covariates, and
model configuration required for simulation. If no data object is
supplied, one is created internally from the provided inputs.

Tag trajectories are then simulated from generated or user-supplied
release events. The full simulated trajectories are retained for
data-storage tags. For mark-resight tags, only the release and a subset
of subsequent observations are kept. For conventional mark-recapture
tags, only the release and final observation are retained.

The returned object contains both the simulated tags and the associated
simulation setup, including the grid, covariates, simulation parameters,
and internally constructed `admove_data` object.

## Examples

``` r
data(skjepo)
sim <- sim_tags("ctags", skjepo$grid)
#> Consider providing time 'units' (tref$units) and/or 'period' (tref$period).
#> Consider providing time 'units' (tref$units) and/or 'period' (tref$period).
```
