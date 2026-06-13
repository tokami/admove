# Simulate a complete admove data set

Simulates a complete data set for `admove`, including a spatial grid,
optional covariate fields, simulated tags, and the corresponding
`admove_data` object required for model fitting.

## Usage

``` r
sim_data(
  grid = NULL,
  cov = NULL,
  par = NULL,
  dat = NULL,
  conf = NULL,
  fit = NULL,
  trange = NULL,
  dt = NULL,
  simulate_cov = FALSE,
  simple = FALSE,
  nt = NULL,
  rho_t = 0.85,
  sd = 2,
  h = 0.2,
  nu = 2,
  rho_s = 0.8,
  delta = 0.1,
  zrange = c(20, 28),
  matern = TRUE,
  sim_buffer = TRUE,
  knots_tax = NULL,
  knots_dif = NULL,
  release_events = NULL,
  n_release_events = 10,
  trange_rel = NULL,
  xrange_rel = NULL,
  yrange_rel = NULL,
  trange_rec = NULL,
  use_dtags = TRUE,
  use_ctags = TRUE,
  use_stags = FALSE,
  n_dtags = 10,
  n_stags = 0,
  n_ctags = 100,
  n_resightings = c(1, 5),
  sim_engine = 1,
  use_reject = FALSE,
  n_reject = 100,
  target_dif_frac = 1/500,
  target_tax_frac = 1/10,
  target_sdO_frac = 1/30,
  plot = FALSE,
  verbose = TRUE
)
```

## Arguments

- grid:

  Optional spatial grid used for simulation. If not already an
  `admove_grid`, it is converted using
  [`create_grid()`](https://tokami.github.io/admove/reference/create_grid.md).

- cov:

  Optional covariate fields. If `NULL`, or if `simulate_cov = TRUE`,
  covariates are simulated internally using
  [`sim_cov()`](https://tokami.github.io/admove/reference/sim_cov.md).

- par:

  Optional named list of simulation parameters. Missing parameters are
  filled using
  [`default_sim_par()`](https://tokami.github.io/admove/reference/default_sim_par.md).

- dat:

  Optional `admove_data` object used as a template for the simulation.
  If supplied, missing inputs such as grid, covariates, time range, and
  spline knots are extracted from it.

- conf:

  Optional configuration list for the movement model.

- fit:

  Optional fitted `admove` model. If supplied, the simulation can reuse
  information from the fitted model, such as the original data object.

- trange:

  Numeric vector of length 2 giving the simulation time range. If
  `NULL`, a default range is used or extracted from `dat`.

- dt:

  Optional model time-step size.

- simulate_cov:

  Logical; if `TRUE`, covariate fields are simulated even when
  covariates are already available from `dat` or `fit`.

- simple:

  Logical; if `TRUE`, simple deterministic covariate fields are
  simulated instead of random fields.

- nt:

  Number of covariate fields (time steps) to simulate.

- rho_t:

  Temporal autocorrelation coefficient for simulated covariate fields.

- sd:

  Standard deviation used in the simulation of covariate fields.

- h:

  Parameter controlling the covariate precision or covariance structure.

- nu:

  Smoothness parameter of the Matérn covariance structure.

- rho_s:

  Spatial range parameter of the Matérn covariance structure.

- delta:

  Small positive value added for numerical stability in the precision
  matrix.

- zrange:

  Numeric vector of length 2 giving the target range of the simulated
  covariate values.

- matern:

  Logical; if `TRUE`, a Matérn-based covariance structure is used for
  covariate simulation.

- sim_buffer:

  Logical; if `TRUE`, the simulation grid is extended by a one-cell
  buffer before covariates are simulated.

- knots_tax:

  Optional knot locations for the taxis component.

- knots_dif:

  Optional knot locations for the diffusion component.

- release_events:

  Optional data frame or matrix of release events. If `NULL`, release
  events are simulated internally using
  [`sim_release_events()`](https://tokami.github.io/admove/reference/sim_release_events.md).

- n_release_events:

  Number of release events to simulate if `release_events` is not
  supplied.

- trange_rel:

  Optional time range within which releases occur.

- xrange_rel:

  Optional x-range within which release locations are drawn.

- yrange_rel:

  Optional y-range within which release locations are drawn.

- trange_rec:

  Optional time range within which recapture or final observation times
  are drawn.

- use_dtags:

  Logical; if `TRUE`, data-storage tags are simulated.

- use_ctags:

  Logical; if `TRUE`, conventional mark-recapture tags are simulated.

- use_stags:

  Logical; if `TRUE`, mark-resight tags are simulated.

- n_dtags:

  Number of data-storage tags to simulate.

- n_stags:

  Number of mark-resight tags to simulate.

- n_ctags:

  Number of conventional mark-recapture tags to simulate.

- n_resightings:

  Integer vector giving the minimum and maximum number of resightings
  for mark-resight tags.

- sim_engine:

  Integer specifying the simulation engine: `1` for continuous-space
  simulation and `2` for CTMC-based grid simulation.

- use_reject:

  Logical; if `TRUE`, invalid simulated locations are rejected and
  redrawn where relevant.

- n_reject:

  Maximum number of rejection attempts used in rejection-based
  simulation steps.

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

  Logical; if `TRUE`, a summary plot of the simulated data is produced.

- verbose:

  Logical; if `TRUE`, informative messages are printed.

## Value

An object of class `"admove_sim"` containing the simulated grid,
covariates, simulated tags, simulation parameters, `admove_data` object,
and default fitting components.

## Details

The function provides a convenient way to generate a complete simulation
setup for `admove`. Depending on the supplied inputs, it can reuse an
existing grid, covariates, or fitted model, or simulate these components
from scratch.

Covariates are simulated with
[`sim_cov()`](https://tokami.github.io/admove/reference/sim_cov.md) if
needed. Tag release events are either provided directly or generated
with
[`sim_release_events()`](https://tokami.github.io/admove/reference/sim_release_events.md).
Tag data for the requested tag types are then simulated using
[`sim_tags()`](https://tokami.github.io/admove/reference/sim_tags.md).
Finally, the simulated tags and covariates are combined into an
`admove_data` object suitable for model fitting.

The returned object also includes default configuration, parameter, and
map objects for downstream fitting with
[`admove()`](https://tokami.github.io/admove/reference/admove.md).

## Examples

``` r
sim <- sim_data()
```
