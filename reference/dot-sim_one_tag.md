# Simulate a single tag track

Simulates the movement trajectory of a single tagged individual under
the `admove` movement model. Depending on `sim_engine`, movement is
simulated either in continuous space using an Euler-type update or on a
discrete spatial grid using a continuous-time Markov chain (CTMC)
formulation.

## Usage

``` r
.sim_one_tag(
  conf,
  funcs,
  par,
  x0 = NULL,
  y0 = NULL,
  t0 = NULL,
  t1 = NULL,
  dt = 0.1,
  id = NULL,
  tag_type = "d",
  xygrid = NULL,
  nextTo = NULL,
  next_dist = NULL,
  xgr = NULL,
  ygr = NULL,
  celltable = NULL,
  xcen = NULL,
  ycen = NULL,
  sim_engine = 1,
  use_reject = FALSE,
  n_reject = 20,
  ctmc_method = 2
)
```

## Arguments

- conf:

  A configuration list controlling which movement components are active,
  typically produced by
  [`default_conf()`](https://tokami.github.io/admove/reference/default_conf.md).

- funcs:

  A named list of simulation functions, typically produced by
  [`default_sim_funcs()`](https://tokami.github.io/admove/reference/default_sim_funcs.md).
  These functions define taxis, diffusion, diffusion gradient, and
  advection.

- par:

  A parameter list containing model parameters used during simulation.

- x0:

  Numeric scalar giving the release x-coordinate.

- y0:

  Numeric scalar giving the release y-coordinate.

- t0:

  Numeric scalar giving the release time.

- t1:

  Numeric scalar giving the final time of the simulated track.

- dt:

  Time step used for simulation.

- id:

  Optional tag identifier. If `NULL`, a random identifier is generated.

- tag_type:

  Character string giving the tag type. Used to select the appropriate
  observation-error component.

- xygrid:

  Matrix or data frame of active grid-cell centres. Required for CTMC
  simulation.

- nextTo:

  Neighbour structure of the grid, as produced for example by
  [`get_neighbours()`](https://tokami.github.io/admove/reference/get_neighbours.md).
  Required for CTMC simulation.

- next_dist:

  Numeric vector of neighbour distances. Required for CTMC simulation.

- xgr:

  Numeric vector of x-direction grid breaks.

- ygr:

  Numeric vector of y-direction grid breaks.

- celltable:

  Matrix linking grid-cell positions to cell indices.

- xcen:

  Numeric vector of x-coordinates of grid-cell centres.

- ycen:

  Numeric vector of y-coordinates of grid-cell centres.

- sim_engine:

  Integer selecting the simulation engine: `1` for continuous-space
  simulation and `2` for CTMC-based grid simulation.

- use_reject:

  Logical; if `TRUE`, rejected diffusion proposals are redrawn when
  simulated moves leave the valid domain.

- n_reject:

  Maximum number of rejection attempts if `use_reject = TRUE`.

- ctmc_method:

  Integer controlling the matrix-exponential method used in CTMC
  simulation.

## Value

A data frame with simulated tag positions over time. The returned data
frame contains at least the columns `t`, `x`, `y`, and `id`.

## Details

For `sim_engine = 1`, movement is simulated in continuous space by
combining taxis, advection, and diffusion increments at each time step.

For `sim_engine = 2`, movement is simulated on the spatial grid by
constructing a transition-rate matrix from diffusion, taxis, and
advection, and then propagating the tag distribution forward over one
time step.

After simulation, observation error is added to all intermediate
positions, while the first and last positions are left unchanged.
