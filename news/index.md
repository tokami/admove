# Changelog

## admove 0.1.0

Initial release.

### Key features

- Fits advection–diffusion movement models to tagging data using
  [RTMB](https://github.com/kaskr/RTMB)-based automatic differentiation;
  no compiled C++ code is required.

- Supports three tag types: archival/data-logging tags
  ([`prep_dtags()`](https://tokami.github.io/admove/reference/prep_tags.md)),
  mark–recapture tags
  ([`prep_ctags()`](https://tokami.github.io/admove/reference/prep_tags.md)),
  and mark–resight tags
  ([`prep_stags()`](https://tokami.github.io/admove/reference/prep_tags.md)).
  Multiple tag types can be fitted jointly.

- Two estimation engines selectable via `conf$engine`:

  - **Kalman filter** (engine 1): continuous space, discrete time.
  - **Continuous-time Markov chain** (engine 2): discrete space,
    continuous time, using matrix exponentiation of the transition-rate
    matrix.

- Habitat preference functions for taxis (directed movement), diffusion
  (random movement intensity), and advection are represented by flexible
  polynomial splines fitted to covariate fields.

- [`setup_data()`](https://tokami.github.io/admove/reference/setup_data.md)
  assembles covariate grids, tag observations, and prediction grids into
  a single data list for model fitting.

- [`default_conf()`](https://tokami.github.io/admove/reference/default_conf.md),
  [`default_par()`](https://tokami.github.io/admove/reference/default_par.md),
  and
  [`default_map()`](https://tokami.github.io/admove/reference/default_map.md)
  provide sensible starting configurations; all can be modified before
  calling
  [`admove()`](https://tokami.github.io/admove/reference/admove.md).

- Post-processing helpers
  [`add_sdreport()`](https://tokami.github.io/admove/reference/add_sdreport.md),
  [`add_report()`](https://tokami.github.io/admove/reference/add_report.md),
  and
  [`add_predictions()`](https://tokami.github.io/admove/reference/add_predictions.md)
  attach uncertainty estimates, RTMB-reported quantities, and grid-level
  predictions to the fitted object.

- Simulation framework
  ([`sim_data()`](https://tokami.github.io/admove/reference/sim_data.md),
  [`sim_tags()`](https://tokami.github.io/admove/reference/sim_tags.md),
  [`sim_cov()`](https://tokami.github.io/admove/reference/sim_cov.md))
  for generating synthetic datasets under specified movement parameters.

- Spatial and temporal reference system attached to all major objects
  via `sref` / `tref` attributes and corresponding S3 generics.

- Visualisation:
  [`plot_fit()`](https://tokami.github.io/admove/reference/plot_fit.md),
  [`plot_compare()`](https://tokami.github.io/admove/reference/plot_compare.md),
  [`plot_pref_func()`](https://tokami.github.io/admove/reference/plot_pref_func.md),
  [`plot_pref_grid()`](https://tokami.github.io/admove/reference/plot_pref_grid.md),
  [`plot_taxis()`](https://tokami.github.io/admove/reference/plot_taxis.md),
  [`plot_diffusion()`](https://tokami.github.io/admove/reference/plot_diffusion.md),
  [`plot_tags()`](https://tokami.github.io/admove/reference/plot_tags.md),
  [`plot_taxis()`](https://tokami.github.io/admove/reference/plot_taxis.md),
  [`plot_cov()`](https://tokami.github.io/admove/reference/plot_cov.md),
  and
  [`plot_land()`](https://tokami.github.io/admove/reference/plot_land.md).

- Introductory vignette and full function documentation.
