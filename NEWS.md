# admove 0.1.0

Initial release.

## Key features

* Fits advection–diffusion movement models to tagging data using
  [RTMB](https://github.com/kaskr/RTMB)-based automatic differentiation;
  no compiled C++ code is required.

* Supports three tag types: archival/data-logging tags (`prep_dtags()`),
  mark–recapture tags (`prep_ctags()`), and mark–resight tags (`prep_stags()`).
  Multiple tag types can be fitted jointly.

* Two estimation engines selectable via `conf$engine`:
  - **Kalman filter** (engine 1): continuous space, discrete time.
  - **Continuous-time Markov chain** (engine 2): discrete space, continuous
    time, using matrix exponentiation of the transition-rate matrix.

* Habitat preference functions for taxis (directed movement), diffusion
  (random movement intensity), and advection are represented by flexible
  polynomial splines fitted to covariate fields.

* `setup_data()` assembles covariate grids, tag observations, and prediction
  grids into a single data list for model fitting.

* `default_conf()`, `default_par()`, and `default_map()` provide sensible
  starting configurations; all can be modified before calling `admove()`.

* Post-processing helpers `add_sdreport()`, `add_report()`, and
  `add_predictions()` attach uncertainty estimates, RTMB-reported quantities,
  and grid-level predictions to the fitted object.

* Simulation framework (`sim_data()`, `sim_tags()`, `sim_cov()`) for
  generating synthetic datasets under specified movement parameters.

* Spatial and temporal reference system attached to all major objects via
  `sref` / `tref` attributes and corresponding S3 generics.

* Visualisation: `plot_fit()`, `plot_compare()`, `plot_pref_func()`,
  `plot_pref_grid()`, `plot_taxis()`, `plot_diffusion()`, `plot_tags()`,
  `plot_taxis()`, `plot_cov()`, and `plot_land()`.

* Introductory vignette and full function documentation.
