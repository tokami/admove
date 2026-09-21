# admove 0.1.7

## Behaviour changes

* `units_space(x) <- `, `crs(x) <- ` and `crs_scale(x) <- ` now **error** on
  objects that carry coordinates (`admove_grid`, `admove_cov`, `admove_tags`,
  `admove_data`). They used to overwrite a single field of the spatial
  reference, leaving the other two fields *and the coordinates themselves*
  stale, e.g. `units_space(grid) <- "km"` relabelled a metre grid as km while
  `crs_scale` stayed 1, and nothing downstream noticed. The error names the
  function that moves everything together: `scale_sref()`, `add_sref()` or
  `transform_sref()`.

  On a bare `admove_sref`, which has no coordinates, the same replacement
  functions now recompute the dependent fields instead of desynchronising:
  `units_space(sp) <- "km"` on a metre CRS sets `crs_scale` to 0.001.

* `create_sref()` no longer overrides an explicitly supplied `crs_scale`. Its
  default is now `NULL`, meaning "derive from the CRS unit and `units`";
  `units` and `crs_scale` are two views of the same thing, so supplying either
  one is enough and the other follows. Supplying both keeps both, with a
  warning when they disagree. Previously `create_sref(32631, units = "m",
  crs_scale = 0.001)` silently returned `crs_scale = 1`, which also discarded
  the value `add_sref()` and `transform_sref()` had just derived.

  Supplying a `crs_scale` without `units` now labels the units with the
  composite form `"metre_x_1e-06"` that `.in_m()` already understood, rather
  than taking the CRS's own unit.

## New features

* `print()` method for `admove_sref`, and `summary()` of an `admove_grid` now
  shows the full spatial reference. Both make explicit that the CRS unit and
  the stored-coordinate unit are allowed to differ:

  ```
  crs:           Azimuthal Equidistant [custom]
  datum:         World Geodetic System 1984
  crs units:     metre
  stored units:  km
  crs scale:     1 metre = 0.001 km
  ```

  `crs()` continues to return the CRS of the *unscaled* coordinates and
  deliberately does not reflect `units_space()`, stored coordinates are CRS
  coordinates times `crs_scale`, and everything that hands a coordinate to sf
  divides by `crs_scale` first. See `?crs`.

## Bug fixes


# admove 0.1.5

## New features

* Tags may now carry an **ambiguous final position**: a set of candidate
  locations, exactly one of which is the true one, with known probabilities.
  The usual case is a mark-recapture tag whose recapturing vessel is known but
  whose individual set is not, so that any of that vessel's fishing sets could
  be the recapture location, weighted by effort. Two optional columns express
  this — `event` (rows sharing an event are mutually exclusive alternatives)
  and `prob` (their probabilities, which must sum to 1 within an event).
  Candidates may differ in time as well as position.

  The likelihood contribution becomes the finite mixture
  `log sum_k prob_k * f(x_k, t_k | release)`, which is **exact in both
  engines**: because only the final observation may be ambiguous, no Kalman or
  CTMC update has to be propagated through the mixture, so there is no
  Gaussian-mixture posterior to approximate. The sum is accumulated with
  `RTMB::logspace_add()` so that a distant candidate cannot underflow the
  objective or its gradient.

  Omitting the columns means "no ambiguity" and reproduces the previous
  behaviour exactly.

* `prep_tags()` / `prep_ctags()` gain a `candidates` argument and accept `names`
  as a list, for reading the repeated-column layout that recapture-uncertainty
  tables come in (`date1`, `lat1`, `lon1`, `per1`, `date2`, ...):

  ```r
  prep_ctags(unc,
             names = c(id = "fish_id", t0 = "release_date",
                       x0 = "release_lon", y0 = "release_lat",
                       t1 = "date", x1 = "lon", y1 = "lat", p1 = "per"),
             candidates = 1:9, date_origin = "1899-12-30")
  ```

* New `add_candidates()` attaches candidate positions that arrive as a separate
  long table (one row per fishing set, keyed by tag) to tags that already have a
  single recapture.

* `sim_tags()` and `sim_data()` gain `n_candidates` and `candidate_sd` for
  simulating ambiguous recaptures.

* `summarise_tags()` reports how many tags have an ambiguous final position and
  the average number of candidates.

* `plot_tags()` draws the alternative recapture positions as a fan from the
  release, shaded and weighted by probability.

## Bug fixes

* `plot_tags()` indexed the start positions of the release→recovery segments
  without the mark-recapture subset, which drew the wrong segments whenever a
  panel mixed tag types (masked by the `by_tag_type = TRUE` default).

* `build_time()` could leave a zero-length time step in the integration grid
  when two observations shared a time.

* `plot_tag_dist()` failed with an opaque `invalid value specified for graphical
  parameter "mfrow"` when `fit$tag_dist` was an empty list rather than `NULL` —
  which is what `add_tag_dist()` leaves behind when it skips every tag it was
  given (for example when a tag is recaptured within the first time step). It
  now says so and suggests a finer `dt`. Selecting no tags via `select` /
  `n_tags` is reported too.

## Behaviour changes

* `conf$obs_var_type = 1` ("all but the last observation") and the duplicated-id
  warning for mark-recapture tags now count observation **events** rather than
  rows, so candidate positions of one ambiguous observation are treated as the
  single observation they represent.

* Starting values derived from tag displacements (`logKappa`, diffusion) skip
  the alternatives within an event, which would otherwise contribute
  zero-length time steps and displacements between candidates.


# admove 0.1.3

## New features

* The habitat preference functions (taxis, diffusion) are now built as **natural
  cubic splines** by default (`conf$smooth_method = "natural"`), replacing the
  previous global interpolating polynomial. The spline has local support, is
  twice continuously differentiable, and extrapolates linearly beyond the outer
  knots — far more robust in the covariate tails than the polynomial. The
  parameters are still the function values at the knots, so parameterisation and
  interpretation are unchanged. Set `conf$smooth_method = "poly"` to recover the
  legacy polynomial behaviour.

* `add_predictions()` gains `grid` and `time` arguments, allowing a fitted
  model to be predicted onto a new spatial grid or set of times without
  refitting. The supplied grid/time is validated (matching spatial reference,
  pruned to the fitted covariate coverage) and stored in the returned object.

* `admove()` gains a `do_tag_dist` argument. When `TRUE`, predicted location
  distributions are precomputed for all tags via `add_tag_dist()` and stored
  in `fit$tag_dist` for use by `plot_tag_dist()`.

* `summary()` now reports the model's seasonality configuration (seasonal
  period, number of seasons, and seasonal covariate/spline settings).

## Improvements and fixes

* Improved multi-panel layout of `plot_tag_dist()`, with a new `asp` argument
  controlling the target panel aspect ratio.

* Plotting fixes.

* R CMD check fixes.


# admove 0.1.2

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
