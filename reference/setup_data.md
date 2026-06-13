# Set up input data for an admove model

Combine the main input components for an `admove` analysis into a single
object of class `admove_data`. The function checks and harmonises
spatial and temporal references, prepares covariates and tags, defines
spline knots, and constructs default prediction grids and time points
used during fitting and plotting.

## Usage

``` r
setup_data(
  grid = NULL,
  cov = NULL,
  tags = NULL,
  trange = NULL,
  knots_tax = NULL,
  knots_dif = NULL,
  sref = NULL,
  tref = NULL,
  transform_sref = FALSE,
  shift_tref = FALSE,
  verbose = TRUE
)
```

## Arguments

- grid:

  Optional grid object, typically of class `admove_grid`, as returned by
  [`create_grid()`](https://tokami.github.io/admove/reference/create_grid.md).
  A grid is required for some model engines and for spatial prediction
  and plotting.

- cov:

  Optional covariate object or list of covariates. Covariates are
  typically prepared with
  [`prep_cov()`](https://tokami.github.io/admove/reference/prep_cov.md).
  If a single covariate is supplied, it is coerced internally to a list.

- tags:

  Optional tag data, typically as returned by one or more of
  [`prep_ctags()`](https://tokami.github.io/admove/reference/prep_tags.md),
  [`prep_dtags()`](https://tokami.github.io/admove/reference/prep_tags.md),
  or
  [`prep_stags()`](https://tokami.github.io/admove/reference/prep_tags.md).

- trange:

  Optional numeric vector of length two giving the model time range. If
  `NULL`, the time range is inferred from available tags and covariates.

- knots_tax:

  Optional matrix of spline knots for the taxis preference functions. If
  `NULL`, default knots are chosen from covariate quantiles.

- knots_dif:

  Optional matrix of spline knots for the diffusion preference
  functions. If `NULL`, default knots are chosen from covariate
  quantiles.

- sref:

  Optional spatial reference to use as the target spatial reference for
  all inputs. If supplied, it should be coercible to an `admove_sref`
  object.

- tref:

  Optional time reference to use as the target time reference for all
  inputs. If supplied, it should be coercible to an `admove_tref`
  object.

- transform_sref:

  Logical; if `TRUE`, spatial components are transformed or rescaled to
  a common spatial reference where possible. If `FALSE` (default), all
  spatial references must already match.

- shift_tref:

  Logical; if `TRUE`, temporal components are shifted or harmonised to a
  common time reference where possible. If `FALSE` (default), all time
  references must already match.

- verbose:

  Logical; if `TRUE`, informative messages are printed during
  processing. Default is `TRUE`.

## Value

An object of class `admove_data`, containing the processed grid,
covariates, tags, spline knots, prediction settings, and associated
spatial and temporal reference information.

## Details

The function first determines a common spatial reference (`sref`) and
time reference (`tref`) either from the user-supplied targets or from
the input objects. If `transform_sref = FALSE` or `shift_tref = FALSE`,
all inputs must already be compatible. Otherwise, inputs are harmonised
to the chosen target references.

If `trange` is not supplied, it is inferred from the union of tag and
covariate time ranges. If no valid time range can be determined, a
default range of `c(0, 1)` is used.

If spline knots are not supplied, default knots are generated from the
marginal covariate distributions. Taxis knots use three quantiles per
covariate, while diffusion knots use one quantile per covariate.

The returned object also contains default prediction components in
`dat$pred`, including:

- `pred$time`:

  A sequence of 10 prediction time points over `trange`.

- `pred$cov`:

  Prediction ranges for each covariate.

- `pred$grid`:

  The prediction grid, if a grid was supplied.

Additional defaults used later during fitting are also stored in the
output, including `eps`, `var_init_kf`, `log2steps`, `min_dt`, `dt`, and
`p_init`.

## Examples

``` r
ctags <- prep_tags(
  skjepo$ctags,
  tag_type = "c",
  names = c(
    t0 = "date_time", t1 = "date_caught",
    x0 = "rel_lon",   x1 = "recap_lon",
    y0 = "rel_lat",   y1 = "recap_lat"
  ),
  date_origin = "1899-12-30")
#> tref (time origin and units) was inferred from dates. Please check and adjust if needed.

## prepare data-storage tags
dtags <- prep_tags(
  skjepo$dtags,
  tag_type = "d",
  names = c(t = "time", x = "mptlon", y = "mptlat"),
  date_origin = "1899-12-30")
#> ID not specified, using order of list elements.
#> tref (time origin and units) was inferred from dates. Please check and adjust if needed.

grid <- create_grid(x = skjepo$grid, cellsize = c(10, 10))

cov <- prep_cov(skjepo$cov)

dat <- setup_data(
  grid = grid,
  cov = cov,
  tags = c(dtags, ctags),
  transform_sref = TRUE,
  shift_tref = TRUE
)
#> tref origins differ; shifting stored time values to match origin of 'tref'.
```
