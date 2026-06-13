# Prepare covariate fields for admove

`prep_cov()` converts covariate data into the standard array-based
format used by *admove*. Covariates can be supplied as a list of
matrices, as a 2D or 3D array, or as raster-like objects. When the input
has multiple layers, the `layers` argument controls whether they are
treated as separate covariate fields (default, returns a named
`admove_cov_list`) or as consecutive time steps of a single covariate
(returns a single `admove_cov`).

The function can also attach spatial and temporal reference information,
convert date-like time labels to numeric model time, and optionally plot
the prepared covariate field.

## Usage

``` r
prep_cov(
  x,
  x_centers = NULL,
  y_centers = NULL,
  times = NULL,
  date_decimal = FALSE,
  date_format = NULL,
  date_origin = NULL,
  tz = "UTC",
  sref = NULL,
  tref = NULL,
  layers = NULL,
  plot = FALSE,
  plot_land = FALSE,
  strict = FALSE,
  verbose = TRUE
)
```

## Arguments

- x:

  Covariate data. Supported inputs include:

  - a `data.frame` with columns `x` and `y` (coordinates) and one
    additional column per covariate; each row is one grid cell. The
    covariate column names become the names of the returned
    `admove_cov_list`. This is the natural output of GIS exports and
    spatial model pipelines.

  - a list of matrices, typically one matrix per time step,

  - a 2D array or matrix, interpreted as a single time slice,

  - a 3D array with dimensions x, y, and time,

  - `RasterLayer`, `RasterBrick`, or `RasterStack` objects, and

  - `SpatRaster` objects.

- x_centers:

  Optional numeric vector giving x coordinates of cell centres.

- y_centers:

  Optional numeric vector giving y coordinates of cell centres.

- times:

  Optional vector giving the time values associated with the third
  dimension.

- date_decimal:

  Logical; if `TRUE`, interpret time labels as decimal years and convert
  them to model time. Default: `FALSE`.

- date_format:

  Optional character string passed to
  [`as.Date()`](https://rdrr.io/r/base/as.Date.html) to parse character
  time labels. Default: `NULL`.

- date_origin:

  Optional origin passed to
  [`as.Date()`](https://rdrr.io/r/base/as.Date.html) when time labels
  are stored numerically. Default: `NULL`.

- tz:

  Time zone used when converting dates. Default: `"UTC"`.

- sref:

  Optional spatial reference information to attach to the returned
  object.

- tref:

  Optional time reference information to attach to the returned object.

- layers:

  Character string controlling how layers in multi-layer inputs are
  interpreted. Applies to `SpatRaster`, `RasterBrick`, `RasterStack`,
  lists of matrices, and 3-D arrays. Use `"covariates"` to treat each
  layer as a separate covariate field; the function then returns a named
  list of class `admove_cov_list`. Use `"time"` to treat layers as
  consecutive time steps of a single covariate; the function then
  returns a single `admove_cov`. Default (`NULL`) is `"covariates"` for
  `SpatRaster`/`Raster*` objects and `"time"` for lists and arrays,
  which preserves backward-compatible behaviour for those types. Ignored
  for single-layer inputs.

- plot:

  Logical; if `TRUE`, plot the prepared covariate field. Default:
  `FALSE`.

- plot_land:

  Logical; passed to the plotting method when `plot = TRUE`. Default:
  `FALSE`.

- strict:

  Logical; if `TRUE`, require stricter matching and validation of
  dimension names. Default: `FALSE`.

- verbose:

  Logical; if `TRUE`, print informative messages. Default: `TRUE`.

## Value

When `layers = "time"` or when the input has only one layer, an object
of class `admove_cov`, stored as a 3D array with dimensions
corresponding to x, y, and time. When `layers = "covariates"` and the
input has multiple layers, a named list of class `admove_cov_list` with
one `admove_cov` element per layer, ready to pass directly to
[`setup_data()`](https://tokami.github.io/admove/reference/setup_data.md).

## Details

If `x` is two-dimensional, it is converted to a 3D array with a single
time slice. Dimension names are validated and, where possible, inferred
from the input or from the optional `x_centers`, `y_centers`, and
`times` arguments.

If date-like time labels are supplied, they can be converted to numeric
model time using `date_format`, `date_origin`, or `date_decimal`.
Spatial and temporal reference information are preserved from the input
where available and supplemented by `sref` and `tref` if provided. A
user-supplied `sref` takes precedence over any CRS derived from a Raster
or SpatRaster object.

## Examples

``` r
cov <- prep_cov(skjepo$cov)
```
