# Prepare tagging data for admove

`prep_tags()` converts raw tagging data into a standardised object of
class `admove_tags` used throughout *admove*. It supports data-storage
(archival) tags, mark-resight tags, and mark-recapture tags, and accepts
input either in long format, wide format, or as a list split by tag.

The function harmonises variable names, optionally converts date
variables to numeric model time, adds missing tag identifiers when
needed, and attaches spatial and temporal reference information.

## Usage

``` r
prep_tags(
  x,
  tag_type = NULL,
  names = NULL,
  date_decimal = FALSE,
  date_format = NULL,
  date_origin = NULL,
  keep_only_recaptured = TRUE,
  tz = "UTC",
  sref = NULL,
  tref = NULL,
  transform_sref = FALSE,
  shift_tref = FALSE,
  verbose = TRUE
)

prep_dtags(
  x,
  names = NULL,
  date_decimal = FALSE,
  date_format = NULL,
  date_origin = NULL,
  tz = "UTC",
  sref = NULL,
  tref = NULL,
  transform_sref = FALSE,
  shift_tref = FALSE,
  verbose = TRUE
)

prep_stags(
  x,
  names = NULL,
  date_decimal = FALSE,
  date_format = NULL,
  date_origin = NULL,
  tz = "UTC",
  sref = NULL,
  tref = NULL,
  transform_sref = FALSE,
  shift_tref = FALSE,
  verbose = TRUE
)

prep_ctags(
  x,
  names = NULL,
  date_decimal = FALSE,
  date_format = NULL,
  date_origin = NULL,
  tz = "UTC",
  sref = NULL,
  tref = NULL,
  transform_sref = FALSE,
  shift_tref = FALSE,
  verbose = TRUE
)
```

## Arguments

- x:

  Tagging data. Can be a data frame, a list of data frames (typically
  one per tag), or an existing object of class `admove_tags`.

- tag_type:

  Character string specifying the tag type: `"d"` for data-storage
  (archival) tags, `"s"` for mark-resight tags, and `"c"` for
  mark-recapture tags.

- names:

  Named character vector giving the column names in `x`. For long
  format, provide at least
  `c(t = "...", x = "...", y = "...", id = "...")`. For wide format,
  provide
  `c(t0 = "...", t1 = "...", x0 = "...", y0 = "...", x1 = "...", y1 = "...")`.
  The order does not matter, but the vector must be named.

- date_decimal:

  Logical; if `TRUE`, interpret the time variable as a decimal year and
  convert it to model time. Default: `FALSE`.

- date_format:

  Optional character string passed to
  [`as.Date()`](https://rdrr.io/r/base/as.Date.html) to parse character
  dates. Default: `NULL`.

- date_origin:

  Optional origin passed to
  [`as.Date()`](https://rdrr.io/r/base/as.Date.html) when times are
  stored numerically. Default: `NULL`.

- keep_only_recaptured:

  Logical; if `TRUE`, only keep tags with at least two observations.
  Currently mainly relevant for mark-recapture-style data. Default:
  `TRUE`.

- tz:

  Time zone used when converting dates. Default: `"UTC"`.

- sref:

  Optional spatial reference information to attach to the returned
  object.

- tref:

  Optional time reference information to attach to the returned object.

- transform_sref:

  Logical; if `TRUE`, transform coordinates to `sref` when possible.
  Default: `FALSE`.

- shift_tref:

  Logical; if `TRUE`, shift times to `tref` when possible. Default:
  `FALSE`.

- verbose:

  Logical; if `TRUE`, print informative messages. Default: `TRUE`.

## Value

A data frame of class `admove_tags` with standardised columns such as
`t`, `x`, `y`, `id`, `tag_type`, and `use`, plus any additional columns
provided in the input.

## Details

`prep_tags()` accepts three common input structures:

- a long-format data frame with one row per observation,

- a wide-format data frame with release and recapture columns (`t0`,
  `t1`, `x0`, `y0`, `x1`, `y1`), or

- a list of data frames, usually one element per tag.

If no tag identifier is supplied for list input, the list order is used
to create an `id` column automatically.

## Examples

``` r
## prepare data-storage tags
dtags <- prep_tags(
  skjepo$dtags,
  tag_type = "d",
  names = c(t = "time", x = "mptlon", y = "mptlat"),
  date_origin = "1899-12-30"
)
#> ID not specified, using order of list elements.
#> tref (time origin and units) was inferred from dates. Please check and adjust if needed.

## prepare mark-recapture tags
ctags <- prep_tags(
  skjepo$ctags,
  tag_type = "c",
  names = c(
    t0 = "date_time", t1 = "date_caught",
    x0 = "rel_lon",   x1 = "recap_lon",
    y0 = "rel_lat",   y1 = "recap_lat"
  ),
  date_origin = "1899-12-30"
)
#> tref (time origin and units) was inferred from dates. Please check and adjust if needed.
```
