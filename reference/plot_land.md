# Plot land masses in the current plotting region

Add land polygons to an existing plot, using the spatial reference
stored in `sref`. Land is transformed to the requested coordinate
reference system, optionally rescaled to match the plotting units,
cropped to the current plot extent, and then added to the active
graphics device.

## Usage

``` r
plot_land(
  sref = NULL,
  col = grDevices::adjustcolor(grey(0.7), 0.5),
  border = grey(0.5),
  download_map = FALSE,
  scale = 110,
  verbose = TRUE,
  warn_once = TRUE
)
```

## Arguments

- sref:

  Optional spatial reference object, typically as returned by
  [`sref()`](https://tokami.github.io/admove/reference/sref.md). It
  should contain at least a valid CRS in `sref$crs`, and may also
  contain plotting units in `sref$units` and a scaling factor in
  `sref$crs_scale`.

- col:

  Fill colour for land polygons. Default is
  `grDevices::adjustcolor(grey(0.7), 0.5)`.

- border:

  Border colour for land polygons. Default is `grey(0.5)`.

- download_map:

  Logical; if `TRUE`, the land map is downloaded if needed. Otherwise, a
  locally available map is used when possible. Default is `FALSE`.

- scale:

  Numeric map scale passed to the internal land-data loader. Default is
  `110`.

- verbose:

  Logical; if `TRUE`, informative messages are printed. Default is
  `TRUE`.

- warn_once:

  Logical; if `TRUE`, warnings about missing CRS information are printed
  only once per session. Default is `TRUE`.

## Value

Invisibly returns `NULL`. Called for its side effect of adding land
masses to an existing plot.

## Details

The function requires package sf. If `sref` does not contain a valid
CRS, no land is plotted. The land polygons are transformed to the
requested CRS, optionally multiplied by `sref$crs_scale`, and cropped to
the current plotting region defined by `par("usr")`.

For geographic coordinate systems, the plotting extent is truncated to
valid longitude and latitude ranges before cropping.
