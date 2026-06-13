# Transform stored coordinates to a new CRS while preserving locations

Transform the stored coordinates of an object from its current spatial
reference to a new CRS, and update `sref(x)` accordingly.

## Usage

``` r
transform_sref(
  x,
  sref = NULL,
  crs = NULL,
  units = NULL,
  crs_scale = NULL,
  verbose = TRUE
)
```

## Arguments

- x:

  An object with an existing valid `sref` and supported stored
  coordinates.

- sref:

  Optional object or sref-like specification providing target
  spatial-reference fields. This can be:

  - an object of class `"admove_sref"`,

  - a named list containing one or more of `crs`, `units`, and
    `crs_scale`,

  - another object from which `sref(sref)` can be extracted.

- crs:

  Optional target CRS. Can be an
  [`sf::crs`](https://r-spatial.github.io/sf/reference/coerce-methods.html)
  object, an EPSG code, or a character string. If supplied, this takes
  precedence over `sref$crs`.

- units:

  Optional target stored coordinate units, e.g. `"m"`, `"km"`, or
  `"degree"`. If `NULL`, they are inferred from the target CRS where
  possible.

- crs_scale:

  Optional numeric scalar giving the conversion factor from target CRS
  units to target stored units. If `NULL`, it is inferred from `crs` and
  `units` where possible.

- verbose:

  Logical; if `TRUE`, informative messages are printed.

## Value

`x` with transformed coordinates and updated `sref(x)`.

## Details

This function should be used when an object already has a valid spatial
reference and you want to adopt a different CRS without changing the
underlying locations represented by the coordinates. In other words, it
changes the spatial coordinate system, not the represented points.

The source coordinates are interpreted using the current `sref(x)`. The
current `crs_scale` is used to convert stored coordinates back to CRS
units before transformation. After transformation, coordinates are
converted to the requested stored units using the target `crs_scale`.

The new spatial reference can be supplied directly via `crs`, `units`,
and `crs_scale`, or indirectly via `sref` as a sref-like object or
another object from which `sref(sref)` can be extracted.

At present, this function is implemented for point-based objects with
coordinate columns named `x`/`y`, `x0`/`y0`, `x1`/`y1`, and so on. It is
not implemented for `admove_grid` or `admove_cov`, because reprojection
of regular grids and covariate arrays requires resampling rather than
simple coordinate transformation.

Let \\c\_{\mathrm{old}}\\ denote stored coordinates,
\\s\_{\mathrm{old}}\\ the old `crs_scale`, and \\\mathcal{T}\\ the CRS
transformation from the old CRS to the new CRS. The transformation is
applied as: \$\$ c\_{\mathrm{new}} = \mathcal{T}\left(c\_{\mathrm{old}}
/ s\_{\mathrm{old}}\right) \times s\_{\mathrm{new}}, \$\$ where
coordinates are first converted from stored units to source CRS units,
transformed to the target CRS, and then converted from target CRS units
to target stored units.

If the target CRS is the same as the current CRS, the function falls
back to
[`scale_sref`](https://tokami.github.io/admove/reference/scale_sref.md)
if only the stored coordinate units differ.

## See also

[`add_sref`](https://tokami.github.io/admove/reference/add_sref.md),
[`scale_sref`](https://tokami.github.io/admove/reference/scale_sref.md),
[`create_sref`](https://tokami.github.io/admove/reference/create_sref.md),
[`sref`](https://tokami.github.io/admove/reference/sref.md),
[`.guess_crs_scale`](https://tokami.github.io/admove/reference/dot-guess_crs_scale.md)

## Examples

``` r
## Not run:
## Reproject point coordinates from lon/lat to a projected CRS
## tags <- transform_sref(tags, crs = 3035, units = "km", crs_scale = 0.001)

## Use the spatial reference of another object
## tags <- transform_sref(tags, sref = other_tags)
## End(Not run)
```
