# Add or harmonise a spatial reference on an object

Attach a spatial reference (`sref`) to an object, or harmonize an
existing spatial reference with that of another object or sref-like
specification.

## Usage

``` r
add_sref(x, sref = NULL, verbose = TRUE, transform_crs = FALSE)
```

## Arguments

- x:

  An object to which spatial reference information should be attached or
  whose spatial reference should be harmonized.

- sref:

  Optional spatial reference specification. This can be:

  - an object of class `"admove_sref"`,

  - a named list containing one or more of `crs`, `units`, and
    `crs_scale`,

  - another object from which `sref(sref)` can be extracted.

  If `NULL`, existing spatial reference information on `x` is kept where
  possible.

- verbose:

  Logical; if `TRUE`, informative messages are printed when CRS, units,
  or scaling information are inferred, retained, or transformed.

- transform_crs:

  Logical; if `TRUE` and the current CRS of `x` differs from the
  requested CRS, stored coordinates are transformed with
  [`transform_sref`](https://tokami.github.io/admove/reference/transform_sref.md)
  so that represented locations remain the same. If `FALSE` (default),
  differing non-missing CRS values cause an error.

## Value

`x` with updated `sref(x)`. Depending on the input and requested spatial
reference, the stored coordinates may also be rescaled (same CRS,
different `crs_scale`) or transformed (different CRS and
`transform_crs = TRUE`).

## Details

The spatial reference bundles a coordinate reference system (CRS), the
units of the stored coordinates, and a unit scaling factor linking CRS
units to the stored coordinate units.

This function is primarily intended to add or update *metadata*. If the
object already has a spatial reference and the requested CRS is the same
but `crs_scale` differs, the stored coordinates are rescaled via
[`scale_sref`](https://tokami.github.io/admove/reference/scale_sref.md).
If the requested CRS differs, changing the CRS metadata alone would
change the meaning of the stored coordinates. Therefore, by default the
function throws an error in this case. To preserve the represented
locations while adopting the new CRS, set `transform_crs = TRUE`; this
calls
[`transform_sref`](https://tokami.github.io/admove/reference/transform_sref.md)
internally.

If a CRS is provided and the sf package is available, the CRS is
normalized using
[`st_crs`](https://r-spatial.github.io/sf/reference/st_crs.html) and
stored as WKT for stability. If `units` are not provided, they are
inferred from the target CRS when possible. If `crs_scale` is not
provided, it is inferred with
[`.guess_crs_scale`](https://tokami.github.io/admove/reference/dot-guess_crs_scale.md)
where possible.

The function distinguishes between three kinds of spatial-reference
updates:

1.  **Adding missing spatial metadata**: if `x` has no spatial
    reference, or if parts of it are missing, missing fields are filled
    from `sref` where possible.

2.  **Changing stored coordinate units**: if the CRS is unchanged but
    `crs_scale` differs, stored coordinates are rescaled using
    [`scale_sref`](https://tokami.github.io/admove/reference/scale_sref.md).

3.  **Changing CRS**: if CRS differ, stored coordinates must be
    transformed if the represented locations are to remain unchanged.
    This is not done silently; use `transform_crs = TRUE` or call
    [`transform_sref`](https://tokami.github.io/admove/reference/transform_sref.md)
    directly.

## See also

[`transform_sref`](https://tokami.github.io/admove/reference/transform_sref.md),
[`scale_sref`](https://tokami.github.io/admove/reference/scale_sref.md),
[`create_sref`](https://tokami.github.io/admove/reference/create_sref.md),
[`sref`](https://tokami.github.io/admove/reference/sref.md),
[`.guess_crs_scale`](https://tokami.github.io/admove/reference/dot-guess_crs_scale.md)

## Examples

``` r
## Not run:
## Add a new spatial reference
grid <- create_grid(cellsize = 5e3, xrange = c(0, 5e4), yrange = c(0, 5e4))
grid <- add_sref(grid, sref = list(crs = 32631, units = "m"))

## Same CRS, but store coordinates in km rather than m
grid <- add_sref(grid, sref = list(crs = 32631, units = "km",
                                     crs_scale = 0.001))

## Adopt the CRS of another object while preserving locations
## tags <- add_sref(tags, sref = other_tags, transform_crs = TRUE)
## End(Not run)
```
