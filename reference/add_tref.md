# Add or harmonize a time reference on an object

Attach a time reference (`tref`) to an object, or harmonize an existing
time reference with that of another object or tref-like specification.

## Usage

``` r
add_tref(x, tref = NULL, verbose = TRUE, shift_origin = FALSE)
```

## Arguments

- x:

  An object to which a tref should be added or whose tref should be
  harmonized.

- tref:

  Optional tref specification. This can be:

  - an object of class `"admove_tref"`,

  - a named list containing one or more of `origin`, `units`, and
    `period`,

  - another object from which `tref(tref)` can be extracted.

  If `NULL`, the existing tref on `x` is kept where possible.

- verbose:

  Logical; if `TRUE`, informative messages are printed when tref fields
  are inferred, retained, or converted.

- shift_origin:

  Logical; if `TRUE` and the origin of `x` differs from the requested
  origin, the stored time values are shifted using
  [`shift_tref`](https://tokami.github.io/admove/reference/shift_tref.md)
  so that the represented dates remain unchanged. If `FALSE` (default),
  differing non-missing origins cause an error.

## Value

The input object `x`, with updated tref metadata. Depending on the input
and requested tref, the stored numeric time values may also be rescaled
(if units differ) or shifted (if `shift_origin = TRUE`).

## Details

This function is intended for adding or updating *metadata* describing a
numeric time axis, such as its origin, units, and period. If the object
already has a tref and the requested tref has different time units, the
stored time values are rescaled via
[`scale_tref`](https://tokami.github.io/admove/reference/scale_tref.md)
where possible.

If the object already has a tref and the requested origin differs,
changing the origin without changing the stored time values would change
the meaning of the time series. Therefore, by default, the function
throws an error in this case. To preserve the represented dates while
adopting the new origin, set `shift_origin = TRUE`; this calls
[`shift_tref`](https://tokami.github.io/admove/reference/shift_tref.md)
internally to shift the stored time values before updating the tref
metadata.

The argument `tref` may be a tref-like object, or any object for which
`tref(tref)` can be extracted.

The function distinguishes between three different kinds of tref
changes:

1.  **Adding missing tref metadata**: if `x` has no tref, or if parts of
    its tref are missing, the missing fields are filled from `tref`
    where available.

2.  **Changing units**: if time units differ between `x` and `tref`, the
    stored numeric time values are converted using
    [`scale_tref`](https://tokami.github.io/admove/reference/scale_tref.md)
    where possible.

3.  **Changing origin**: if origins differ, the stored time values must
    also change if the represented dates are to remain the same. This is
    not done silently; use `shift_origin = TRUE` or call
    [`shift_tref`](https://tokami.github.io/admove/reference/shift_tref.md)
    directly.

Period is treated strictly: if both `x` and `tref` define a finite
period and the two differ, an error is raised.

## See also

[`shift_tref`](https://tokami.github.io/admove/reference/shift_tref.md),
[`scale_tref`](https://tokami.github.io/admove/reference/scale_tref.md),
[`create_tref`](https://tokami.github.io/admove/reference/create_tref.md),
[`tref`](https://tokami.github.io/admove/reference/tref.md)

## Examples

``` r
## Assume x and y are objects with numeric time and tref metadata
## x <- add_tref(x, list(origin = as.Date("2020-01-01"), units = "day"))

## Copy tref fields from another object
## x <- add_tref(x, y)

## Convert units if needed
## x <- add_tref(x, list(units = "month"))

## Adopt the origin of y while preserving represented dates
## x <- add_tref(x, y, shift_origin = TRUE)
```
