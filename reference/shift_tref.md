# Shift the origin of a time reference while preserving represented dates

Change the origin of an object's time reference (`tref`) and shift the
stored numeric time values accordingly, so that the represented dates or
date-times remain unchanged.

## Usage

``` r
shift_tref(x, tref = NULL, origin = NULL, verbose = TRUE)
```

## Arguments

- x:

  An object with an existing valid tref and numeric time values stored
  in a supported format.

- tref:

  Optional object or tref-like specification providing the new origin.
  Ignored if `origin` is supplied. This can be:

  - an object of class `"admove_tref"`,

  - a named list containing at least `origin`,

  - another object from which `tref(tref)` can be extracted.

- origin:

  Optional new origin as a `Date` or `POSIXct`-like object. If supplied,
  this takes precedence over `tref`.

- verbose:

  Logical; if `TRUE`, informative messages are printed about the shift
  applied.

## Value

The input object `x`, with numeric time values shifted and
`tref(x)$origin` replaced by the new origin.

## Details

This function should be used when an object already has a valid tref and
you want to adopt a different origin without changing the underlying
time points that the numeric values represent. In other words, it
changes the coordinate system of the time axis, not the actual dates.

The new origin can be supplied directly via `origin`, or indirectly via
`tref` as a tref-like object or another object from which `tref(tref)`
can be extracted.

Let \\t\\ denote the stored numeric time values, \\o_0\\ the old origin,
and \\o_1\\ the new origin. To preserve the represented dates, the new
time values \\t'\\ must satisfy: \$\$o_0 + t = o_1 + t'\$\$ and
therefore \$\$t' = t - (o_1 - o_0).\$\$

The offset \\o_1 - o_0\\ is computed in the current tref units. For
`"second"`, `"minute"`, `"hour"`, `"day"`, and `"week"`, fixed-duration
differences are used. For `"month"` and `"year"`, the offset is computed
in a calendar-aware way using
[`lubridate::time_length()`](https://lubridate.tidyverse.org/reference/time_length.html)
on the interval between the two origins.

This function requires that `x` already has a valid tref with
non-missing `origin` and `units`. It does not change tref units or
period.

## See also

[`add_tref`](https://tokami.github.io/admove/reference/add_tref.md),
[`scale_tref`](https://tokami.github.io/admove/reference/scale_tref.md),
[`create_tref`](https://tokami.github.io/admove/reference/create_tref.md),
[`tref`](https://tokami.github.io/admove/reference/tref.md)

## Examples

``` r
## Shift an object to a new explicit origin
x <- shift_tref(skjepo$sim$dat, origin = as.Date("1995-01-01"))
#> Shifted time values by -300 month to match new tref origin.

## Shift x to use the same origin as y
## x <- shift_tref(x, tref = y)

## Equivalent workflow via add_tref()
## x <- shift_tref(x, tref = tref(skjepo$sim$dat))
```
