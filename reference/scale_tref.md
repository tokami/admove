# Change time units by scaling an admove object

Rescale the stored time axis of an `admove_*` object and update its time
reference (`tref(x)`). This is a purely numeric scaling of the stored
times and is appropriate when the time axis is represented in
fixed-length model units (e.g. "month" meaning 1/12 year, "quarters"
meaning 1/4 year, or a custom discretization with `tref$period` steps
per year).

## Usage

``` r
scale_tref(x, scale = 1, units = NULL, verbose = TRUE)
```

## Arguments

- x:

  An `admove_tags`, `admove_cov`, `admove_data` object, or a numeric
  vector/matrix/array of times.

- scale:

  Numeric scalar. Stored times are multiplied by this factor. For
  example, to convert years to months (numerically), use `scale = 12`.

- units:

  Optional character string giving the new units label to store in
  `tref(x)$units`. If provided, `scale` is inferred from
  `guess_t_crs_scale(tref(x), units)`.

- verbose:

  Logical; if `TRUE`, prints informational messages.

## Value

`x` with rescaled stored times and updated `tref(x)`.

## Details

The time reference is updated by:

- multiplying stored times by `scale`

- multiplying `tref(x)$period` by `scale` (if defined)

- updating `tref(x)$units` (either from `units` or by guessing)

## See also

[`create_tref`](https://tokami.github.io/admove/reference/create_tref.md),
[`add_tref`](https://tokami.github.io/admove/reference/add_tref.md),
[`tref`](https://tokami.github.io/admove/reference/tref.md)

## Examples

``` r
## Scale an object to a new time step using scale directly (e.g. from month to week)
x <- scale_tref(skjepo$sim$dat, scale = 1/12*52)

## Scale by providing units directly
x <- scale_tref(skjepo$sim$tags, units = "year")
```
