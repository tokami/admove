# Compare two time references for equality

Test whether two `admove_tref` objects describe the same time reference.
By default, equality requires matching `units` and (when defined)
matching `period`. Optionally, the `origin` can be checked as well.

## Usage

``` r
tref_equal(a, b, tol = 1e-12, check_origin = FALSE)
```

## Arguments

- a, b:

  Objects to compare.

- tol:

  Numeric tolerance used for comparing `period` and (if
  `check_origin = TRUE`) `origin`. For `origin`, the tolerance is
  interpreted in seconds.

- check_origin:

  Logical; if `TRUE`, also compare `origin`.

## Value

Logical scalar; `TRUE` if equal, otherwise `FALSE`.

## Details

The comparison is designed to be robust to missing values:

- `units` must match exactly (character comparison after normalization).

- `period` is considered equal if both are `NA`, or if both are finite
  and their absolute difference is `<= tol`.

- If `check_origin = TRUE`, `origin` is considered equal if both are
  `NA`, or if both are defined and their absolute time difference is
  `<= tol` seconds.

## Examples

``` r
tr1 <- create_tref(origin = as.Date("2025-01-01"), units = "month")
tr2 <- create_tref(origin = as.Date("2025-01-01"), units = "month", period = 12)
tref_equal(tr1, tr2)
#> [1] TRUE

tr3 <- create_tref(origin = as.Date("2025-01-02"), units = "month")
tref_equal(tr1, tr3, check_origin = TRUE)
#> [1] FALSE
```
