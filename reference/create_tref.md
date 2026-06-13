# Time reference for admove objects

A lightweight container that stores the time reference information used
throughout admove. It bundles:

- **origin**: a POSIXct date-time used as the temporal origin

- **units**: the units of the numeric time axis stored in objects

- **period**: optional seasonal cycle length in the same units as the
  stored times

## Usage

``` r
create_tref(origin = NA, units = NA_character_, period = NULL)
```

## Arguments

- origin:

  A date-time corresponding to the temporal origin (e.g. a release
  date). Will be converted to `POSIXct`. Use `NA` to leave undefined.

- units:

  Character string describing the units of the stored numeric time axis,
  e.g. `"year"`, `"quarter"`, `"month"`, `"week"`, `"day"`. Use
  `list_units_time()` for supported options.

- period:

  Optional numeric seasonal period (cycle length) in the same units as
  the stored time axis. If `NULL`, a default is inferred for common
  annual discretizations when `units` is known.

## Value

An object of class `admove_tref`.

## Details

The intention is that *objects store their times on a simple numeric
axis* (e.g. months since `origin`), while `origin` and `units` define
how these values should be interpreted and displayed.

If `period` is provided (or can be inferred), it can be used to create
seasonally repeating covariate fields or spline bases by wrapping time
with modulo arithmetic. For example, if times are stored in months, a
natural annual seasonal cycle is `period = 12`. For custom
discretizations (e.g. 10 time steps per year), set `period = 10` and
`units = "custom"`.

If `units` is missing but `period` is supplied, `units` can be inferred
for common annual discretizations (`period = 1, 2, 4, 12, 52`) and
otherwise defaults to `"custom"`.

## Examples

``` r
## Monthly time axis with annual seasonality (12 months)
tr <- create_tref(origin = as.Date("2025-01-01"), units = "month")
tr
#> $origin
#> [1] "2025-01-01 UTC"
#> 
#> $units
#> [1] "month"
#> 
#> $period
#> [1] 12
#> 
#> attr(,"class")
#> [1] "admove_tref"

## Infer units from period
tr2 <- create_tref(origin = as.Date("2025-01-01"), period = 12)
tr2
#> $origin
#> [1] "2025-01-01 UTC"
#> 
#> $units
#> [1] "month"
#> 
#> $period
#> [1] 12
#> 
#> attr(,"class")
#> [1] "admove_tref"

## Custom discretization: 10 time steps per year
tr10 <- create_tref(origin = as.Date("2025-01-01"), period = 10)
tr10
#> $origin
#> [1] "2025-01-01 UTC"
#> 
#> $units
#> [1] "custom"
#> 
#> $period
#> [1] 10
#> 
#> attr(,"class")
#> [1] "admove_tref"
```
