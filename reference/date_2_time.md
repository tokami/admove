# Convert dates to numeric time since an origin

Convert a vector of valid dates or date-times to a numeric time scale
measured since a reference origin.

## Usage

``` r
date_2_time(dates, tref = NULL)
```

## Arguments

- dates:

  A vector of class `Date`, `POSIXct`, or `POSIXlt`.

- tref:

  Optional time reference. One of:

  - `NULL`: infer origin and units from `dates`.

  - a single `Date` or `POSIXct`/`POSIXlt`: use as origin and infer
    units from `dates`.

  - a list with element `origin` and optional element `units`, e.g.
    `list(origin = as.Date("2020-01-01"), units = "day")`.

## Value

A numeric vector of the same length as `dates`. The vector has an
attribute `"tref"`, a list with elements:

- origin:

  The reference origin used.

- units:

  The time units used for the numeric scale.

- floor_unit:

  The unit used to floor the inferred origin, if relevant.

- inferred:

  Logical; whether the reference was inferred.

## Details

If `tref` is supplied, the function uses that reference. If
`tref = NULL`, the function infers a sensible reference from the input:
the origin is based on the earliest non-missing date (optionally floored
to a convenient boundary), and the time units are guessed from the
overall time range.

The returned value is a numeric vector with an attribute `"tref"`
containing the reference used. This makes it easy to reuse the same
reference later, e.g. by calling
`date_2_time(new_dates, tref = attr(x, "tref"))`.

For units `"month"` and `"year"`, the conversion is calendar-aware and
uses
[`lubridate::time_length()`](https://lubridate.tidyverse.org/reference/time_length.html)
on an interval rather than fixed-day approximations.

## Examples

``` r
d <- as.Date(c("2020-01-15", "2020-02-01", "2020-03-10"))

x <- date_2_time(d)
x
#> [1]  0 17 55
#> attr(,"tref")
#> attr(,"tref")$origin
#> [1] "2020-01-15"
#> 
#> attr(,"tref")$units
#> [1] "day"
#> 
#> attr(,"tref")$floor_unit
#> [1] "day"
#> 
#> attr(,"tref")$inferred
#> [1] TRUE
#> 
attr(x, "tref")
#> $origin
#> [1] "2020-01-15"
#> 
#> $units
#> [1] "day"
#> 
#> $floor_unit
#> [1] "day"
#> 
#> $inferred
#> [1] TRUE
#> 

## Reuse the same reference
date_2_time(d, tref = attr(x, "tref"))
#> [1]  0 17 55
#> attr(,"tref")
#> attr(,"tref")$origin
#> [1] "2020-01-15"
#> 
#> attr(,"tref")$units
#> [1] "day"
#> 
#> attr(,"tref")$floor_unit
#> [1] "day"
#> 
#> attr(,"tref")$inferred
#> [1] TRUE
#> 

## Explicit origin and units
date_2_time(d, tref = list(origin = as.Date("2020-01-01"), units = "day"))
#> [1] 14 31 69
#> attr(,"tref")
#> attr(,"tref")$origin
#> [1] "2020-01-01"
#> 
#> attr(,"tref")$units
#> [1] "day"
#> 
#> attr(,"tref")$floor_unit
#> [1] "day"
#> 
#> attr(,"tref")$inferred
#> [1] FALSE
#> 

z <- as.POSIXct(c("2020-01-01 00:00:00",
                  "2020-01-01 12:00:00",
                  "2020-01-02 06:00:00"),
                tz = "UTC")
date_2_time(z, tref = list(origin = z[1], units = "hour"))
#> [1]  0 12 30
#> attr(,"tref")
#> attr(,"tref")$origin
#> [1] "2020-01-01 UTC"
#> 
#> attr(,"tref")$units
#> [1] "hour"
#> 
#> attr(,"tref")$floor_unit
#> [1] "hour"
#> 
#> attr(,"tref")$inferred
#> [1] FALSE
#> 
```
