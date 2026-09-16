
##' @importFrom grDevices adjustcolor col2rgb grey hcl.colors n2mfrow rgb terrain.colors
##' @importFrom graphics abline arrows axis box contour identify image layout legend lines mtext par plot.new points polygon rect segments text title grconvertX grconvertY
##' @importFrom stats approx dist median qnorm quantile rnorm runif setNames
##' @importFrom utils capture.output head packageDescription tail
##' @importFrom RTMB ADoverload REPORT ADREPORT
##' @importFrom lubridate interval time_length floor_date
##' @keywords internal
"_PACKAGE"

utils::globalVariables(c("get_cov", "get_sim_par", "get_sim_funcs"))

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

.is_na_scalar <- function(x) {
  length(x) == 1L && is.na(x)
}


.add_class <- function(x, class) {
  if (is.null(x)) return(x)
  if (!inherits(x, class)) {
    class(x) <- c(class, class(x))
  }
  return(x)
}


.check_class <- function(x, class) {
  if (!inherits(x, class)) {
    stop(paste0("The object ", deparse(substitute(x)),
                " does not inherit class ", class,
                ". Please check your code."))
  }
  return(invisible(NULL))
}


.get_land <- function(download_map = FALSE, scale = 110, make_valid = TRUE) {

  if (!download_map) {
    path <- system.file("extdata", "land_110m.rds", package = "admove")
    if (nzchar(path)) return(readRDS(path))
  }

  if (!requireNamespace("rnaturalearth", quietly = TRUE)) {
    stop("To download land polygons, install 'rnaturalearth'.")
  }

  land <- rnaturalearth::ne_download(
    scale = scale, type = "land", category = "physical", returnclass = "sf"
  )

  if (isTRUE(make_valid)) {
    if (!requireNamespace("sf", quietly = TRUE)) {
      stop("To make the spatial object valid, install 'sf'.")
    }
    land <- sf::st_make_valid(land)
  }

  land
}


.admove_cols <- function(n = 1, alpha = 1, type = NULL){

  if(is.null(type) || is.na(type)){
    adjustcolor(c("dodgerblue3","goldenrod2",
                  "darkgreen","purple",
                  "hotpink")[1:n], alpha)
  }else if(type == "true"){
    adjustcolor("darkorange", alpha)
  }else if(type == "pos"){
    adjustcolor("dodgerblue3", alpha)
  }else if(type == "neg"){
    adjustcolor("firebrick3", alpha)
  }else if(type == "notsig"){
    adjustcolor("chartreuse4", alpha)
  }else if(type == "sig"){
    adjustcolor("firebrick3", alpha)
  }else{
    adjustcolor(c("dodgerblue3","goldenrod2",
                  "darkgreen","purple",
                  "hotpink")[1:n], alpha)
  }
}

.is_empty <- function(x){
  is.null(x) || length(x) == 0
}


.pfx <- function(x, eps=0.01){
  return(eps * log(exp(x/eps) + 1))
}

.smooth_identity <- function(x, from=0, to=1){
  z <- (x - from) / (to - from)
  y <- 1 - .pfx(1 - .pfx(z))
  return(y * (to - from) + from)
}


.get_adv <- function(dat, par, conf, funcs = NULL){
  par <- get_sim_par(par)
  cov <- .make_cov_list(dat$cov)
  dat$cov <- cov
  dat$pred$cov <- NULL
  conf <- default_conf(dat)
  funcs <- get_sim_funcs(funcs, dat, conf, cov, par)
  hTdx.true <- sapply(dat$pred$time,
                      function(t) apply(dat$pred$xygrid, 1, function(x) funcs$tax(x,t)[1]))
  hTdy.true <- sapply(dat$pred$time,
                      function(t) apply(dat$pred$xygrid, 1, function(x) funcs$tax(x,t)[2]))
  uv.true <- matrix(NA, nrow(dat$xygrid), 2)
  uv.true[,1] <- rowMeans(hTdx.true)
  uv.true[,2] <- rowMeans(hTdy.true)

  return(uv.true)
}

.get_par_names <- function(fit){
  tab <- table(names(fit$opt$par))
  res <- unlist(sapply(seq_along(tab), function(x) if(tab[x] > 1) paste0(names(tab)[x],1:tab[x]) else names(tab)[x]))
  return(res)
}

.date_2_decimal_year <- function(dates) {
  leap_year <- function(year) {
    (year %% 4 == 0 & year %% 100 != 0) | (year %% 400 == 0)
  }
  year <- as.numeric(format(dates, "%Y"))
  doy <- as.numeric(format(dates, "%j")) +
    as.numeric(format(dates, "%H")) / 24 +
    as.numeric(format(dates, "%M")) / 60 / 24 +
    as.numeric(format(dates, "%S")) / 60 / 60 / 24
  is_leap <- leap_year(year)
  days_in_year <- ifelse(is_leap, 366, 365)
  dec_year <- year + (doy - 1) / days_in_year
  return(dec_year)
}

.decimal_year_2_date <- function(dec_year, tz = "UTC") {
  leap_year <- function(year) {
    (year %% 4 == 0 & year %% 100 != 0) | (year %% 400 == 0)
  }

  year <- floor(dec_year)
  frac <- dec_year - year

  is_leap <- leap_year(year)
  days_in_year <- ifelse(is_leap, 366, 365)

  start_of_year <- as.POSIXct(
    sprintf("%04d-01-01 00:00:00", year),
    tz = tz
  )

  date <- start_of_year + frac * days_in_year * 24 * 60 * 60
  return(date)
}


#' Convert dates to numeric time since an origin
#'
#' Convert a vector of valid dates or date-times to a numeric time scale measured
#' since a reference origin.
#'
#' If \code{tref} is supplied, the function uses that reference. If
#' \code{tref = NULL}, the function infers a sensible reference from the input:
#' the origin is based on the earliest non-missing date (optionally floored to a
#' convenient boundary), and the time units are guessed from the overall time
#' range.
#'
#' The returned value is a numeric vector with an attribute \code{"tref"}
#' containing the reference used. This makes it easy to reuse the same reference
#' later, e.g. by calling \code{date_2_time(new_dates, tref = attr(x, "tref"))}.
#'
#' For units \code{"month"} and \code{"year"}, the conversion is calendar-aware
#' and uses \code{lubridate::time_length()} on an interval rather than fixed-day
#' approximations.
#'
#' @param dates A vector of class \code{Date}, \code{POSIXct}, or
#'   \code{POSIXlt}.
#' @param tref Optional time reference. One of:
#'   \itemize{
#'   \item \code{NULL}: infer origin and units from \code{dates}.
#'   \item a single \code{Date} or \code{POSIXct}/\code{POSIXlt}: use as origin
#'     and infer units from \code{dates}.
#'   \item a list with element \code{origin} and optional element \code{units},
#'     e.g. \code{list(origin = as.Date("2020-01-01"), units = "day")}.
#'   }
#'
#' @return A numeric vector of the same length as \code{dates}. The vector has an
#'   attribute \code{"tref"}, a list with elements:
#'   \describe{
#'   \item{origin}{The reference origin used.}
#'   \item{units}{The time units used for the numeric scale.}
#'   \item{floor_unit}{The unit used to floor the inferred origin, if relevant.}
#'   \item{inferred}{Logical; whether the reference was inferred.}
#'   }
#'
#' @examples
#' d <- as.Date(c("2020-01-15", "2020-02-01", "2020-03-10"))
#'
#' x <- date_2_time(d)
#' x
#' attr(x, "tref")
#'
#' ## Reuse the same reference
#' date_2_time(d, tref = attr(x, "tref"))
#'
#' ## Explicit origin and units
#' date_2_time(d, tref = list(origin = as.Date("2020-01-01"), units = "day"))
#'
#' z <- as.POSIXct(c("2020-01-01 00:00:00",
#'                   "2020-01-01 12:00:00",
#'                   "2020-01-02 06:00:00"),
#'                 tz = "UTC")
#' date_2_time(z, tref = list(origin = z[1], units = "hour"))
#'
#' @export
date_2_time <- function(dates, tref = NULL) {

  `%||%` <- function(x, y) {
    if (is.null(x)) y else x
  }

  ## check input class
  if (!inherits(dates, "Date") && !inherits(dates, "POSIXt")) {
    stop("'dates' must be of class 'Date' or 'POSIXt'.")
  }

  ## time zone of 'dates' (a Date is a calendar day, i.e. UTC by convention)
  tz <- "UTC"
  if (inherits(dates, "POSIXt")) {
    tzd <- attr(dates, "tzone", exact = TRUE)
    if (!is.null(tzd) && length(tzd) && nzchar(tzd[1L])) tz <- tzd[1L]
  }

  ## helper to coerce to same broad class as 'dates'
  if (inherits(dates, "POSIXt")) {
    coerce_time <- function(x) {
      as.POSIXct(x, tz = tz)
    }
  } else {
    coerce_time <- function(x) {
      as.Date(x)
    }
  }

  ## the origin may carry a zone of its own (e.g. "2003-01-01 UTC" or a POSIXct
  ## in another zone), which as.POSIXct()/as.Date() would silently drop
  coerce_origin <- function(x) {
    x <- .origin_2_posix(x, tz = tz)
    if (inherits(dates, "POSIXt")) x else as.Date(x, tz = attr(x, "tzone"))
  }

  dates <- coerce_time(dates)

  ## helper to guess units from overall range
  guess_units <- function(x) {

    x <- x[!is.na(x)]
    if (!length(x)) return("day")
    if (length(x) == 1L) return("day")

    span_sec <- as.numeric(max(x) - min(x), units = "secs")
    span_day <- span_sec / 86400

    if (inherits(x, "Date")) {
      if (span_day <= 180) {
        return("day")
      } else if (span_day <= 365.25 * 3) {
        return("week")
      } else if (span_day <= 365.25 * 15) {
        return("month")
      } else {
        return("year")
      }
    }

    if (span_sec <= 120) {
      "second"
    } else if (span_sec <= 7200) {
      "minute"
    } else if (span_sec <= 86400 * 3) {
      "hour"
    } else if (span_day <= 180) {
      "day"
    } else if (span_day <= 365.25 * 3) {
      "week"
    } else if (span_day <= 365.25 * 15) {
      "month"
    } else {
      "year"
    }
  }

  ## helper to choose flooring unit for inferred origin
  guess_floor_unit <- function(units) {
    units <- .normalise_time_unit(units)
    if (is.na(units)) return("day")
    switch(units,
           second = "second",
           minute = "minute",
           hour = "hour",
           day = "day",
           week = "week",
           month = "month",
           year = "year",
           "day")
  }

  units <- NULL
  origin <- NULL
  inferred <- FALSE
  floor_unit <- NULL

  ## parse tref
  if (is.null(tref)) {

    inferred <- TRUE
    units <- guess_units(dates)
    floor_unit <- guess_floor_unit(units)

    dmin <- min(dates, na.rm = TRUE)
    origin <- lubridate::floor_date(dmin, unit = floor_unit)
    origin <- coerce_time(origin)

  } else if (inherits(tref, "Date") || inherits(tref, "POSIXt")) {

    inferred <- FALSE
    units <- guess_units(dates)
    origin <- coerce_origin(tref)

  } else if (is.list(tref)) {

    inferred <- isTRUE(tref$inferred)
    units <- tref$units %||% guess_units(dates)
    floor_unit <- tref$floor_unit %||% guess_floor_unit(units)

    origin <- tref$origin %||% tref$t0 %||% tref$start
    if (is.null(origin)) {
      stop("If 'tref' is a list, it must contain an element named 'origin'.")
    }
    origin <- coerce_origin(origin)

  } else {
    stop("'tref' must be NULL, a Date/POSIXt object, or a list with element 'origin'.")
  }

  ## normalise unit aliases (e.g. "months" -> "month", "secs" -> "second")
  units_norm <- .normalise_time_unit(units)
  if (!is.na(units_norm)) units <- units_norm

  ## compute numeric time since origin
  out <- rep(NA_real_, length(dates))
  ok <- !is.na(dates)

  if (units %in% c("second", "minute", "hour", "day", "week")) {

    difftime_units <- switch(units,
                             second = "secs",
                             minute = "mins",
                             hour = "hours",
                             day = "days",
                             week = "weeks")

    out[ok] <- as.numeric(difftime(dates[ok], origin, units = difftime_units))

  } else if (units %in% c("month", "year")) {

    out[ok] <- lubridate::time_length(
      lubridate::interval(origin, dates[ok]),
      unit = units
    )

  } else {
    stop("Unsupported units: ", units, ". Supported units are 'second', ",
         "'minute', 'hour', 'day', 'week', 'month' and 'year'.")
  }

  attr(out, "tref") <- list(
    origin = origin,
    units = units,
    floor_unit = floor_unit,
    inferred = inferred
  )

  out
}


## -- time_2_date ---------------------------------------------------------------

## Resolve a flexible time-reference argument to a list(origin, units).
## Accepts an admove_tref, an admove object (its tref() is used), a list with
## 'origin'/'units', or a bare Date/POSIXt origin (units default to "day").
.resolve_tref <- function(ref) {
  if (is.null(ref)) {
    stop("'tref' must be provided to convert numeric time to dates.")
  }
  if (inherits(ref, c("admove_tags", "admove_cov", "admove_data",
                      "admove_grid", "admove"))) {
    ref <- tref(ref)
  }
  if (inherits(ref, "admove_tref")) {
    return(list(origin = ref$origin, units = ref$units))
  }
  if (inherits(ref, "Date") || inherits(ref, "POSIXt")) {
    return(list(origin = ref, units = "day"))
  }
  if (is.list(ref)) {
    origin <- ref$origin
    if (is.null(origin)) origin <- ref$t0
    if (is.null(origin)) origin <- ref$start
    if (is.null(origin)) {
      stop("If 'tref' is a list, it must contain an element named 'origin'.")
    }
    units <- ref$units
    if (is.null(units)) {
      stop("If 'tref' is a list, it must contain an element named 'units'.")
    }
    return(list(origin = origin, units = units))
  }
  stop("'tref' must be an admove_tref, an admove object, a list with ",
       "'origin'/'units', or a Date/POSIXt origin.")
}

## Invert a calendar-aware (month/year) numeric time to dates.
## date_2_time() defines these units via lubridate::time_length() on an
## interval, whose fractional part depends on the (variable) length of the
## containing month/year and can even interact with DST. Rather than reproduce
## that piecewise, DST-aware map in closed form, we invert date_2_time()
## itself by bisection: it is strictly increasing in date, so this recovers the
## unique date whenever date_2_time() is injective and is exact to machine
## precision in the round-trip sense date_2_time(time_2_date(t)) == t.
.calendar_time_2_date <- function(t, origin, unit) {
  tz <- attr(origin, "tzone")
  if (is.null(tz) || !nzchar(tz)) tz <- "UTC"

  fwd <- function(sec) {
    as.numeric(date_2_time(.POSIXct(sec, tz = tz),
                           tref = list(origin = origin, units = unit)))
  }

  sec_unit <- if (unit == "year") 365.25 * 86400 else 365.25 / 12 * 86400
  os <- as.numeric(origin)
  lo <- os + (t - 2) * sec_unit
  hi <- os + (t + 2) * sec_unit

  ## widen the bracket in the (rare) event the root falls outside it
  for (k in seq_len(64)) {
    b <- which(fwd(lo) > t)
    if (!length(b)) break
    lo[b] <- lo[b] - sec_unit
  }
  for (k in seq_len(64)) {
    b <- which(fwd(hi) < t)
    if (!length(b)) break
    hi[b] <- hi[b] + sec_unit
  }

  ## bisection (date_2_time is monotone increasing in date)
  for (i in seq_len(60)) {
    mid <- (lo + hi) / 2
    higher <- fwd(mid) < t
    lo[higher] <- mid[higher]
    hi[!higher] <- mid[!higher]
  }

  .POSIXct((lo + hi) / 2, tz = tz)
}

## Core converter: numeric time (+ origin/units) -> Date/POSIXct. Mirrors the
## unit handling of date_2_time() exactly (fixed durations for sub-monthly
## units, calendar-aware for month/year), so the two are inverses.
.time_2_date_core <- function(t, origin, units) {
  u <- .normalise_time_unit(units)
  if (is.na(u) || u == "custom") stop("Unsupported time unit: ", units)

  is_date <- inherits(origin, "Date") && !inherits(origin, "POSIXt")

  ## honour a zone carried by the origin itself; zone-less origins are UTC
  origin_ct <- .origin_2_posix(origin)
  tz <- attr(origin_ct, "tzone", exact = TRUE) %||% "UTC"
  if (length(tz) != 1L || is.na(tz) || !nzchar(tz)) tz <- "UTC"

  out <- .POSIXct(rep(NA_real_, length(t)), tz = tz)
  ok <- is.finite(t)

  if (any(ok)) {
    if (u %in% c("second", "minute", "hour", "day", "week")) {
      sec_per <- c(second = 1, minute = 60, hour = 3600,
                   day = 86400, week = 604800)[[u]]
      out[ok] <- origin_ct + as.difftime(t[ok] * sec_per, units = "secs")
    } else {
      out[ok] <- .calendar_time_2_date(t[ok], origin_ct, u)
    }
  }

  if (is_date) as.Date(out, tz = tz) else out
}

#' Convert numeric model time back to dates
#'
#' Inverse of [date_2_time()]: given a numeric time scale measured since a
#' reference origin (for example the `t` column of an `admove_tags` object),
#' return the corresponding dates / date-times.
#'
#' The conversion mirrors [date_2_time()] exactly. Fixed-length units
#' (`second`, `minute`, `hour`, `day`, `week`) use a constant number of seconds
#' per unit. Units `month` and `year` are calendar-aware (their length varies),
#' so the inverse is obtained by numerically inverting [date_2_time()] itself;
#' the round trip `date_2_time(time_2_date(t, tref), tref)` recovers `t` to
#' machine precision, and original dates are recovered exactly wherever
#' [date_2_time()] is injective (e.g. day-resolution data).
#'
#' @param x A numeric vector of times since the origin, or an `admove` object
#'   carrying a time vector: `admove_tags` (its `t`), `admove_data` (its
#'   `tags$t`), or `admove_cov` (its layer times).
#' @param tref A time reference supplying `origin` and `units`. One of an
#'   `admove_tref`, an `admove` object (its [tref()] is used), a list with
#'   elements `origin` and `units`, or a bare `Date`/`POSIXt` origin (units then
#'   default to `"day"`). When `x` is an `admove` object and `tref` is `NULL`,
#'   `tref(x)` is used.
#' @param ... Further arguments passed to methods.
#'
#' @return A `POSIXct` vector (or `Date` when `origin` is a `Date`) the same
#'   length as the input time vector, with `NA` preserved.
#'
#' @seealso [date_2_time()]
#'
#' @examples
#' tr <- create_tref("2003-01-01 UTC", units = "months")
#' time_2_date(c(9.33, 138.16), tref = tr)
#'
#' ## round-trips date_2_time()
#' d <- as.Date(c("2003-10-11", "2014-07-07"))
#' time_2_date(date_2_time(d, tref = tr), tref = tr)
#'
#' @name time_2_date
#' @export
time_2_date <- function(x, tref = NULL, ...) UseMethod("time_2_date")

#' @rdname time_2_date
#' @export
time_2_date.numeric <- function(x, tref = NULL, ...) {
  ref <- .resolve_tref(tref)
  .time_2_date_core(x, ref$origin, ref$units)
}

#' @rdname time_2_date
#' @export
time_2_date.default <- function(x, tref = NULL, ...) {
  time_2_date(as.numeric(x), tref = tref, ...)
}

#' @rdname time_2_date
#' @export
time_2_date.admove_tags <- function(x, tref = NULL, ...) {
  if (is.null(tref)) tref <- tref(x)
  time_2_date(as.numeric(x$t), tref = tref, ...)
}

#' @rdname time_2_date
#' @export
time_2_date.admove_data <- function(x, tref = NULL, ...) {
  if (is.null(tref)) tref <- tref(x)
  time_2_date(as.numeric(x$tags$t), tref = tref, ...)
}

#' @rdname time_2_date
#' @export
time_2_date.admove_cov <- function(x, tref = NULL, ...) {
  if (is.null(tref)) tref <- tref(x)
  time_2_date(as.numeric(dimnames(x)[[3]]), tref = tref, ...)
}


## Collapse a vector of identifiers (tag ids, grid cells, time slices, ...) into
## a short comma-separated string for messages, showing at most 'max_show'
## entries followed by a count of the remainder.
.format_ids <- function(ids, max_show = 5) {
  ids <- as.character(ids)
  if (length(ids) > max_show) {
    paste0(paste(ids[seq_len(max_show)], collapse = ", "),
           ", ... and ", length(ids) - max_show, " more")
  } else {
    paste(ids, collapse = ", ")
  }
}


group_consecutive_ranges <- function(x){
  x <- sort(unique(x))
  breaks <- c(0, which(diff(x) != 1), length(x))

  ranges <- mapply(function(i, j) {
    start <- x[i + 1]
    end <- x[j]
    if (start == end) {
      as.character(start)
    } else {
      paste0(start, ":", end)
    }
  }, breaks[-length(breaks)], breaks[-1])

  res <- paste0("c(",paste(ranges,collapse = ","), ")")

  return(res)
}

t2index <- function(time, time_vec, period = NULL, seasonal = FALSE){
  if (seasonal) {
    if (is.null(period)) stop("period has to be defined in tref(x) for seasonal splines to be used.")
    time <- time %% period
  }
  findInterval(time, time_vec, rightmost.closed = TRUE, left.open = TRUE)
}


## Number of spline slices (seasons) used by each covariate.
##
## alpha/beta/gamma share a single third dimension, sized by the covariate with
## the most breaks, but each covariate is indexed through its own
## dat$time_spline[[i]] (see t2index()). A covariate with fewer breaks therefore
## only ever evaluates its first .get_nsea(dat)[i] slices; the remaining slices
## never enter the likelihood and must stay fixed in the map.
.get_nsea <- function(dat) {
  ts <- dat$time_spline
  if (is.null(ts) || length(ts) == 0L) return(1L)
  n <- vapply(ts, length, integer(1L))
  n[!is.finite(n) | n < 1L] <- 1L
  n
}

build_time <- function(t_obs,
                          mode = c("fill_gaps", "fixed_dt"),
                          dt_min,
                          dt = dt_min,
                          eps = 1e-3) {
  mode <- match.arg(mode)

  t_obs <- as.numeric(t_obs)
  if (length(t_obs) < 2) stop("Need at least 2 observation times.")
  if (any(!is.finite(t_obs))) stop("Non-finite times in t_obs.")
  if (is.unsorted(t_obs, strictly = FALSE)) t_obs <- sort(t_obs)

  t0 <- t_obs[1]
  t1 <- t_obs[length(t_obs)]

  if (mode == "fixed_dt") {

  if (!is.finite(dt) || dt <= 0) stop("dt must be > 0.")

    ts <- seq(t0, t1, by = dt)
    if (tail(ts, 1) < t1 - eps) ts <- c(ts, t1)

    ## map each observed time to nearest grid point (or error if too far)
    ## If you prefer: findInterval + 1 (next), but that shifts forward.
    idx <- vapply(t_obs[-1], function(to) which.min(abs(ts - to)), integer(1))
    if (any(abs(ts[idx] - t_obs[-1]) > eps)) {
      stop("Some observation times do not fall on the fixed_dt grid (increase eps or change dt/mode).")
    }
    observed <- idx

  } else { ## fill_gaps
    ## start with all observed times; add interior points in big gaps only.
    ## unique(): candidate locations of one ambiguous observation may share a
    ## time, and a repeated time would otherwise leave a dt = 0 step in ts.
    ts <- sort(unique(t_obs))

  if (!is.finite(dt_min) || dt_min <= 0) stop("dt_min must be > 0.")

    tu <- ts
    gaps <- diff(tu)
    if (any(gaps > dt_min + eps)) {
      extra <- unlist(lapply(seq_along(gaps), function(i) {
        if (gaps[i] <= dt_min + eps) return(numeric(0))
        ## insert points: tu[i] + dt_min, ..., strictly before tu[i+1]
        if ((tu[i] + dt_min) < (tu[i + 1] - dt_min)) {
          seq(tu[i] + dt_min, tu[i + 1] - dt_min, by = dt_min)
        } else {
          return(numeric(0))
        }
      }), use.names = FALSE)

      ts <- sort(unique(c(ts, extra)))
    }

    ## observed times are in ts exactly (up to floating eps), so map by matching
    ## Use nearest-match to avoid floating-point equality headaches
    observed <- vapply(t_obs[-1], function(to) which.min(abs(ts - to)), integer(1))
    if (any(abs(ts[observed] - t_obs[-1]) > eps)) {
      stop("Failed to match some observation times onto ts (tolerance too small?).")
    }
  }

  dts <- diff(ts)
  list(ts = ts, dts = dts, nts = length(ts), observed = observed)
}

##' Continuous-time Markov chain generator matrices
##'
##' @description
##' Builds the continuous-time Markov chain (CTMC) generator matrices that
##' describe fine-scale movement between grid cells, one per prediction time
##' slice in \code{fit$dat$pred$time}.
##'
##' @param fit A fitted model object of class \code{"admove"} that already
##'   carries prediction quantities in \code{fit$pred} (see
##'   [add_predictions()]).
##'
##' @details
##' The result is a list of length \code{nt = length(fit$dat$pred$time) - 1},
##' one sparse generator per prediction time slice. Element \code{t} is an
##' \code{nc x nc} sparse matrix (\code{Matrix} \code{"dgCMatrix"}, where
##' \code{nc} is the number of grid cells) holding the generator \eqn{Q_t}
##' evaluated from the diffusion, taxis, and (optionally) advection fields at
##' prediction time \code{t}. A sparse representation is used because each row
##' has at most five non-zero entries (the cell itself plus its four
##' neighbours), so storage is \eqn{O(nc)} rather than \eqn{O(nc^2)}; this
##' mirrors the sparse assembly used inside the likelihood.
##'
##' Off-diagonal entries \code{mstar[[t]][i, j]} are the instantaneous rates of
##' moving from cell \code{i} to a neighbouring cell \code{j}, in units of
##' \strong{1 / time} (the reciprocal of the model time unit, e.g. per day).
##' Diagonal entries are the negative row sums, so each slice is a valid CTMC
##' generator. The rates do \emph{not} include the prediction time step: to
##' obtain the movement (transition) probability matrix over a time step
##' \code{dt}, exponentiate the generator, e.g.
##' \code{Matrix::expm(mstar[[t]] * dt)}.
##'
##' The generator is assembled exactly as in the likelihood ([nll()]): diffusion
##' contributes \code{D / next_dist^2} to each neighbour, while taxis and
##' advection are added with the same upwind scheme via \code{fill_inst_mat()} (so
##' off-diagonal rates are guaranteed non-negative). The only difference from
##' \code{nll()} is that the prediction time step \code{dt} is \emph{not} folded
##' in, leaving a pure per-time-unit generator. The taxis drift velocity is
##' \eqn{\kappa \nabla h} (already stored, with \eqn{\kappa} folded in, as
##' \code{fit$pred$hTdx} / \code{hTdy}); it is \emph{not} scaled by the diffusion
##' coefficient.
##'
##' @return A list of length \code{nt} of sparse (\code{"dgCMatrix"}) CTMC
##'   generator matrices, each \code{nc x nc} in units of 1 / time.
##'
##' @seealso [add_predictions()], [add_tag_dist()]
##'
##' @keywords internal
calc_mstar <- function(fit) {

  grid <- fit$dat$pred$grid
  nc <- nrow(grid$xygrid)
  nt <- length(fit$dat$pred$time) - 1
  nextTo <- get_neighbours(grid)
  next_dist <- c(grid$cellsize[1], grid$cellsize[1],
                 grid$cellsize[2], grid$cellsize[2])

  mstar_template <- make_mstar_template(nextTo)

  mstar <- vector("list", nt)
  neg_slices <- integer(0)

  for (t in 1:nt) {

    Zstar <- Astar <- Dstar <- mstar_template
    Zstar@x[] <- Astar@x[] <- Dstar@x[] <- 0

    ## taxis
    if (fit$conf$use_taxis) {
      move <- cbind(fit$pred$hTdx[, t], fit$pred$hTdy[, t])  ## distance / time
      Zstar <- fill_inst_mat(Zstar, move, nextTo, next_dist, fit$conf$drift_scheme)
    }

    ## advection
    if (fit$conf$use_advection) {
      move <- cbind(fit$pred$hAx[, t], fit$pred$hAy[, t])  ## distance / time
      Astar <- fill_inst_mat(Astar, move, nextTo, next_dist, fit$conf$drift_scheme)
    }

    ## diffusion
    D <- exp(fit$pred$hD[, t])  ## distance^2 / time
    for (k in 1:4) {
      j <- k + 1
      ind <- which(!is.na(nextTo[, j]))
      Dstar[cbind(ind, nextTo[ind, j])] <- D[ind] / next_dist[k]^2
    }

    ## movement rates (per unit time)
    Mstar <- Zstar + Astar + Dstar

     ## fill_inst_mat is upwind
     ## (na.rm: a prediction grid/time outside covariate coverage yields NA
     ## rates; those propagate into the generator rather than crashing this check)
    if (any(Mstar@x < 0, na.rm = TRUE)) neg_slices <- c(neg_slices, t)

    ## mass balance on the diagonal
    Mstar[cbind(1:nc, 1:nc)] <- 0
    Mstar[cbind(1:nc, 1:nc)] <- -Matrix::rowSums(Mstar)

    mstar[[t]] <- Mstar
  }

  if (length(neg_slices)) {
    warning("calc_mstar(): negative off-diagonal generator rates in ", length(neg_slices),
            " time slice(s) (", .format_ids(neg_slices), ")",
            "; the CTMC generator is invalid there and expm() may yield negative ",
            "probabilities. This usually means drift dominates diffusion at the ",
            "current grid resolution (grid-Peclet > 1); consider a finer grid.",
            call. = FALSE)
  }

  return(mstar)
}


fill_inst_mat <- function(mat, move, nextTo, next_dist, scheme = "upwind") {
  xyind <- c(2, 2, 1, 1)
  dirsign <- c(+1, -1, -1, +1)
  ## NULL-safe (e.g. conf from objects built before drift_scheme existed):
  ## anything other than "central" is treated as the default upwind scheme.
  central <- identical(scheme, "central")
  ## AD-safe pmax (upwind keeps only outflow in the drift direction)
  pos <- function(x) 0.5 * (x + abs(x))
  ## 4 neighbours
  for (k in 1:4) {
    j <- k + 1
    ind <- which(!is.na(nextTo[, j]))
    v <- dirsign[k] * move[ind, xyind[k]]
    ## upwind: pos(v); central: v/2 split symmetrically (can be negative)
    rate <- if (central) 0.5 * v else pos(v)
    mat[cbind(ind, nextTo[ind, j])] <- rate / next_dist[k]
  }
  return(mat)
}

get_par_est <- function(par, map, opt) {

  ## basic checks
  stopifnot(is.list(par), is.list(map), is.list(opt), !is.null(opt$par))
  res <- par

  ## helper: convert map to integer vector with NAs preserved
  map_to_int <- function(m) {
    if (is.null(m)) return(NULL)

    if (is.factor(m)) {
      ## try to interpret factor labels as integers (more robust than codes)
      m_chr <- as.character(m)
      mi <- suppressWarnings(as.integer(m_chr))

      ## fallback: use factor codes if labels aren't numeric
      if (any(!is.na(m_chr) & is.na(mi))) mi <- as.integer(m)

      return(mi)
    }

    suppressWarnings(as.integer(m))
  }

  for (nm in names(res)) {

    if (!nm %in% names(map)) next
    mi <- map_to_int(map[[nm]])
    if (is.null(mi)) next

    arr <- res[[nm]]

    ## ensure map length matches parameter length (linear indexing)
    if (length(mi) != length(arr)) {
      stop(sprintf(
        "Parameter '%s': length(map)=%d does not match length(par)=%d.",
        nm, length(mi), length(arr)
      ))
    }

    idx <- which(!is.na(mi))
    if (length(idx) == 0L) next

    ## pull the estimated values for this parameter name
    est <- opt$par[names(opt$par) == nm]
    if (length(est) == 0L) next

    levs <- sort(unique(mi[idx]))

    ## align estimates to map levels
    if (length(est) == length(levs)) {
      ## if levels aren't 1..k, map them explicitly to the available estimates
      est_by_level <- setNames(as.numeric(est), as.character(levs))
    } else if (max(mi[idx]) <= length(est)) {
      ## assume levels are 1..k and estimates are in that order
      est_by_level <- setNames(as.numeric(est), as.character(seq_along(est)))
    } else {
      stop(sprintf(
        "Parameter '%s': map refers up to level %d but only %d estimate(s) found in opt$par.",
        nm, max(mi[idx]), length(est)
      ))
    }

    vals <- est_by_level[as.character(mi[idx])]
    if (anyNA(vals)) {
      stop(sprintf("Parameter '%s': could not match some map levels to estimates.", nm))
    }

    arr[idx] <- vals
    res[[nm]] <- arr
  }

  res
}



get_pretty_probs <- function(n) {
  if (n == 1) {
    0.5
  } else if (n == 2) {
    c(0.25, 0.75)
  } else if (n == 3) {
    c(0.05, 0.5, 0.95)
  } else if (n == 4) {
    c(0.05, 0.3, 0.7, 0.95)
  } else stop("pretty_probs only implemented for upt to 4 knots")
}


slice_first2 <- function(a, at = 1L) {
  d <- dim(a)
  if (is.null(d) || length(d) < 2L) {
    stop("`a` must have at least 2 dimensions.")
  }

  if (length(d) == 2L) return(a)

  if (length(at) == 1L) at <- rep.int(at, length(d) - 2L)
  if (length(at) != (length(d) - 2L)) {
    stop("`at` must have length 1 or ndims(a)-2.")
  }

  out <- do.call(
    `[`,
    c(list(a), list(TRUE), list(TRUE), as.list(at), list(drop = FALSE))
  )

  dim(out) <- d[1:2]
  return(out)
}

list_2_3Darray <- function(x) {

  if (!inherits(x, "list")) stop("x is not a list!")

  n <- length(x)

  if (length(dim(x[[1]])) > 2 && any(dim(x[[1]])[-c(1,2)] > 1)) stop("More than 2 dimensions in the list elements. Not sure how to simplify that.")

  if (any(apply(sapply(x, dim), 1, diff) > 0)) stop("The dimensions of the list elements are not consistent. Not sure how to simplify that.")

  if (!is.null(names(x))) {
    time_names <- names(x)
  } else time_names <- rep(NA, n)

  res <- array(NA, dim = c(dim(x[[1]])[c(1,2)],n))
  for (i in 1:n) {
    res[,,i] <- slice_first2(x[[i]])
    if (length(dimnames(x[[i]])) > 2 && !is.null(dimnames(x[[i]])[[3]]) && is.na(time_names[i])) time_names[i] <- dimnames(x[[i]])[[3]]
  }

  dimnames(res) <- list(dimnames(x[[1]])[[1]],
                        dimnames(x[[1]])[[2]],
                        time_names)

  return(res)
}

guard_neg <- function(x, eps = 1e-12) {
  0.5 * (x + sqrt(x*x + eps*eps))
}

make_mstar_template <- function(nextTo, ad = FALSE) {

  nc <- nrow(nextTo)

  from <- as.vector(matrix(1:nc, nc, 5))
  to <- as.vector(nextTo)
  keep <- !is.na(to)
  Iall <- from[keep]
  Jall <- to[keep]

  template <- Matrix::sparseMatrix(i = Iall, j = Jall, x = 1, dims = c(nc, nc))

  if (ad) {
    template <- RTMB::AD(template)
  }

  return(template)
}



##' Create x- and y-coordinate covariate fields
##'
##' @description
##' `make_x_y_cov()` creates two simple covariate fields representing the spatial
##' x and y coordinates of a grid. The function returns these as a list of
##' covariate arrays that can be used as input to \emph{admove}, for example as
##' simple spatial trend covariates or for testing and demonstration purposes.
##'
##' The first covariate varies along the x dimension and the second varies along
##' the y dimension. Both covariates are created on the spatial domain and cell
##' centres of the supplied grid.
##'
##' @param grid A grid of class `admove_grid`, as returned by [create_grid()].
##' @param tref Optional time reference information to attach to the returned
##'   covariates. Default: `NULL`.
##'
##' @return
##' An object of class `admove_cov_list` with two covariate fields, `x` (varies
##' along the x dimension) and `y` (varies along the y dimension).
##'
##' @details
##' Both covariates are created as single-time-slice fields with `times = 0`.
##' The spatial reference is always taken from `grid`, so that the covariates
##' match the grid in [setup_data()]; if the grid has no spatial reference, add
##' one with [add_sref()] first. A grid carries no time reference, so it can be
##' attached via `tref`.
##'
##' @examples
##' grid <- create_grid(verbose = FALSE)
##' xy_cov <- make_x_y_cov(grid)
##'
##' @export
make_x_y_cov <- function(grid, tref = NULL) {

  if (!inherits(grid, "admove_grid")) {
    stop("'grid' must be an 'admove_grid' (see create_grid()), not an object of class '",
         paste(class(grid), collapse = "/"), "'.")
  }

  nx <- dim(grid)[1]
  ny <- dim(grid)[2]

  cov1 <- prep_cov(matrix(1:nx, nx, ny),
                   ## matrix(seq(0.5,nx,1) - nx/2, nx, ny), ## not working
                   x_centers = x_centers(grid),
                   y_centers = y_centers(grid),
                   times = 0,
                   sref = sref(grid),
                   tref = tref)
  cov2 <- prep_cov(matrix(1:ny, nx, ny, byrow = TRUE),
                   ## matrix(seq(0.5,ny,1) - ny/2, nx, ny, byrow = TRUE),
                   x_centers = x_centers(grid),
                   y_centers = y_centers(grid),
                   times = 0,
                   sref = sref(grid),
                   tref = tref)
  cov <- list(x = cov1, y = cov2)
  cov <- .add_class(cov, "admove_cov_list")
  cov <- add_sref(cov, sref(grid))
  cov <- add_tref(cov, tref)

  cov
}



.get_non_na_from_map <- function(par, map) {
  nms <- intersect(names(par), names(map))

  unlist(lapply(nms, function(nm) {
    x <- as.vector(par[[nm]])
    keep <- !is.na(map[[nm]])

    if (length(x) != length(keep)) {
      stop(sprintf(
        "Length mismatch for '%s': par has %d values, map has %d",
        nm, length(x), length(keep)
      ))
    }
    y <- x[keep]
    names(y) <- rep(nm, length(y))
    y
  }), use.names = TRUE)
}



## NA -> 0, for picking a representative candidate by probability
.na_zero <- function(x) {
  if (is.null(x)) return(0)
  x[is.na(x)] <- 0
  x
}


## Observation-event index of one tag, as used inside nll().
##
## Rows sharing an event are mutually exclusive candidate positions for a single
## observation. Tags without the column (the common case, and any hand-built
## test fixture) get one event per row, which is what the likelihood assumed
## before ambiguous positions existed.
.tag_events <- function(tag) {
  ev <- tag[["event"]]
  if (is.null(ev)) return(seq_len(nrow(tag)))
  as.integer(match(ev, unique(ev)))
}


## Numerically stable log(sum(exp(x))) for a short AD vector.
##
## Used for the finite mixture over candidate observation locations (see nll()).
## A plain log(sum(exp(x))) underflows to -Inf as soon as one candidate is far
## enough away that its log-density is very negative, which would poison the
## objective and every gradient. RTMB::logspace_add(a, b) = log(exp(a) + exp(b))
## is both AD-differentiable and stable, so folding over it keeps the whole
## accumulation on the tape. max()/which.max() are deliberately avoided: they are
## non-smooth and would need the values off the tape.
.logsumexp_ad <- function(x) {
  n <- length(x)
  if (n == 0L) return(NULL)
  if (n == 1L) return(x[1])
  out <- x[1]
  for (k in 2:n) out <- RTMB::logspace_add(out, x[k])
  out
}


## Natural cubic spline through the knots (xp, yp), returning either the value
## function (deriv = FALSE) or its analytic first derivative (deriv = TRUE).
##
## `yp` are the estimated function values at the knots (the parameters), so the
## spline interpolates them exactly. Boundary conditions are "natural" (second
## derivative zero at the outer knots), which makes the spline extrapolate
## LINEARLY beyond the knot range - far more robust in the covariate tails than a
## global polynomial.
##
## Implemented with plain automatic-differentiation operations only - no
## comparisons or branching on the (possibly AD) evaluation point, which RTMB
## forbids. The spline is written in a truncated-power basis
##
##   S(x) = a + b * x + sum_i c_i * pos(x - t_i)^3,   pos(z) = max(0, z),
##
## subject to the two "natural" constraints  sum_i c_i = 0  and  sum_i c_i t_i = 0,
## which force the second derivative to vanish at the ends and hence make S linear
## beyond the outer knots. pos(z) = 0.5 * (z + abs(z)) is AD-safe (uses abs, not a
## comparison), so the whole basis is a single smooth formula valid for every x.
##
## The coefficients theta = (a, b, c_1..c_n) solve a fixed (n+2) linear system
## whose matrix depends only on the knots (numeric); the right-hand side is linear
## in the knot values yp, so RTMB::solve(numeric matrix, AD rhs) is AD-safe -
## exactly the pattern used by the legacy Vandermonde solve.
##
## Retained as the "natural" method. The default is now "rtmb", which delegates
## to RTMB::splinefun and agrees with this construction to machine precision;
## this pure-R version needs no atomic and stays available as a cross-check and
## as a fallback (see .poly_fun).
.natural_spline_fun <- function(xp, yp, deriv = FALSE) {

  "c" <- RTMB::ADoverload("c")

  n <- length(xp)

  pos <- function(z) 0.5 * (z + abs(z))              # max(0, z), AD-safe

  ## Two knots -> straight line (natural spline degenerates to linear)
  if (n == 2) {
    slope <- (yp[2] - yp[1]) / (xp[2] - xp[1])
    if (!deriv) {
      return(function(x) yp[1] + slope * (x - xp[1]))
    } else {
      return(function(x) rep(slope, length(x)))
    }
  }

  ## Linear system for theta = (a, b, c_1..c_n):
  ##   interpolation:  a + b t_j + sum_i c_i (t_j - t_i)_+^3 = y_j   (j = 1..n)
  ##   natural end 1:  sum_i c_i       = 0
  ##   natural end 2:  sum_i c_i t_i   = 0
  G <- matrix(0, n + 2, n + 2)
  G[1:n, 1]           <- 1                            # a
  G[1:n, 2]           <- xp                           # b * t_j
  for (j in seq_len(n)) {                             # c_i * (t_j - t_i)_+^3
    G[j, 2 + seq_len(n)] <- pmax(0, xp[j] - xp)^3
  }
  G[n + 1, 2 + seq_len(n)] <- 1                       # sum c_i = 0
  G[n + 2, 2 + seq_len(n)] <- xp                      # sum c_i t_i = 0

  rhs <- c(yp, 0, 0)
  theta <- RTMB::solve(G, rhs)                        # AD vector, length n + 2
  a  <- theta[1]
  b  <- theta[2]
  cc <- theta[2 + seq_len(n)]                         # cubic coefficients

  if (!deriv) {
    function(x) {
      out <- a + b * x
      for (i in seq_len(n)) out <- out + cc[i] * pos(x - xp[i])^3
      out
    }
  } else {
    function(x) {
      out <- b + 0 * x                                # keep length(x)
      for (i in seq_len(n)) out <- out + cc[i] * 3 * pos(x - xp[i])^2
      out
    }
  }
}


## Build a preference smooth (and its first derivative) from knot locations `xp`
## and knot values `yp`. `yp` are the estimated function values at the knots, so
## the smooth always interpolates them exactly. `method` selects the construction:
##   "rtmb"    - RTMB::splinefun natural cubic spline. Default. Same spline as
##               "natural" (they agree to machine precision), but evaluated by
##               RTMB's atomic, whose cost is independent of the number of knots.
##   "natural" - the same natural cubic spline written out in a truncated-power
##               basis in plain R (see .natural_spline_fun); no atomic involved.
##   "poly"    - legacy single global interpolating polynomial (Vandermonde solve);
##               retained for reproducibility of older fits.
## The `adv` and single-knot cases are method-independent (linear-through-origin and
## constant respectively).
.poly_fun <- function(xp, yp, deriv = FALSE, adv = FALSE, method = "rtmb") {

  if (is.null(method)) method <- "rtmb"

  if (!adv && length(xp) > 1 && all(diff(xp) == 0)) return(NULL)

  if (adv) {
    ## Simple linear case
    val <- yp[1]
    f <- function(x) x * val
    df <- function(x) rep(val, length(x))
  } else {
    if (length(xp) == 1) {
      ## For one knot return parameter (e.g. constant diffusion)
      val <- yp[1]
      f <- function(x) val
      df <- function(x) rep(0, length(x))
    } else if (method == "rtmb") {
      ## Same natural cubic spline as "natural", evaluated by RTMB's atomic.
      S <- RTMB::splinefun(xp, yp, method = "natural")
      f <- function(x) S(x)
      ## S(x, deriv = 1) has no analytic derivative: it tapes the spline again
      ## and differentiates that nested tape, on EVERY call. That tape depends
      ## only on length(x), never on the values, so build it once per smooth and
      ## replay it - ~5x faster over a whole fit. Taping it here also keeps the
      ## evaluation point a bare tape variable rather than an inline expression,
      ## which is what lets this branch work on RTMB builds predating the
      ## MakeTape promise fix (kaskr/RTMB "Fix #90"). as.vector() because a
      ## replayed tape returns a 1 x n matrix.
      D <- NULL
      df <- function(x) {
        ## Off-tape (plotting, simulation) nothing needs caching: the spline
        ## falls back to stats::splinefun and differentiates directly, whereas
        ## taping would build one tape the size of the whole covariate field.
        ## That fallback errors on NA/NaN in its deriv > 0 branch though, and
        ## the covariate interpolant returns NaN off the field, so evaluate the
        ## finite points only and leave the gaps missing.
        if (!inherits(x, "advector") && !inherits(yp, "advector")) {
          out <- rep(NA_real_, length(x))
          ok <- !is.na(x)
          if (any(ok)) out[ok] <- as.vector(S(x[ok], deriv = 1L))
          return(out)
        }
        if (is.null(D) || length(D$par()) != length(x))
          D <<- RTMB::MakeTape(function(z) S(z, deriv = 1L), numeric(length(x)))
        as.vector(D(x))
      }
    } else if (method == "natural") {
      ## Natural cubic spline through the knot values (see .natural_spline_fun).
      f <- .natural_spline_fun(xp, yp, deriv = FALSE)
      df <- .natural_spline_fun(xp, yp, deriv = TRUE)
    } else if (method == "poly") {
      ## Solve for polynomial coefficients
      n <- length(xp)
      A <- outer(xp, 0:(n-1), "^")
      alpha <- RTMB::solve(A, yp)

      f <- function(x){
        ## Evaluate polynomial: sum(alpha[j+1] * x^j)
        v <- outer(x, 0:(n-1), "^")
        as.vector(v %*% alpha)
        ## as.vector(alpha[1] + sum(alpha[-1] * x^(1:(n-1))))
      }

      df <- function(x){
        ## Evaluate derivative: sum(j * alpha[j+1] * x^(j-1))
        if (n == 2) {
          rep(alpha[2], length(x))  # Linear case
        }else{
          j <- 1:(n-1)
          v <- outer(x, j - 1, "^")
          as.vector(v %*% (j * alpha[-1]))
        }
        ## as.vector(alpha[2] +
        ##           sum(alpha[-c(1:2)] *
        ##               (2:(n-1)) * x^(1:(n-2))))
      }
    } else {
      stop("Unknown smooth method '", method,
           "' (expected 'rtmb', 'natural' or 'poly').")
    }
  }

  if (!deriv) return(f) else return(df)
}



.make_pref_funcs <- function(alpha, beta, gamma,
                            knots_tax, knots_dif,
                            method = "rtmb") {

  if (is.null(method)) method <- "rtmb"

  ncov <- dim(knots_tax)[2]

  pref_funcs <- vector("list", 8)
  for (i in 1:8) pref_funcs[[i]] <- vector("list", ncov)
  names(pref_funcs) <- c("dif", "ddif",
                         "tax", "dtax",
                         "adv_x", "dadv_x",
                         "adv_y", "dadv_y")

  for(i in 1:ncov){

    ## advection
    pref_funcs$tax[[i]] <- pref_funcs$dtax[[i]] <- vector("list", dim(alpha)[3])
    for(j in 1:dim(alpha)[3]){
      pref_funcs$tax[[i]][[j]] <- .poly_fun(knots_tax[,i], alpha[,i,j],
                                          method = method)
      pref_funcs$dtax[[i]][[j]] <- .poly_fun(knots_tax[,i],
                                          alpha[,i,j],
                                          deriv = TRUE,
                                          method = method)
    }

    ## diffusion
    pref_funcs$dif[[i]] <- pref_funcs$ddif[[i]] <- vector("list", dim(beta)[3])
    for(j in 1:dim(beta)[3]){
      pref_funcs$dif[[i]][[j]] <- .poly_fun(knots_dif[,i], beta[,i,j],
                                          method = method)
      pref_funcs$ddif[[i]][[j]] <- .poly_fun(knots_dif[,i], beta[,i,j],
                                          deriv=TRUE,
                                          method = method)
    }

    ## advection (x)
    pref_funcs$adv_x[[i]] <- pref_funcs$dadv_x[[i]] <- vector("list", dim(gamma)[3])
    for(j in 1:dim(gamma)[3]){
      pref_funcs$adv_x[[i]][[j]] <- .poly_fun(NULL, gamma[1,i,j], adv = TRUE)
    }

    ## advection (y)
    pref_funcs$adv_y[[i]] <- pref_funcs$dadv_y[[i]] <- vector("list", dim(gamma)[3])
    for(j in 1:dim(gamma)[3]){
      pref_funcs$adv_y[[i]][[j]] <- .poly_fun(NULL, gamma[2,i,j], adv = TRUE)
    }
  }

  return(pref_funcs)
}


.get_tag_type_integer <- function(x) {
  if (is.numeric(x)) return(as.integer(x))
  tag_types <- c("d","s","c","a")
  as.integer(factor(x, levels = tag_types))
}


.get_engine_integer <- function(x) {
  if (is.numeric(x)) return(as.integer(x))
  engine_types <- c("KF","CTMC")
  as.integer(factor(x, levels = engine_types))
}
