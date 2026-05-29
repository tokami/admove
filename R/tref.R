
## Functions ------------------------------------------------------------------------


##' Time reference for admove objects
##'
##' A lightweight container that stores the time reference information used
##' throughout \pkg{admove}. It bundles:
##' \itemize{
##'   \item \strong{origin}: a POSIXct date-time used as the temporal origin
##'   \item \strong{units}: the units of the numeric time axis stored in objects
##'   \item \strong{period}: optional seasonal cycle length in the same units as
##'     the stored times
##' }
##'
##' The intention is that \emph{objects store their times on a simple numeric
##' axis} (e.g. months since \code{origin}), while \code{origin} and \code{units}
##' define how these values should be interpreted and displayed.
##'
##' If \code{period} is provided (or can be inferred), it can be used to create
##' seasonally repeating covariate fields or spline bases by wrapping time with
##' modulo arithmetic. For example, if times are stored in months, a natural
##' annual seasonal cycle is \code{period = 12}. For custom discretizations
##' (e.g. 10 time steps per year), set \code{period = 10} and \code{units = "custom"}.
##'
##' If \code{units} is missing but \code{period} is supplied, \code{units} can be
##' inferred for common annual discretizations (\code{period = 1, 2, 4, 12, 52})
##' and otherwise defaults to \code{"custom"}.
##'
##' @param origin A date-time corresponding to the temporal origin (e.g. a release
##'   date). Will be converted to \code{POSIXct}. Use \code{NA} to leave undefined.
##' @param units Character string describing the units of the stored numeric time
##'   axis, e.g. \code{"year"}, \code{"quarter"}, \code{"month"}, \code{"week"},
##'   \code{"day"}. Use \code{list_units_time()} for supported options.
##' @param period Optional numeric seasonal period (cycle length) in the same
##'   units as the stored time axis. If \code{NULL}, a default is inferred for
##'   common annual discretizations when \code{units} is known.
##'
##' @return An object of class \code{admove_tref}.
##'
##' @examples
##' ## Monthly time axis with annual seasonality (12 months)
##' tr <- create_tref(origin = as.Date("2025-01-01"), units = "month")
##' tr
##'
##' ## Infer units from period
##' tr2 <- create_tref(origin = as.Date("2025-01-01"), period = 12)
##' tr2
##'
##' ## Custom discretization: 10 time steps per year
##' tr10 <- create_tref(origin = as.Date("2025-01-01"), period = 10)
##' tr10
##'
##' @export
create_tref <- function(origin = NA, units = NA_character_, period = NULL) {

  ## origin handling
  if (is.null(origin) || (length(origin) == 1 && is.na(origin))) {
    origin <- as.POSIXct(NA)
  } else {
    origin <- as.POSIXct(origin)
  }

  ## validate / infer period
  if (!is.null(period)) {
    period <- as.numeric(period)
    if (length(period) != 1 || !is.finite(period) || period <= 0) {
      stop("period must be a single positive finite number.")
    }
  }

  ## infer units if missing but period given
  u_norm <- .normalise_time_unit(units)
  if ((is.na(u_norm) || u_norm == "") && !is.null(period)) {
    units <- .infer_units_from_period(period)
    u_norm <- .normalise_time_unit(units)
  }

  ## infer period if missing but units known
  if (is.null(period)) {
    period <- .default_period_from_units(u_norm)
  }

  ## optional: warn if both given but look inconsistent for the classic cases
  if (!is.null(period) && !is.na(u_norm)) {
    p_def <- .default_period_from_units(u_norm)
    if (!is.na(p_def) && !isTRUE(all.equal(period, p_def))) {
      warning("create_tref(): units='", units, "' usually implies period=", p_def,
              ", but period=", period, " was provided. Using period=", period, ".")
    }
  }

  structure(
    list(origin = origin,
         units = units,
         period = period),
    class = "admove_tref"
  )
}



##' Add or harmonize a time reference on an object
##'
##' Attach a time reference (\code{tref}) to an object, or harmonize an existing
##' time reference with that of another object or tref-like specification.
##'
##' This function is intended for adding or updating \emph{metadata} describing a
##' numeric time axis, such as its origin, units, and period. If the object
##' already has a tref and the requested tref has different time units, the stored
##' time values are rescaled via \code{\link{scale_tref}} where possible.
##'
##' If the object already has a tref and the requested origin differs, changing
##' the origin without changing the stored time values would change the meaning of
##' the time series. Therefore, by default, the function throws an error in this
##' case. To preserve the represented dates while adopting the new origin, set
##' \code{shift_origin = TRUE}; this calls \code{\link{shift_tref}} internally to
##' shift the stored time values before updating the tref metadata.
##'
##' The argument \code{tref} may be a tref-like object, or any object for which
##' \code{tref(tref)} can be extracted.
##'
##' @param x An object to which a tref should be added or whose tref should be
##'   harmonized.
##' @param tref Optional tref specification. This can be:
##'   \itemize{
##'   \item an object of class \code{"admove_tref"},
##'   \item a named list containing one or more of \code{origin}, \code{units},
##'     and \code{period},
##'   \item another object from which \code{tref(tref)} can be extracted.
##'   }
##'   If \code{NULL}, the existing tref on \code{x} is kept where possible.
##' @param verbose Logical; if \code{TRUE}, informative messages are printed when
##'   tref fields are inferred, retained, or converted.
##' @param shift_origin Logical; if \code{TRUE} and the origin of \code{x} differs
##'   from the requested origin, the stored time values are shifted using
##'   \code{\link{shift_tref}} so that the represented dates remain unchanged.
##'   If \code{FALSE} (default), differing non-missing origins cause an error.
##'
##' @return The input object \code{x}, with updated tref metadata. Depending on
##'   the input and requested tref, the stored numeric time values may also be
##'   rescaled (if units differ) or shifted (if \code{shift_origin = TRUE}).
##'
##' @details
##' The function distinguishes between three different kinds of tref changes:
##' \enumerate{
##'   \item \strong{Adding missing tref metadata}: if \code{x} has no tref, or if
##'     parts of its tref are missing, the missing fields are filled from
##'     \code{tref} where available.
##'   \item \strong{Changing units}: if time units differ between \code{x} and
##'     \code{tref}, the stored numeric time values are converted using
##'     \code{\link{scale_tref}} where possible.
##'   \item \strong{Changing origin}: if origins differ, the stored time values
##'     must also change if the represented dates are to remain the same. This is
##'     not done silently; use \code{shift_origin = TRUE} or call
##'     \code{\link{shift_tref}} directly.
##' }
##'
##' Period is treated strictly: if both \code{x} and \code{tref} define a finite
##' period and the two differ, an error is raised.
##'
##' @examples
##' ## Assume x and y are objects with numeric time and tref metadata
##' ## x <- add_tref(x, list(origin = as.Date("2020-01-01"), units = "day"))
##'
##' ## Copy tref fields from another object
##' ## x <- add_tref(x, y)
##'
##' ## Convert units if needed
##' ## x <- add_tref(x, list(units = "month"))
##'
##' ## Adopt the origin of y while preserving represented dates
##' ## x <- add_tref(x, y, shift_origin = TRUE)
##'
##' @seealso \code{\link{shift_tref}}, \code{\link{scale_tref}},
##'   \code{\link{create_tref}}, \code{\link{tref}}
##'
##' @export
add_tref <- function(x, tref = NULL, verbose = TRUE, shift_origin = FALSE) {

  tref0 <- try(tref(x), silent = TRUE)

  ## extract requested tref components
  origin_new <- units_new <- NULL
  period_new <- NULL

  tref_target <- .get_tref_target(tref)
  if (!is.null(tref_target)) {
    origin_new <- tref_target$origin
    units_new <- tref_target$units
    period_new <- tref_target$period
  }

  if (is.null(origin_new)) origin_new <- NA
  if (is.null(units_new))  units_new <- NA_character_

  ## if already present, enforce consistency / convert if needed
  if (!is.null(tref0) && !inherits(tref0, "try-error")) {

    ## ---- origin consistency ----
    origin0 <- origin(tref0)

    origin0_na <- .is_na_scalar(origin0)
    originN_na <- .is_na_scalar(origin_new)

    ## if new is NA and old is defined -> keep old
    if (originN_na && !origin0_na) {
      origin_new <- origin0
    }

    ## if old is NA and new is defined -> set it
    if (origin0_na && !originN_na && isTRUE(verbose)) {
      message("Setting tref$origin on object (was NA).")
    }

    ## if both defined but differ
    if (!origin0_na && !originN_na && !.same_origin(origin0, origin_new)) {

      if (isTRUE(shift_origin)) {
        if (isTRUE(verbose)) {
          message("tref origins differ; shifting stored time values to match origin of 'tref'.")
        }
        x <- shift_tref(x, origin = origin_new, verbose = FALSE)
      } else {
        stop(
          "tref origins differ. ",
          "Changing origin without shifting the stored time values would change the meaning of the series. ",
          "Use shift_tref(x, tref = ...) or add_tref(..., shift_origin = TRUE)."
        )
      }
    }

    ## ---- units consistency (convert if needed) ----
    units0 <- units_time(tref(x))

    if (.is_na_scalar(units_new) && !.is_na_scalar(units0)) {
      units_new <- units0
    }

    if (!.is_na_scalar(units0) &&
        !.is_na_scalar(units_new) &&
        units0 != units_new) {

      if (exists("scale_tref", mode = "function")) {
        x <- scale_tref(x, units = units_new, verbose = verbose)
      } else {
        stop("tref units differ (", units0, " -> ", units_new,
             "), but scale_tref() is not available.")
      }
    }

    ## ---- period consistency ----
    tref1 <- try(tref(x), silent = TRUE)
    period0 <- tref1$period
    if (is.null(period0)) period0 <- NA_real_

    if (is.null(period_new) && !is.na(period0)) {
      period_new <- period0
    }

    if (!is.null(period_new)) {
      period_new_num <- suppressWarnings(as.numeric(period_new))
      if (length(period_new_num) == 1L &&
          is.finite(period_new_num) &&
          length(period0) == 1L &&
          is.finite(as.numeric(period0)) &&
          !isTRUE(all.equal(as.numeric(period0), period_new_num))) {
        stop("tref period not the same! Please check tref$period of the input object and the specified period.")
      }
    }
  }

  if ((is.null(units_new) || .is_na_scalar(units_new)) && isTRUE(verbose)) {
    message("Consider providing time 'units' (tref$units) and/or 'period' (tref$period).")
  }

  ## only pass period if valid
  period_arg <- NULL
  if (!is.null(period_new)) {
    p <- suppressWarnings(as.numeric(period_new))
    if (length(p) == 1L && is.finite(p) && p > 0) period_arg <- p
  }

  tref(x) <- create_tref(origin = origin_new,
                                 units = units_new, period = period_arg)
  x
}



##' Shift the origin of a time reference while preserving represented dates
##'
##' Change the origin of an object's time reference (\code{tref}) and shift the
##' stored numeric time values accordingly, so that the represented dates or
##' date-times remain unchanged.
##'
##' This function should be used when an object already has a valid tref and you
##' want to adopt a different origin without changing the underlying time points
##' that the numeric values represent. In other words, it changes the coordinate
##' system of the time axis, not the actual dates.
##'
##' The new origin can be supplied directly via \code{origin}, or indirectly via
##' \code{tref} as a tref-like object or another object from which
##' \code{tref(tref)} can be extracted.
##'
##' @param x An object with an existing valid tref and numeric time values stored
##'   in a supported format.
##' @param tref Optional object or tref-like specification providing the new
##'   origin. Ignored if \code{origin} is supplied. This can be:
##'   \itemize{
##'   \item an object of class \code{"admove_tref"},
##'   \item a named list containing at least \code{origin},
##'   \item another object from which \code{tref(tref)} can be extracted.
##'   }
##' @param origin Optional new origin as a \code{Date} or \code{POSIXct}-like
##'   object. If supplied, this takes precedence over \code{tref}.
##' @param verbose Logical; if \code{TRUE}, informative messages are printed about
##'   the shift applied.
##'
##' @return The input object \code{x}, with numeric time values shifted and
##'   \code{tref(x)$origin} replaced by the new origin.
##'
##' @details
##' Let \eqn{t} denote the stored numeric time values, \eqn{o_0} the old origin,
##' and \eqn{o_1} the new origin. To preserve the represented dates, the new time
##' values \eqn{t'} must satisfy:
##' \deqn{o_0 + t = o_1 + t'}
##' and therefore
##' \deqn{t' = t - (o_1 - o_0).}
##'
##' The offset \eqn{o_1 - o_0} is computed in the current tref units. For
##' \code{"second"}, \code{"minute"}, \code{"hour"}, \code{"day"}, and
##' \code{"week"}, fixed-duration differences are used. For \code{"month"} and
##' \code{"year"}, the offset is computed in a calendar-aware way using
##' \code{lubridate::time_length()} on the interval between the two origins.
##'
##' This function requires that \code{x} already has a valid tref with non-missing
##' \code{origin} and \code{units}. It does not change tref units or period.
##'
##' @examples
##' ## Shift an object to a new explicit origin
##' x <- shift_tref(skjepo$sim$dat, origin = as.Date("1995-01-01"))
##'
##' ## Shift x to use the same origin as y
##' ## x <- shift_tref(x, tref = y)
##'
##' ## Equivalent workflow via add_tref()
##' ## x <- shift_tref(x, tref = tref(skjepo$sim$dat))
##'
##' @seealso \code{\link{add_tref}}, \code{\link{scale_tref}},
##'   \code{\link{create_tref}}, \code{\link{tref}}
##'
##' @export
shift_tref <- function(x, tref = NULL, origin = NULL, verbose = TRUE) {

  tr0 <- try(tref(x), silent = TRUE)
  if (inherits(tr0, "try-error") || is.null(tr0)) {
    stop("Input object 'x' must already have a valid tref.")
  }

  origin0 <- origin(tr0)
  if (.is_na_scalar(origin0)) {
    stop("Cannot shift tref because tref(x)$origin is NA.")
  }

  if (is.null(origin)) {
    tr_target <- .get_tref_target(tref)
    if (is.null(tr_target) || is.null(tr_target$origin)) {
      stop("Provide either 'origin=' or an object/value with tref$origin.")
    }
    origin_new <- tr_target$origin
  } else {
    origin_new <- origin
  }

  if (.is_na_scalar(origin_new)) {
    stop("Target origin is NA.")
  }

  if (.same_origin(origin0, origin_new)) {
    if (isTRUE(verbose)) message("Origin unchanged; nothing to shift.")
    return(x)
  }

  delta <- .tref_origin_offset(origin0, origin_new, tr0)

  ## preserve represented dates:
  ## origin_old + t_old == origin_new + t_new
  ## => t_new = t_old - (origin_new - origin_old)
  x <- .shift_time_values(x, delta = delta)

  tr1 <- tr0
  tr1$origin <- origin_new
  tref(x) <- tr1

  if (isTRUE(verbose)) {
    message("Shifted time values by ", signif(delta, 8), " ", units_time(tr0),
            " to match new tref origin.")
  }

  x
}



##' Compare two time references for equality
##'
##' Test whether two \code{admove_tref} objects describe the same time
##' reference. By default, equality requires matching \code{units} and (when
##' defined) matching \code{period}. Optionally, the \code{origin} can be checked
##' as well.
##'
##' The comparison is designed to be robust to missing values:
##' \itemize{
##'   \item \code{units} must match exactly (character comparison after
##'   normalization).
##'   \item \code{period} is considered equal if both are \code{NA}, or if both
##'   are finite and their absolute difference is \code{<= tol}.
##'   \item If \code{check_origin = TRUE}, \code{origin} is considered equal if
##'   both are \code{NA}, or if both are defined and their absolute time
##'   difference is \code{<= tol} seconds.
##' }
##'
##' @param a,b Objects to compare.
##' @param tol Numeric tolerance used for comparing \code{period} and (if
##'   \code{check_origin = TRUE}) \code{origin}. For \code{origin}, the tolerance
##'   is interpreted in seconds.
##' @param check_origin Logical; if \code{TRUE}, also compare \code{origin}.
##'
##' @return Logical scalar; \code{TRUE} if equal, otherwise \code{FALSE}.
##'
##' @examples
##' tr1 <- create_tref(origin = as.Date("2025-01-01"), units = "month")
##' tr2 <- create_tref(origin = as.Date("2025-01-01"), units = "month", period = 12)
##' tref_equal(tr1, tr2)
##'
##' tr3 <- create_tref(origin = as.Date("2025-01-02"), units = "month")
##' tref_equal(tr1, tr3, check_origin = TRUE)
##'
##' @export
tref_equal <- function(a, b, tol = 1e-12, check_origin = FALSE) {

  if (!inherits(a, "admove_tref") || !inherits(b, "admove_tref")) return(FALSE)

  ## ---- normalize units ------------------------------------------------
  ua <- .normalise_time_unit(a$units)
  ub <- .normalise_time_unit(b$units)
  if (!identical(ua, ub)) return(FALSE)

  ## ---- compare period --------------------------------------------------
  pa <- a$period; if (is.null(pa)) pa <- NA_real_
  pb <- b$period; if (is.null(pb)) pb <- NA_real_

  pa <- suppressWarnings(as.numeric(pa))
  pb <- suppressWarnings(as.numeric(pb))

  if (!(length(pa) == 1 && length(pb) == 1)) return(FALSE)

  if (is.na(pa) && is.na(pb)) {
    ## ok
  } else {
    if (is.na(pa) || is.na(pb) || !is.finite(pa) || !is.finite(pb) || abs(pa - pb) > tol) {
      return(FALSE)
    }
  }

  ## ---- compare origin --------------------------------------------------
  if (isTRUE(check_origin)) {
    oa <- a$origin
    ob <- b$origin

    oa_na <- (length(oa) == 1 && is.na(oa))
    ob_na <- (length(ob) == 1 && is.na(ob))

    if (oa_na && ob_na) {
      ## ok
    } else if (oa_na || ob_na) {
      return(FALSE)
    } else {
      ## POSIXct diffs are in seconds
      diff_abs <- abs(as.numeric(difftime(oa, ob, units = "secs")))
      if (is.na(diff_abs) || diff_abs > tol) return(FALSE)
    }
  }

  TRUE
}


##' Change time units by scaling an admove object
##'
##' Rescale the stored time axis of an \code{admove_*} object and update its
##' time reference (\code{tref(x)}). This is a purely numeric scaling of the
##' stored times and is appropriate when the time axis is represented in
##' fixed-length model units (e.g. "month" meaning 1/12 year, "quarters" meaning
##' 1/4 year, or a custom discretization with \code{tref$period} steps per year).
##'
##' The time reference is updated by:
##' \itemize{
##'   \item multiplying stored times by \code{scale}
##'   \item multiplying \code{tref(x)$period} by \code{scale} (if defined)
##'   \item updating \code{tref(x)$units} (either from \code{units} or by guessing)
##' }
##'
##' @param x An \code{admove_tags}, \code{admove_cov}, \code{admove_data} object,
##'   or a numeric vector/matrix/array of times.
##' @param scale Numeric scalar. Stored times are multiplied by this factor.
##'   For example, to convert years to months (numerically), use \code{scale = 12}.
##' @param units Optional character string giving the new units label to store in
##'   \code{tref(x)$units}. If provided, \code{scale} is inferred from
##'   \code{guess_t_crs_scale(tref(x), units)}.
##' @param verbose Logical; if \code{TRUE}, prints informational messages.
##'
##' @return \code{x} with rescaled stored times and updated \code{tref(x)}.
##'
##' @examples
##' ## Scale an object to a new time step using scale directly (e.g. from month to week)
##' x <- scale_tref(skjepo$sim$dat, scale = 1/12*52)
##'
##' ## Scale by providing units directly
##' x <- scale_tref(skjepo$sim$tags, units = "year")
##'
##' @seealso \code{\link{create_tref}}, \code{\link{add_tref}}, \code{\link{tref}}
##'
##' @export
scale_tref <- function(x, scale = 1, units = NULL, verbose = TRUE) {

  use_units <- !is.null(units)

  if (!use_units && (!is.numeric(scale) || length(scale) != 1 ||
                     !is.finite(scale) || scale <= 0)) {
    stop("'scale' must be a single positive finite number.")
  }

  if (use_units) {
    scale <- .guess_t_crs_scale(tref(x), units)
    if (!is.numeric(scale) || length(scale) != 1 || !is.finite(scale) || scale <= 0) {
      stop("Could not infer a valid positive 'scale' from units='", units, "'.")
    }
  }

  ## ---- rescale stored time -----------------------------------------------
  if (inherits(x, "admove_tags")) {

    x$t <- x$t * scale

  } else if (inherits(x, "admove_cov")) {

    times <- as.numeric(attributes(x)$dimnames[[3]])
    times <- times * scale
    attributes(x)$dimnames[[3]] <- format(times, nsmall = 3)

  } else if (inherits(x, "admove_data")) {

    x$trange <- x$trange * scale
    x$time_cov <- lapply(x$time_cov, function(x) x * scale)
    x$time_spline <- lapply(x$time_spline, function(x) x * scale)
    x$pred$time <- x$pred$time * scale

  } else if (is.numeric(x) || is.matrix(x) || is.array(x)) {

    x <- x * scale

  } else {
    stop("Only works for admove_*tags, admove_cov, admove_data, or numeric/matrix/array.")
  }

  ## ---- update time reference --------------------------------------------
  tr <- tref(x)

  ## update period (cycle length in stored time units)
  if (!is.null(tr$period) && length(tr$period) == 1 && is.finite(tr$period) && !is.na(tr$period)) {
    tr$period <- tr$period * scale
  }

  ## update units label where possible
  if (!use_units) {
    tr$units <- .guess_units_time(tr$units, scale, eps = 1e-10,
                                 days_per_year = 365, verbose = verbose)
  } else {
    tr$units <- units
  }

  tref(x) <- tr
  x
}




## Internal functions -----------------------------------------------------------

.normalise_time_unit <- function(u) {
  if (is.null(u) || length(u) != 1L || is.na(u)) return(NA_character_)
  u0 <- tolower(trimws(u))

  ## Canonical names are singular
  map <- list(
    year = c("years", "year", "y", "yr", "yrs"),
    semester = c("semesters", "semester", "sem", "sems"),
    quarter = c("quarters", "quarter", "qtr", "qtrs"),
    month = c("months", "month", "mo", "mos"),
    week = c("weeks", "week", "wk", "wks"),
    day = c("days", "day", "d"),
    hour = c("hours", "hour", "h", "hr", "hrs"),
    minute = c("minutes", "minute", "min", "mins"),
    second = c("seconds", "second", "sec", "secs", "s"),
    custom = c("custom")
  )

  hit <- vapply(names(map), function(nm) u0 %in% map[[nm]], logical(1))
  if (!any(hit)) return(NA_character_)
  names(map)[which(hit)[1]]
}

.time_unit_fractions <- function(days_per_year = 365.25) {
  ## Fraction of a year represented by one unit
  c(
    year = 1,
    semester = 1/2,
    quarter = 1/4,
    month = 1/12,
    week = 1/52,
    day = 1/days_per_year,
    hour = 1/(days_per_year * 24),
    minute = 1/(days_per_year * 24 * 60),
    second = 1/(days_per_year * 24 * 60 * 60)
  )
}

.guess_t_crs_scale <- function(tref, units, days_per_year = 365.25) {
  ## Returns 'scale' such that: time_new = time_old * scale
  if (is.null(tref)) stop("'tref' must be provided.")
  if (is.null(units) || !length(units)) stop("'units' must be provided.")

  u_target <- .normalise_time_unit(units)
  u_stored <- .normalise_time_unit(tref$units)

  if (is.na(u_target) || is.na(u_stored)) {
    stop("Could not normalise units (stored='", tref$units, "', target='", units, "').")
  }
  if (u_target == "custom" || u_stored == "custom") {
    stop("Cannot infer scale for 'custom' time units. Please provide an explicit numeric scale.")
  }

  frac <- .time_unit_fractions(days_per_year = days_per_year)
  frac_stored <- frac[[u_stored]]
  frac_target <- frac[[u_target]]

  frac_stored / frac_target
}

.guess_units_time <- function(units, scale, eps = 1e-12, days_per_year = 365.25, verbose = TRUE) {
  u <- .normalise_time_unit(units)
  if (is.na(u)) return(units)

  if (is.null(scale) || length(scale) != 1L || is.na(scale)) return(u)
  if (isTRUE(all.equal(scale, 1, tolerance = eps))) return(u)

  if (u == "custom") return("custom")

  frac <- .time_unit_fractions(days_per_year = days_per_year)
  ratios <- frac[[u]] / frac
  diffs <- abs(ratios - scale)
  j <- which.min(diffs)

  if (isTRUE(all.equal(unname(scale), unname(ratios[j]), tolerance = eps))) {
    return(names(frac)[j])
  }

  ## Fallback: keep info but don't pretend it's a known unit
  paste0(u, "_x_", format(scale, scientific = TRUE))
}

.add_rel_time_posix <- function(origin, x, units, days_per_year = 365.25) {
  stopifnot(inherits(origin, "POSIXt"))

  u <- .normalise_time_unit(units)
  if (is.na(u) || u == "custom") stop("Unsupported time unit: ", units)

  frac <- .time_unit_fractions(days_per_year = days_per_year)
  sec_per_unit <- frac[[u]] * days_per_year * 24 * 60 * 60

  origin + as.difftime(x * sec_per_unit, units = "secs")
}

.infer_units_from_period <- function(period, eps = 1e-12) {
  ## Heuristic only for the classic annual discretizations
  if (length(period) != 1 || is.na(period)) return(NA_character_)
  if (isTRUE(all.equal(period, 1,  tolerance = eps))) return("year")
  if (isTRUE(all.equal(period, 2,  tolerance = eps))) return("semester")
  if (isTRUE(all.equal(period, 4,  tolerance = eps))) return("quarter")
  if (isTRUE(all.equal(period, 12, tolerance = eps))) return("month")
  if (isTRUE(all.equal(period, 52, tolerance = eps))) return("week")
  "custom"
}

.default_period_from_units <- function(units) {
  u <- .normalise_time_unit(units)
  switch(u,
         year = 1,
         semester = 2,
         quarter = 4,
         month = 12,
         week = 52,
         ## No obvious defaults for these
         day = NA_real_,
         hour = NA_real_,
         minute = NA_real_,
         second = NA_real_,
         custom = NA_real_,
         NA_real_)
}

.same_origin <- function(x, y) {
  if (.is_na_scalar(x) && .is_na_scalar(y)) return(TRUE)
  if (.is_na_scalar(x) || .is_na_scalar(y)) return(FALSE)
  isTRUE(all.equal(x, y))
}

.get_tref_target <- function(tref) {

  if (is.null(tref)) return(NULL)

  if (inherits(tref, "admove_tref") ||
      any(names(tref) %in% c("origin", "units", "period"))) {
    return(tref)
  }

  tr <- try(tref(tref), silent = TRUE)
  if (inherits(tr, "try-error")) return(NULL)

  tr
}

.tref_origin_offset <- function(origin_from, origin_to, tr) {

  u <- units_time(tr)

  if (length(u) != 1L || is.na(u)) {
    stop("Cannot shift tref origin because tref$units are missing.")
  }

  if (.is_na_scalar(origin_from) || .is_na_scalar(origin_to)) {
    stop("Cannot shift tref origin because one of the origins is NA.")
  }

  ## harmonise classes
  if (inherits(origin_from, "POSIXt") || inherits(origin_to, "POSIXt")) {
    tz <- attr(origin_to, "tzone", exact = TRUE) %||%
      attr(origin_from, "tzone", exact = TRUE) %||%
      "UTC"
    origin_from <- as.POSIXct(origin_from, tz = tz)
    origin_to <- as.POSIXct(origin_to,   tz = tz)
  } else {
    origin_from <- as.Date(origin_from)
    origin_to <- as.Date(origin_to)
  }

  if (u %in% c("second", "minute", "hour", "day", "week")) {
    return(
      switch(u,
             second = as.numeric(difftime(origin_to, origin_from, units = "secs")),
             minute = as.numeric(difftime(origin_to, origin_from, units = "mins")),
             hour = as.numeric(difftime(origin_to, origin_from, units = "hours")),
             day = as.numeric(difftime(origin_to, origin_from, units = "days")),
             week = as.numeric(difftime(origin_to, origin_from, units = "days")) / 7)
    )
  }

if (u %in% c("month", "year")) {
  return(
    lubridate::time_length(
      lubridate::interval(origin_from, origin_to),
      unit = u
    )
  )
}

if (u == "quarter") {
  return(
    lubridate::time_length(
      lubridate::interval(origin_from, origin_to),
      unit = "month"
    ) / 3
  )
}

  stop("Unsupported tref units: ", u)
}

.shift_time_values <- function(x, delta) {

  if (inherits(x, "admove_tags")) {

    x$t <- x$t - delta

  } else if (inherits(x, "admove_cov")) {

    times <- as.numeric(attributes(x)$dimnames[[3]])
    times <- times - delta
    attributes(x)$dimnames[[3]] <- format(times, trim = TRUE, scientific = FALSE)

  } else if (inherits(x, "admove_data")) {

    x$trange <- x$trange - delta
    x$time_cov <- lapply(x$time_cov, function(x) x - delta)
    x$time_spline <- lapply(x$time_spline, function(x) x - delta)
    x$pred$time <- x$pred$time - delta

  } else if (is.numeric(x) || is.matrix(x) || is.array(x)) {

    x <- x - delta

  } else {
    stop("Only works for admove_*tags, admove_cov, admove_data, or numeric/matrix/array.")
  }

  x
}


.format_tref_short <- function(x, labw = 15) {
  stopifnot(inherits(x, "admove_tref"))

  ## Format origin nicely
  origin_txt <- NULL
  if (!is.null(x$origin) && length(x$origin) == 1L && !is.na(x$origin)) {
    if (inherits(x$origin, "POSIXt")) {
      origin_txt <- format(x$origin, tz = "UTC", usetz = TRUE)
    } else if (inherits(x$origin, "Date")) {
      origin_txt <- format(x$origin)
    } else {
      origin_txt <- as.character(x$origin)
    }
  }

  ## Format units
  units_txt <- NULL
  if (!is.null(x$units) && length(x$units) == 1L && !is.na(x$units)) {
    units_txt <- as.character(x$units)
  }

  ## Format period
  period_txt <- NULL
  if (!is.null(x$period) && length(x$period) == 1L && !is.na(x$period)) {
    period_txt <- as.character(x$period)
  }

  out <- c(
    if (!is.null(origin_txt))
      sprintf(paste0("  %-", labw, "s %s"), "origin:", origin_txt),
    if (!is.null(units_txt))
      sprintf(paste0("  %-", labw, "s %s"), "units:", units_txt),
    if (!is.null(period_txt))
      sprintf(paste0("  %-", labw, "s %s"), "period:", period_txt)
  )

  out
}
