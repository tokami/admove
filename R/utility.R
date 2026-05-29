
##' @importFrom grDevices adjustcolor col2rgb grey hcl.colors n2mfrow rgb terrain.colors
##' @importFrom graphics abline arrows axis box contour identify image layout legend lines mtext par plot.new points polygon segments text
##' @importFrom stats dist median qnorm quantile rnorm runif setNames
##' @importFrom utils capture.output packageDescription tail
##' @importFrom RTMB ADoverload REPORT ADREPORT
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

  ## helper to coerce to same broad class as 'dates'
  if (inherits(dates, "POSIXt")) {
    tz <- attr(dates, "tzone", exact = TRUE)
    if (is.null(tz) || !length(tz) || identical(tz, "")) tz <- "UTC"

    coerce_time <- function(x) {
      as.POSIXct(x, tz = tz)
    }
  } else {
    coerce_time <- function(x) {
      as.Date(x)
    }
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
    origin <- coerce_time(tref)

  } else if (is.list(tref)) {

    inferred <- isTRUE(tref$inferred)
    units <- tref$units %||% guess_units(dates)
    floor_unit <- tref$floor_unit %||% guess_floor_unit(units)

    origin <- tref$origin %||% tref$t0 %||% tref$start
    if (is.null(origin)) {
      stop("If 'tref' is a list, it must contain an element named 'origin'.")
    }
    origin <- coerce_time(origin)

  } else {
    stop("'tref' must be NULL, a Date/POSIXt object, or a list with element 'origin'.")
  }

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
    stop("Unsupported units: ", units)
  }

  attr(out, "tref") <- list(
    origin = origin,
    units = units,
    floor_unit = floor_unit,
    inferred = inferred
  )

  out
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
    ## start with all observed times; add interior points in big gaps only
    ts <- t_obs

  if (!is.finite(dt_min) || dt_min <= 0) stop("dt_min must be > 0.")

    gaps <- diff(t_obs)
    if (any(gaps > dt_min + eps)) {
      extra <- unlist(lapply(seq_along(gaps), function(i) {
        if (gaps[i] <= dt_min + eps) return(numeric(0))
        ## insert points: t_obs[i] + dt_min, ..., strictly before t_obs[i+1]
        if ((t_obs[i] + dt_min) < (t_obs[i + 1] - dt_min)) {
          seq(t_obs[i] + dt_min, t_obs[i + 1] - dt_min, by = dt_min)
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

calc_mstar <- function(fit) {

  nc <- nrow(fit$dat$pred$grid$xygrid)
  nt <- length(fit$dat$pred$time) - 1
  dts <- diff(fit$dat$pred$time)
  cs <- diff(sort(unique(fit$dat$pred$grid$xygrid[,1])))[1]
  nextTo <- get_neighbours(fit$dat$pred$grid)
  mstar <- array(0, c(nc,nc,nt))
  for(t in 1:nt){
    for(i in 1:nc){
      dif <- exp(fit$pred$hD[i,t])
      if(fit$conf$use_taxis){
        tax <- dif * c(fit$pred$hTdx[i,t], fit$pred$hTdy[i,t])
      }else{
        tax <- c(0,0)
      }
      if(!is.na(nextTo[i,2])){
        mstar[i,nextTo[i,2],t] <- (dif/cs/cs + 0.5 * tax[2]/cs) *  dts[t]
      }
      if(!is.na(nextTo[i,3])){
        mstar[i,nextTo[i,3],t] <- (dif/cs/cs - 0.5 * tax[2]/cs) *  dts[t]
      }
      if(!is.na(nextTo[i,4])){
        mstar[i,nextTo[i,4],t] <- (dif/cs/cs - 0.5 * tax[1]/cs) *  dts[t]
      }
      if(!is.na(nextTo[i,5])){
        mstar[i,nextTo[i,5],t] <- (dif/cs/cs + 0.5 * tax[1]/cs) *  dts[t]
      }
    }
    if(fit$conf$use_advection){
      for(i in 1:nc){
        adv <- c(fit$pred$hAx[i,t], fit$pred$hAy[i,t])
        if(!is.na(nextTo[i,2])){
          mstar[i,nextTo[i,2],t] <- mstar[i,nextTo[i,2],t] + 0.5 * adv[2]/cs *  dts[t]
        }
        if(!is.na(nextTo[i,3])){
          mstar[i,nextTo[i,3],t] <- mstar[i,nextTo[i,3],t] - 0.5 * adv[2]/cs *  dts[t]
        }
        if(!is.na(nextTo[i,4])){
          mstar[i,nextTo[i,4],t] <- mstar[i,nextTo[i,4],t] - 0.5 * adv[1]/cs *  dts[t]
        }
        if(!is.na(nextTo[i,5])){
          mstar[i,nextTo[i,5],t] <- mstar[i,nextTo[i,5],t] + 0.5 * adv[1]/cs *  dts[t]
        }
      }
    }
    diag(mstar[,,t]) <- -rowSums(mstar[,,t])
  }
  return(mstar)
}


fill_inst_mat <- function(mat, move, nextTo, next_dist) {
  xyind <- c(2, 2, 1, 1)
  dirsign <- c(+1, -1, -1, +1)
  ## AD-safe pmax
  pos <- function(x) 0.5 * (x + abs(x))
  ## 4 neighbours
  for (k in 1:4) {
    j <- k + 1
    ind <- which(!is.na(nextTo[, j]))
    v <- dirsign[k] * move[ind, xyind[k]]
    mat[cbind(ind, nextTo[ind, j])] <- pos(v) / next_dist[k]
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
##' centres of the supplied object.
##'
##' @param x An object providing spatial dimensions and cell centres, typically
##'   an `admove_grid` or another object for which [x_centers()], [y_centers()],
##'   and [sref()] are defined.
##' @param tref Optional time reference information to attach to the returned
##'   covariates. Default: `NULL`.
##'
##' @return
##' A list-like object of class `admove_cov` containing two covariate fields:
##' one for the x coordinate and one for the y coordinate.
##'
##' @details
##' Both covariates are created as single-time-slice fields with `times = 0`.
##' Spatial reference information is copied from `x`, and optional temporal
##' reference information can be attached via `tref`.
##'
##' @examples
##' ## xy_cov <- make_x_y_cov(grid)
##'
##' @export
make_x_y_cov <- function(x, tref = NULL) {

  nx <- dim(x)[1]
  ny <- dim(x)[2]

  cov1 <- prep_cov(matrix(1:nx, nx, ny),
                   ## matrix(seq(0.5,nx,1) - nx/2, nx, ny), ## not working
                   x_centers = x_centers(x),
                   y_centers = y_centers(x),
                   times = 0,
                   sref = sref(x),
                   tref = tref)
  cov2 <- prep_cov(matrix(1:ny, nx, ny, byrow = TRUE),
                   ## matrix(seq(0.5,ny,1) - ny/2, nx, ny, byrow = TRUE),
                   x_centers = x_centers(x),
                   y_centers = y_centers(x),
                   times = 0,
                   sref = sref(x),
                   tref = tref)
  cov <- list(cov1, cov2)
  cov <- .add_class(cov, "admove_cov")
  cov <- add_sref(cov, sref(x))
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



.poly_fun <- function(xp, yp, deriv = FALSE, adv = FALSE) {

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
    } else {
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
    }
  }

  if (!deriv) return(f) else return(df)
}



.make_pref_funcs <- function(alpha, beta, gamma,
                            knots_tax, knots_dif) {

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
      pref_funcs$tax[[i]][[j]] <- .poly_fun(knots_tax[,i], alpha[,i,j])
      pref_funcs$dtax[[i]][[j]] <- .poly_fun(knots_tax[,i],
                                          alpha[,i,j],
                                          deriv = TRUE)
    }

    ## diffusion
    pref_funcs$dif[[i]] <- pref_funcs$ddif[[i]] <- vector("list", dim(beta)[3])
    for(j in 1:dim(beta)[3]){
      pref_funcs$dif[[i]][[j]] <- .poly_fun(knots_dif[,i], beta[,i,j])
      pref_funcs$ddif[[i]][[j]] <- .poly_fun(knots_dif[,i], beta[,i,j],
                                          deriv=TRUE)
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
