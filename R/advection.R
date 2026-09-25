##' Prepare an advection field
##'
##' @description
##' `prep_adv()` prepares a vector field that moves the tagged animals
##' passively, typically ocean currents (or wind), from its two components: `u`
##' towards geographic east and `v` towards geographic north. The result is
##' passed to [setup_data()] through its `adv` argument, separately from the
##' habitat covariates (`cov`) that inform taxis and diffusion.
##'
##' @param u,v The east and north components of the field. Either `admove_cov`
##'   objects (see [prep_cov()]) or any input [prep_cov()] accepts, in which
##'   case `...` is passed on to it. Both components must share the grid, the
##'   time steps and the spatial and temporal reference.
##' @param units Physical units of `u` and `v`: one of `"m/s"` (default),
##'   `"cm/s"`, `"mm/s"`, `"km/h"`, `"km/day"` or `"knots"`. [setup_data()]
##'   converts the field into model units (space units of the grid per time
##'   unit of the model) and onto the axes of the grid, which needs a
##'   coordinate reference system on the field. `NULL` means that `u` and `v`
##'   are already in model units and already along the grid's x and y axes, so
##'   no conversion is done.
##' @param verbose Logical; passed to [prep_cov()].
##' @param ... Further arguments passed to [prep_cov()] when `u` and `v` are
##'   not yet `admove_cov` objects.
##'
##' @details
##' Advection adds `gamma * (u, v)` to the movement of a tag, where `gamma` is
##' an estimated entrainment coefficient (see [default_par()] and
##' [default_map()]). After the conversion in [setup_data()] `gamma` is
##' dimensionless: `gamma = 1` is a passive drifter moving with the field, and
##' `gamma = 0.1` means the animal is carried at a tenth of its speed.
##'
##' The conversion differentiates the map projection at every cell of the field
##' (the local Jacobian), so it handles longitude/latitude grids, where a degree
##' of longitude shrinks with the cosine of latitude, as well as projected
##' grids, where grid north is rotated away from true north and distances are
##' stretched away from the projection centre. `u` and `v` must therefore be
##' geographic east and north components. Output of ocean models on their native
##' curvilinear grid is often grid-relative and has to be rotated to east/north
##' first.
##'
##' Several fields (e.g. currents and wind) can be passed to [setup_data()] as a
##' named list of `prep_adv()` objects; each gets its own `gamma`.
##'
##' Taxis and advection both explain directed movement, and from the default
##' start (`alpha` and `gamma` at 0) the optimizer can end in a local optimum
##' where the taxis takes the drift of the field. If the fitted taxis looks
##' implausible, fit `gamma` first with the taxis coefficients fixed
##' (`map$alpha <- factor(rep(NA, length(par$alpha)))`), refit all parameters
##' from those estimates (`par = fit$pl`) and compare the objective values.
##'
##' @return An object of class `admove_adv`: a list with the `admove_cov`
##'   components `u` and `v`, and the physical units in `attr(, "units")`.
##'
##' @seealso [setup_data()], [plot_advection()]
##'
##' @examples
##' \dontrun{
##' cur <- prep_adv(cur_u, cur_v, units = "m/s")
##' dat <- setup_data(grid, cov = list(temp = temp), tags = tags, adv = cur)
##' }
##'
##' @export
prep_adv <- function(u, v, units = "m/s", verbose = TRUE, ...) {

  u <- .as_adv_component(u, "u", verbose, ...)
  v <- .as_adv_component(v, "v", verbose, ...)

  if (!identical(dim(u), dim(v)) ||
        !isTRUE(all.equal(lapply(dimnames(u), as.numeric),
                          lapply(dimnames(v), as.numeric)))) {
    stop("'u' and 'v' must share the grid and the time steps (same dimensions ",
         "and dimnames).", call. = FALSE)
  }
  if (!sref_equal(sref(u), sref(v))) {
    stop("'u' and 'v' have different spatial references.", call. = FALSE)
  }
  if (!tref_equal(tref(u), tref(v))) {
    stop("'u' and 'v' have different time references.", call. = FALSE)
  }
  if (all(is.na(u)) || all(is.na(v))) {
    stop("All values of '", if (all(is.na(u))) "u" else "v", "' are NA.",
         call. = FALSE)
  }

  if (!is.null(units)) .adv_speed_factor(units)

  structure(list(u = u, v = v), class = "admove_adv", units = units)
}


.as_adv_component <- function(x, nm, verbose, ...) {
  if (inherits(x, "admove_cov_list")) {
    if (length(x) != 1L) {
      stop("'", nm, "' must be a single field, not a list of ", length(x),
           ".", call. = FALSE)
    }
    x <- x[[1L]]
  }
  if (!inherits(x, "admove_cov")) {
    x <- prep_cov(x, verbose = verbose, ...)
    if (inherits(x, "admove_cov_list")) x <- x[[1L]]
  }
  x
}


##' @export
print.admove_adv <- function(x, ...) {
  u <- x$u
  times <- as.numeric(dimnames(u)[[3]])
  spd <- sqrt(unclass(x$u)^2 + unclass(x$v)^2)
  un <- attr(x, "units")
  cat("<admove_adv> advection field\n")
  cat("  grid:     ", dim(u)[1], "x", dim(u)[2], "cells\n")
  cat("  times:    ", length(times), "step(s)",
      if (length(times)) paste0("[", signif(min(times), 6), ", ",
                                signif(max(times), 6), "]"), "\n")
  cat("  units:    ",
      if (is.null(un)) "model units (space unit per time unit, grid axes)"
      else paste0(un, " (east/north)"), "\n")
  cat("  speed:    ",
      paste(signif(range(spd, na.rm = TRUE), 4), collapse = " - "), "\n")
  invisible(x)
}


##' Plot advection fields
##'
##' @description
##' `plot_adv_field()` plots an advection field as it enters the model: the
##' speed of the field as a colour image with arrows for its direction, one panel
##' per time step. It shows the *input* field (see [prep_adv()]); the advection
##' *estimated* by a fit (the field scaled by `gamma`, plus a constant drift) is
##' drawn by [plot_advection()].
##'
##' @param x An advection field created by [prep_adv()], a list of them, or an
##'   object carrying them: `admove_data` (see [setup_data()]), `admove_sim` or
##'   `admove`.
##' @param i The field to plot, by index or name. Default: `1`.
##' @param select Optional time-step indices of the field. Default: `NULL`, all
##'   time steps.
##' @param main Main title, drawn above the panels together with the units of
##'   the speed. Default: `"Advection field"`.
##' @param plot_land Logical; if `TRUE`, add land. Default: `FALSE`.
##' @param auto_layout Logical; if `TRUE` (default), the panel layout is set
##'   and restored afterwards. With several panels the axes are drawn on the
##'   outer panels only.
##' @param n_arrows Approximate number of arrows along each axis. Default: `25`.
##' @param cor Optional scaling factor for the arrow lengths. If `NULL`
##'   (default), the fastest arrow is as long as the spacing between arrows. The
##'   same scale is used in all panels.
##' @param col Colour palette for the speed.
##' @param col_arrows Colour of the arrows. Default: `"black"`.
##' @param zlim Optional speed range mapped onto `col`. Default: `NULL`, the
##'   range over all plotted time steps, so that the panels share one scale.
##' @param legend Logical; if `TRUE`, a colour bar is drawn to the right of each
##'   panel. Default: `NULL`, `TRUE` when `auto_layout = TRUE`.
##' @param titles Optional panel titles, one per time step. Default: `NULL`, the
##'   start time of each time step as a date (e.g. `"Jan 2007"`); a single panel
##'   gets no title.
##' @param xlab,ylab Axis labels. If `NULL` (default), `"x"` and `"y"` with
##'   the spatial units in brackets.
##' @param land_col,land_border Fill and border colour of the land.
##' @param ... Further graphical arguments passed to [plot()].
##'
##' @return Invisibly `NULL`.
##'
##' @seealso [prep_adv()], [plot_advection()], [plot_cov()]
##'
##' @examples
##' \dontrun{
##' cur <- prep_adv(cur_u, cur_v, units = "m/s")
##' plot_adv_field(cur, select = 1:4, plot_land = TRUE)
##' dat <- setup_data(grid, cov = cov, tags = tags, adv = cur)
##' plot_adv_field(dat, select = 1:4)   ## in model units (e.g. km/month)
##' }
##'
##' @export
plot_adv_field <- function(x,
                           i = 1,
                           select = NULL,
                           main = "Advection field",
                           plot_land = FALSE,
                           auto_layout = TRUE,
                           n_arrows = 25,
                           cor = NULL,
                           col = hcl.colors(100, "YlOrRd", rev = TRUE),
                           col_arrows = "black",
                           zlim = NULL,
                           legend = NULL,
                           titles = NULL,
                           xlab = NULL,
                           ylab = NULL,
                           land_col = grey(0.85),
                           land_border = grey(0.5),
                           ...) {

  adv <- if (inherits(x, "admove_adv")) {
    list(x)
  } else if (inherits(x, c("admove", "admove_sim"))) {
    x$dat$adv
  } else if (inherits(x, "admove_data")) {
    x$adv
  } else if (is.list(x) && length(x) > 0L &&
               all(vapply(x, inherits, logical(1), "admove_adv"))) {
    x
  } else NULL
  if (length(adv) == 0L) stop("'x' has no advection field.", call. = FALSE)

  if (is.character(i)) {
    if (!i %in% names(adv)) {
      stop("No advection field '", i, "'. Available: ",
           paste(names(adv), collapse = ", "), call. = FALSE)
    }
  } else if (length(i) != 1L || i < 1 || i > length(adv)) {
    stop("'i' must select one of the ", length(adv), " advection field(s).",
         call. = FALSE)
  }
  fld <- adv[[i]]
  u <- fld$u
  v <- fld$v

  nt_all <- dim(u)[3]
  if (is.null(select)) select <- seq_len(nt_all)
  if (any(select < 1) || any(select > nt_all)) {
    stop("'select' must be time steps of the field (1 to ", nt_all, ").",
         call. = FALSE)
  }
  nt <- length(select)

  map_labs <- .map_labs(u)
  if (is.null(xlab)) xlab <- map_labs[1]
  if (is.null(ylab)) ylab <- map_labs[2]
  if (is.null(legend)) legend <- isTRUE(auto_layout)

  ## units of the speed: physical before setup_data(), model units after
  un <- attr(fld, "units")
  if (is.null(un)) {
    us <- tryCatch(units_space(u), error = function(e) NA)
    ut <- tryCatch(.normalise_time_unit(units_time(u)), error = function(e) NA)
    un <- if (length(us) == 1L && length(ut) == 1L && !is.na(us) && !is.na(ut))
      paste0(us, "/", ut) else NULL
  }
  main_txt <- if (nzchar(main) && !is.null(un)) paste0(main, " [", un, "]") else main

  xs <- as.numeric(dimnames(u)[[1]])
  ys <- as.numeric(dimnames(u)[[2]])
  au <- unclass(u)[, , select, drop = FALSE]
  av <- unclass(v)[, , select, drop = FALSE]
  spd <- sqrt(au^2 + av^2)

  if (is.null(titles)) {
    titles <- if (nt > 1L) {
      .time_labels(as.numeric(dimnames(u)[[3]])[select], tref(u))
    } else ""
  } else if (length(titles) != nt) {
    stop("'titles' must have one entry per plotted time step (", nt, ").",
         call. = FALSE)
  }

  ## one colour scale and one arrow scale for all panels
  if (is.null(zlim)) {
    zlim <- suppressWarnings(range(spd, na.rm = TRUE, finite = TRUE))
  }
  zlim_ok <- all(is.finite(zlim))
  if (zlim_ok && zlim[1] == zlim[2]) zlim <- zlim + c(-0.5, 0.5) * max(abs(zlim[1]), 1)

  ix <- unique(round(seq(1, length(xs), length.out = min(n_arrows, length(xs)))))
  iy <- unique(round(seq(1, length(ys), length.out = min(n_arrows, length(ys)))))
  g <- expand.grid(i = ix, j = iy)
  spacing <- min(diff(range(xs)) / max(1, length(ix) - 1),
                 diff(range(ys)) / max(1, length(iy) - 1))
  if (is.null(cor)) {
    smax <- suppressWarnings(max(spd[cbind(rep(g$i, nt), rep(g$j, nt),
                                           rep(seq_len(nt), each = nrow(g)))],
                                 na.rm = TRUE))
    cor <- if (is.finite(smax) && smax > 0) spacing / smax else 1
  }

  dots <- list(...)
  xaxt <- if (is.null(dots$xaxt)) "s" else dots$xaxt
  yaxt <- if (is.null(dots$yaxt)) "s" else dots$yaxt
  dots$xaxt <- dots$yaxt <- NULL
  shared <- auto_layout && nt > 1L

  if (auto_layout) {
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar))
    mar_right <- if (legend) 4.5 else 1.5
    par(mfrow = n2mfrow(nt, asp = 2),
        mar = c(if (shared) 0.3 else 1.5, if (shared) 0.3 else 1.5,
                if (shared) 1.4 else 1.5, mar_right),
        oma = c(3, 3.5, if (nzchar(main_txt)) 1.5 else 0, 0),
        mgp = c(2, 0.5, 0),
        tcl = -0.3)
  }

  for (k in seq_len(nt)) {
    do.call(plot, c(list(1, 1, type = "n", xlim = range(xs), ylim = range(ys),
                         xlab = if (auto_layout) "" else xlab,
                         ylab = if (auto_layout) "" else ylab,
                         asp = 1,
                         xaxt = if (shared) "n" else xaxt,
                         yaxt = if (shared) "n" else yaxt),
                    dots))
    if (shared) .outer_panel_axes(k, nt, xaxt, yaxt)
    if (zlim_ok) {
      z <- pmin(pmax(spd[, , k], zlim[1]), zlim[2])
      image(xs, ys, z, col = col, zlim = zlim, add = TRUE)
    }
    if (plot_land) plot_land(sref(u), col = land_col, border = land_border)
    du <- au[cbind(g$i, g$j, k)]
    dv <- av[cbind(g$i, g$j, k)]
    ok <- is.finite(du) & is.finite(dv) & (du != 0 | dv != 0)
    if (any(ok)) {
      ## arrows too short to draw at this scale are skipped; drop the warning
      withCallingHandlers(
        arrows(xs[g$i][ok], ys[g$j][ok],
               xs[g$i][ok] + du[ok] * cor, ys[g$j][ok] + dv[ok] * cor,
               length = 0.04, col = col_arrows),
        warning = function(w) {
          if (grepl("zero-length arrow", conditionMessage(w)))
            invokeRestart("muffleWarning")
        })
    }
    if (nzchar(titles[k])) {
      title(main = titles[k], line = if (shared) 0.3 else 0.4,
            font.main = 1, cex.main = 1)
    }
    box(lwd = 1.5)
    if (legend && zlim_ok) .color_bar(col, zlim)
  }

  if (auto_layout) {
    mtext(main_txt, 3, 0, outer = TRUE)
    mtext(xlab, 1, if (shared) 2 else 1, outer = TRUE)
    mtext(ylab, 2, if (shared) 2.2 else 1.5, outer = TRUE)
  }

  invisible(NULL)
}


##' @rdname plot_adv_field
##' @export
plot.admove_adv <- function(x, ...) {
  plot_adv_field(x, ...)
}


## Metres per second per unit of the given speed units.
.adv_speed_factor <- function(units) {
  if (!is.character(units) || length(units) != 1L || is.na(units)) {
    stop("'units' must be a single character string such as \"m/s\".",
         call. = FALSE)
  }
  key <- gsub("\\s+", "", tolower(units))
  key <- sub("^(m|cm|mm|km)(s-1|s\\^-1)$", "\\1/s", key)
  key <- sub("^kmh-1$", "km/h", key)
  fac <- c("m/s" = 1, "cm/s" = 0.01, "mm/s" = 0.001,
           "km/h" = 1000 / 3600, "km/day" = 1000 / 86400, "km/d" = 1000 / 86400,
           "knots" = 1852 / 3600, "knot" = 1852 / 3600, "kn" = 1852 / 3600,
           "kt" = 1852 / 3600)
  if (!key %in% names(fac)) {
    stop("Unknown speed units \"", units, "\". Use one of \"m/s\", \"cm/s\", ",
         "\"mm/s\", \"km/h\", \"km/day\" or \"knots\", or units = NULL for a ",
         "field already in model units.", call. = FALSE)
  }
  unname(fac[[key]])
}


## Normalise the 'adv' input of setup_data() to a named list of admove_adv.
.as_adv_list <- function(adv) {
  if (is.null(adv)) return(NULL)
  if (inherits(adv, "admove_adv")) {
    adv <- list(adv = adv)
  } else if (!is.list(adv) || length(adv) == 0L ||
               !all(vapply(adv, inherits, logical(1), "admove_adv"))) {
    stop("'adv' must be a field created by prep_adv(), or a list of them.",
         call. = FALSE)
  }
  nms <- names(adv)
  if (is.null(nms)) nms <- rep("", length(adv))
  nms[!nzchar(nms)] <- paste0("adv", which(!nzchar(nms)))
  if (anyDuplicated(nms)) stop("The advection fields must have unique names.",
                               call. = FALSE)
  names(adv) <- nms
  adv
}


## The components of all fields as one covariate list, in the order u, v of the
## first field, u, v of the second, ... Everything that indexes the components
## (dat$time_adv, dat$xrange_adv, .make_adv()) relies on this order.
.adv_flatten <- function(adv) {
  if (length(adv) == 0L) return(NULL)
  flat <- unlist(lapply(adv, function(a) list(a$u, a$v)), recursive = FALSE)
  names(flat) <- paste0(rep(names(adv), each = 2L), c("_u", "_v"))
  .add_class(flat, "admove_cov_list")
}


## Inverse of .adv_flatten(), keeping the units of `template`.
.adv_unflatten <- function(flat, template) {
  out <- lapply(seq_along(template), function(f) {
    structure(list(u = flat[[2L * f - 1L]], v = flat[[2L * f]]),
              class = "admove_adv", units = attr(template[[f]], "units"))
  })
  names(out) <- names(template)
  out
}


## Convert every field that carries physical units into model units along the
## grid axes: space units of its sref per time unit of its tref. Called by
## setup_data() after the spatial and temporal references are harmonised.
.adv_to_model_units <- function(adv, verbose = TRUE) {
  for (f in seq_along(adv)) {
    un <- attr(adv[[f]], "units")
    if (is.null(un)) next
    adv[[f]] <- .adv_convert_field(adv[[f]], un, names(adv)[f], verbose)
    attr(adv[[f]], "units") <- NULL
  }
  adv
}


## Geographic (east, north) speeds to velocities along the grid axes. Each cell
## centre is inverse-projected to lon/lat and the projection is differentiated
## there, so the image of a unit east and a unit north velocity is exact: this
## covers the cos(latitude) shrinking of a degree of longitude, the rotation of
## grid north away from true north and the scale distortion of projected grids
## in one expression. See dev/code_notes.org, "Advection fields and
## entrainment".
.adv_convert_field <- function(a, units, name = "adv", verbose = TRUE) {

  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Converting the advection field from ", units, " needs the 'sf' ",
         "package. Install it, or convert the field yourself and use ",
         "prep_adv(..., units = NULL).", call. = FALSE)
  }

  u <- a$u
  v <- a$v
  sr <- sref(u)
  crs_sf <- tryCatch(sf::st_crs(sr$crs), error = function(e) NA)
  if (is.null(sr$crs) || .is_na_scalar(sr$crs) || is.na(crs_sf)) {
    stop("The advection field '", name, "' has no coordinate reference system, ",
         "so it cannot be converted from ", units, " to grid units. Give it an ",
         "sref with a crs (prep_cov(..., sref = )), or convert it yourself and ",
         "use prep_adv(..., units = NULL).", call. = FALSE)
  }
  cs <- sr$crs_scale
  if (is.null(cs) || .is_na_scalar(cs)) cs <- 1

  tu <- .normalise_time_unit(units_time(tref(u)))
  if (is.na(tu) || tu == "custom") {
    stop("The advection field '", name, "' has no standard time unit, so ",
         "speeds in ", units, " cannot be converted to model time. Give it a ",
         "tref with units (e.g. \"month\"), or use prep_adv(..., units = NULL).",
         call. = FALSE)
  }
  frac <- .time_unit_fractions()
  s_per_t <- frac[[tu]] / frac[["second"]]

  ## metres per model time unit, per unit of the input
  fac <- .adv_speed_factor(units) * s_per_t

  xs <- as.numeric(dimnames(u)[[1]])
  ys <- as.numeric(dimnames(u)[[2]])
  xy <- as.matrix(expand.grid(x = xs, y = ys))
  crs_wkt <- crs_sf$wkt
  ll <- sf::sf_project(crs_wkt, "OGC:CRS84", xy / cs)
  lon <- ll[, 1]
  lat <- ll[, 2]

  prj <- function(lo, la) sf::sf_project("OGC:CRS84", crs_wkt, cbind(lo, la),
                                         warn = FALSE) * cs
  dd <- 1e-3                                      ## degrees
  d_rad <- dd * pi / 180
  ## grid displacement per radian of longitude / latitude, in stored units
  j_lon <- (prj(lon + dd, lat) - prj(lon - dd, lat)) / (2 * d_rad)
  j_lat <- (prj(lon, lat + dd) - prj(lon, lat - dd)) / (2 * d_rad)
  ## Radii of curvature of the WGS84 ellipsoid (the lon/lat above are CRS84):
  ## metres per radian of latitude (m_rad) and, times cos(lat), of longitude
  ## (n_rad). The mean Earth radius instead would be off by up to ~0.7%.
  a_wgs <- 6378137
  e2 <- 0.00669437999014
  sin2 <- sin(lat * pi / 180)^2
  n_rad <- a_wgs / sqrt(1 - e2 * sin2)
  m_rad <- a_wgs * (1 - e2) / (1 - e2 * sin2)^1.5
  coslat <- cos(lat * pi / 180)

  au <- unclass(u)
  av <- unclass(v)
  ou <- au
  ov <- av
  for (t in seq_len(dim(au)[3])) {
    dlon <- as.vector(au[, , t]) * fac / (n_rad * coslat)  ## radians / time
    dlat <- as.vector(av[, , t]) * fac / m_rad
    ou[, , t] <- j_lon[, 1] * dlon + j_lat[, 1] * dlat
    ov[, , t] <- j_lon[, 2] * dlon + j_lat[, 2] * dlat
  }
  ou[!is.finite(ou)] <- NA
  ov[!is.finite(ov)] <- NA

  if (verbose) {
    rot <- atan2(j_lat[, 1], j_lat[, 2]) * 180 / pi
    ok <- is.finite(rot)
    su <- units_space(sr)
    message("Advection field '", name, "' converted from ", units, " to ",
            if (!is.null(su) && !.is_na_scalar(su)) su else "space units",
            "/", tu,
            if (any(ok) && max(abs(rot[ok])) > 0.05)
              sprintf("; grid north deviates from true north by %.1f to %.1f deg",
                      min(rot[ok]), max(rot[ok])),
            ".")
  }

  u[] <- ou
  v[] <- ov
  a$u <- u
  a$v <- v
  a
}


## Season of the advection coefficients at time t (1 when not seasonal). Seasons
## are equal divisions of the cycle, anchored at the time origin, as for the
## habitat splines (.season_breaks()).
.adv_season_fun <- function(nsea, period) {
  if (nsea <= 1L) return(function(t) 1L)
  if (is.null(period) || !is.finite(period) || period <= 0) {
    stop("Seasonal advection coefficients need a seasonal period; set one ",
         "with period(dat) <- ...", call. = FALSE)
  }
  breaks <- .season_breaks(period, nsea)
  function(t) max(1L, t2index(t, breaks, period = period, seasonal = TRUE))
}


## Advection velocity at positions xy (n x 2) and time t, as an n x 2 matrix:
## adv_const[, s] + sum_f gamma[, f, s] * (u_f, v_f)(xy, t), with s the season.
## Row 1 of gamma / adv_const acts on x, row 2 on y. Called inside nll(), so it
## must keep working with advector gamma and adv_const.
.make_adv <- function(adv, time_adv, gamma, adv_const, period = NULL) {

  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")

  nfield <- length(adv)
  liv <- if (nfield > 0L) .get_liv(.adv_flatten(adv), deriv = FALSE)$liv else list()

  nsea <- if (!is.null(adv_const)) ncol(adv_const)
          else if (!is.null(gamma)) dim(gamma)[3L] else 1L
  season <- .adv_season_fun(nsea, period)

  val <- function(xy, t) {
    if (is.null(dim(xy))) xy <- matrix(xy, ncol = 2L)
    s <- season(t)
    n <- nrow(xy)
    ax <- ay <- rep(0, n)
    if (!is.null(adv_const)) {
      ax <- ax + adv_const[1L, s]
      ay <- ay + adv_const[2L, s]
    }
    for (f in seq_len(nfield)) {
      iu <- 2L * f - 1L
      it <- t2index(t, time_adv[[iu]])
      if (it > 0L) {
        ax <- ax + gamma[1L, f, s] * liv[[iu]][[it]](xy[, 1], xy[, 2])
        ay <- ay + gamma[2L, f, s] * liv[[iu + 1L]][[it]](xy[, 1], xy[, 2])
      }
    }
    cbind(ax, ay)
  }

  ## the season and field time slices val() reads at t; the CTMC generator cache
  ## keys on it (.ctmc_slice_key()), so it must change whenever val() does
  slice <- function(t) {
    base::c(season(t), vapply(seq_len(nfield), function(f) {
      as.integer(t2index(t, time_adv[[2L * f - 1L]]))
    }, integer(1L)))
  }

  list(val = val, slice = slice, nfield = nfield, nsea = nsea)
}


## Subset or combine covariate lists keeping the class and the list-level sref
## and tref (plain `[` and c() drop them).
.cov_list_keep <- function(x, template) {
  if (length(x) == 0L) return(NULL)
  x <- .add_class(x, "admove_cov_list")
  for (a in c("sref", "tref")) {
    v <- attr(template, a)
    if (!is.null(v)) attr(x, a) <- v
  }
  x
}


## Fill advection settings missing from a configuration (e.g. one written by hand).
.adv_conf <- function(conf) {
  if (is.null(conf$adv_gamma)) conf$adv_gamma <- "shared"
  if (is.null(conf$adv_const)) conf$adv_const <- FALSE
  if (is.null(conf$n_seasons_adv)) conf$n_seasons_adv <- 1L
  if (is.null(conf$use_advection)) conf$use_advection <- FALSE
  conf
}


## Number of advection seasons of a parameter list.
.adv_nsea <- function(par) {
  if (!is.null(par$adv_const)) return(ncol(par$adv_const))
  if (!is.null(par$gamma)) return(dim(par$gamma)[3L])
  1L
}
