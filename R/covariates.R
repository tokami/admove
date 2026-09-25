## Main functions ---------------------------------------------------------------

##' Prepare covariate fields for admove
##'
##' @description
##' `prep_cov()` converts covariate data into the standard array-based format
##' used by \emph{admove}. Covariates can be supplied as a list of matrices, as a
##' 2D or 3D array, or as raster-like objects. When the input has multiple
##' layers, the `layers` argument controls whether they are treated as
##' separate covariate fields (default, returns a named `admove_cov_list`) or
##' as consecutive time steps of a single covariate (returns a single
##' `admove_cov`).
##'
##' The function can also attach spatial and temporal reference information,
##' convert date-like time labels to numeric model time, and optionally plot the
##' prepared covariate field.
##'
##' @param x Covariate data. Supported inputs include:
##'   \itemize{
##'     \item a `data.frame` with columns `x` and `y` (coordinates) and one
##'       additional column per covariate; each row is one grid cell. The
##'       covariate column names become the names of the returned
##'       `admove_cov_list`. This is the natural output of GIS exports and
##'       spatial model pipelines.
##'     \item a list of matrices, typically one matrix per time step,
##'     \item a 2D array or matrix, interpreted as a single time slice,
##'     \item a 3D array with dimensions x, y, and time,
##'     \item `RasterLayer`, `RasterBrick`, or `RasterStack` objects, and
##'     \item `SpatRaster` objects.
##'   }
##' @param x_centers Optional numeric vector giving x coordinates of cell
##'   centres.
##' @param y_centers Optional numeric vector giving y coordinates of cell
##'   centres.
##' @param times Optional vector giving the time values associated with the third
##'   dimension.
##' @param date_decimal Logical; if `TRUE`, interpret time labels as decimal
##'   years and convert them to model time. Default: `FALSE`.
##' @param date_format Optional format string (see [strptime()]) to parse
##'   character time labels, e.g. `"%Y-%m-%d %H:%M"`. The time of day is kept.
##'   Default: `NULL`.
##' @param date_origin Optional origin when time labels are stored numerically
##'   as (possibly fractional) days since `date_origin`. Default: `NULL`.
##' @param tz Time zone used when converting dates. Default: `"UTC"`.
##' @param sref Optional spatial reference information to attach to the returned
##'   object.
##' @param tref Optional time reference information to attach to the returned
##'   object.
##' @param plot Logical; if `TRUE`, plot the prepared covariate field. Default:
##'   `FALSE`.
##' @param plot_land Logical; passed to the plotting method when `plot = TRUE`.
##'   Default: `FALSE`.
##' @param strict Logical; if `TRUE`, require stricter matching and validation of
##'   dimension names. Default: `FALSE`.
##' @param verbose Logical; if `TRUE`, print informative messages. Default:
##'   `TRUE`.
##'
##' @param layers Character string controlling how layers in multi-layer
##'   inputs are interpreted. Applies to `SpatRaster`, `RasterBrick`,
##'   `RasterStack`, lists of matrices, and 3-D arrays. Use `"covariates"` to
##'   treat each layer as a separate covariate field; the function then returns
##'   a named list of class `admove_cov_list`. Use `"time"` to treat layers as
##'   consecutive time steps of a single covariate; the function then returns a
##'   single `admove_cov`. Default (`NULL`) is `"covariates"` for
##'   `SpatRaster`/`Raster*` objects and `"time"` for lists and arrays, which
##'   preserves backward-compatible behaviour for those types. Ignored for
##'   single-layer inputs.
##'
##' @return
##' When `layers = "time"` or when the input has only one layer, an
##' object of class `admove_cov`, stored as a 3D array with dimensions
##' corresponding to x, y, and time. When `layers = "covariates"` and
##' the input has multiple layers, a named list of class `admove_cov_list` with
##' one `admove_cov` element per layer, ready to pass directly to [setup_data()].
##'
##' @details
##' If `x` is two-dimensional, it is converted to a 3D array with a single time
##' slice. Dimension names are validated and, where possible, inferred from the
##' input or from the optional `x_centers`, `y_centers`, and `times` arguments.
##'
##' If date-like time labels are supplied, they can be converted to numeric model
##' time using `date_format`, `date_origin`, or `date_decimal`. Spatial and
##' temporal reference information are preserved from the input where available
##' and supplemented by `sref` and `tref` if provided. A user-supplied `sref`
##' takes precedence over any CRS derived from a Raster or SpatRaster object.
##'
##' @examples
##' cov <- prep_cov(skjepo$cov)
##'
##' @export
prep_cov <- function(x,
                     x_centers = NULL,
                     y_centers = NULL,
                     times = NULL,
                     date_decimal = FALSE,
                     date_format = NULL,
                     date_origin = NULL,
                     tz = "UTC",
                     sref = NULL,
                     tref = NULL,
                     layers = NULL,
                     plot = FALSE,
                     plot_land = FALSE,
                     strict = FALSE,
                     verbose = TRUE) {

  ## data.frame → named list of matrices (one per covariate column) ----------
  if (inherits(x, "data.frame")) {
    nms <- names(x)
    x_col <- if ("x" %in% nms) "x" else if ("X" %in% nms) "X" else {
      num_nms <- nms[vapply(x, is.numeric, logical(1L))]
      if (length(num_nms) >= 1L) num_nms[1L] else
        stop("Cannot identify the x-coordinate column in the data.frame. ",
             "Name a column 'x', or ensure the first numeric column is x.")
    }
    y_col <- if ("y" %in% nms) "y" else if ("Y" %in% nms) "Y" else {
      num_nms <- nms[vapply(x, is.numeric, logical(1L))]
      if (length(num_nms) >= 2L) num_nms[2L] else
        stop("Cannot identify the y-coordinate column in the data.frame. ",
             "Name a column 'y', or ensure the second numeric column is y.")
    }
    cov_cols <- setdiff(nms, c(x_col, y_col))
    if (length(cov_cols) == 0L)
      stop("data.frame has no covariate columns beyond '", x_col,
           "' and '", y_col, "'.")
    xu <- sort(unique(x[[x_col]]))
    yu <- sort(unique(x[[y_col]]))
    ix <- match(x[[x_col]], xu)
    iy <- match(x[[y_col]], yu)
    x_lab <- sprintf("%.2f", xu)
    y_lab <- sprintf("%.2f", yu)
    cov_mats <- lapply(cov_cols, function(nm) {
      m <- matrix(NA_real_, length(xu), length(yu),
                  dimnames = list(x = x_lab, y = y_lab))
      m[cbind(ix, iy)] <- x[[nm]]
      m
    })
    names(cov_mats) <- cov_cols
    x <- cov_mats
    if (is.null(layers)) layers <- "covariates"
  }

  if (is.null(layers)) {
    layers <- if (inherits(x, "SpatRaster") ||
                         inherits(x, "RasterBrick") ||
                         inherits(x, "RasterStack")) "covariates" else "time"
  } else {
    layers <- match.arg(layers, c("covariates", "time"))
  }
  sref_from_raster <- NULL

  dimnams0 <- dimnames(x)
  sref0 <- try(sref(x), silent = TRUE)
  tref0 <- try(tref(x), silent = TRUE)

  ## When layers = "covariates" and the input has multiple layers/elements,
  ## split into individual prep_cov calls and return a named admove_cov_list.
  if (layers == "covariates") {

    nl_check <- NULL
    nms_check <- NULL
    sref_for_layers <- sref  # user-supplied sref; fall back to raster CRS if NULL

    if (inherits(x, "SpatRaster")) {
      nl_check <- terra::nlyr(x)
      nms_check <- names(x)
      if (is.null(nms_check) || any(!nzchar(nms_check))) nms_check <- paste0("layer", seq_len(nl_check))
      if (is.null(sref_for_layers) && requireNamespace("sf", quietly = TRUE)) {
        crs_sf <- try(sf::st_crs(terra::crs(x)), silent = TRUE)
        if (!inherits(crs_sf, "try-error"))
          sref_for_layers <- list(crs = crs_sf, units = crs_sf$units_gdal, crs_scale = 1)
      }
    } else if (inherits(x, "RasterBrick") || inherits(x, "RasterStack")) {
      nl_check <- raster::nlayers(x)
      nms_check <- names(x)
      if (is.null(nms_check) || any(!nzchar(nms_check))) nms_check <- paste0("layer", seq_len(nl_check))
      if (is.null(sref_for_layers) && requireNamespace("sf", quietly = TRUE)) {
        crs_sf <- try(sf::st_crs(as.character(raster::crs(x))), silent = TRUE)
        if (!inherits(crs_sf, "try-error"))
          sref_for_layers <- list(crs = crs_sf, units = crs_sf$units_gdal, crs_scale = 1)
      }
    } else if (is.list(x)) {
      nl_check <- length(x)
      nms_check <- names(x)
      if (is.null(nms_check)) nms_check <- as.character(seq_len(nl_check))
    } else if (is.array(x) && length(dim(x)) == 3L) {
      nl_check <- dim(x)[3L]
      nms_check <- dimnames(x)[[3]]
      if (is.null(nms_check)) nms_check <- as.character(seq_len(nl_check))
    }

    if (!is.null(nl_check) && nl_check > 1L) {
      collected_msgs <- character(0)
      cov_list <- lapply(seq_len(nl_check), function(i) {
        layer_i <- if (inherits(x, "SpatRaster")) {
          terra::subset(x, i)
        } else if (inherits(x, "RasterBrick") || inherits(x, "RasterStack")) {
          raster::subset(x, i)
        } else if (is.list(x)) {
          x[[i]]
        } else {
          x[, , i]
        }
        withCallingHandlers(
          prep_cov(layer_i,
                   x_centers = x_centers, y_centers = y_centers,
                   times = times,
                   date_decimal = date_decimal,
                   date_format = date_format, date_origin = date_origin,
                   tz = tz,
                   sref = sref_for_layers, tref = tref,
                   layers = "time",
                   plot = FALSE, strict = strict, verbose = verbose),
          message = function(m) {
            collected_msgs <<- c(collected_msgs, conditionMessage(m))
            invokeRestart("muffleMessage")
          }
        )
      })
      if (verbose) {
        for (msg in unique(trimws(collected_msgs))) {
          if (nzchar(msg)) message(msg)
        }
      }
      names(cov_list) <- nms_check
      cov_list <- .add_class(cov_list, "admove_cov_list")
      if (plot) lapply(cov_list, function(co) plot_cov(co, plot_land = plot_land))
      return(cov_list)
    }

  }

  if (inherits(x, "RasterLayer") || inherits(x, "RasterBrick") ||
        inherits(x, "RasterStack")) {

    if (!requireNamespace("raster", quietly = TRUE)) {
      stop("Package 'raster' is required to convert Raster* objects. Please install it or convert to an array first. (See vignette)")
    }

    nr <- raster::nrow(x)
    nc <- raster::ncol(x)
    nl <- raster::nlayers(x)

    x_lab <- sprintf("%.2f", raster::xFromCol(x, seq_len(nc)))
    y_lab <- sprintf("%.2f", rev(raster::yFromRow(x, seq_len(nr))))

    if (nl == 1L) {
      m <- raster::as.matrix(x)
      m <- m[nr:1, , drop = FALSE]
      arr <- t(m)
      dimnames(arr) <- list(x = x_lab, y = y_lab)
      x <- arr
    } else {
      a <- raster::as.array(x)
      a <- a[nr:1, , , drop = FALSE]
      arr <- aperm(a, c(2, 1, 3))

      layer_lab <- names(x)
      if (is.null(layer_lab) || any(!nzchar(layer_lab))) {
        layer_lab <- seq_len(nl)
      }

      dimnames(arr) <- list(x = x_lab, y = y_lab, layer = layer_lab)
    }

    p4 <- as.character(raster::crs(x))

    if (!requireNamespace("sf", quietly = TRUE)) {
      stop("Package 'sf' is required to extract information from sf objects. Please install it or provide information to create_grid and set x to NULL.")
    }

    sref_from_raster <- list()
    sref_from_raster$crs <- sf::st_crs(p4)
    sref_from_raster$units <- sref_from_raster$crs$units_gdal
    sref_from_raster$crs_scale <- 1


  } else if (inherits(x, "SpatRaster")) {

    if (!requireNamespace("terra", quietly = TRUE)) {
      stop("Package 'terra' is required to convert SpatRaster* objects. Please install it or convert to an array first. (see vignette)")
    }

    nr <- terra::nrow(x)
    nc <- terra::ncol(x)
    nl <- terra::nlyr(x)

    ext_r <- terra::ext(x)
    res_r <- terra::res(x)
    x_cent <- ext_r[1] + res_r[1]/2 + (0:(nc-1)) * res_r[1]
    y_cent <- ext_r[3] + res_r[2]/2 + (0:(nr-1)) * res_r[2]

    x_lab <- sprintf("%.2f", x_cent)
    y_lab <- sprintf("%.2f", y_cent)

    if (nl == 1L) {

      ## x-y array
      m <- terra::as.matrix(x, wide = TRUE)
      m <- m[nr:1, , drop = FALSE]
      arr <- t(m)
      dimnames(arr) <- list(x = x_lab, y = y_lab)

    } else {

      a <- terra::as.array(x)
      a <- a[nr:1, , , drop = FALSE]
      arr <- aperm(a, c(2, 1, 3))

      layer_lab <- names(x)
      if (is.null(layer_lab) || any(!nzchar(layer_lab))) {
        layer_lab <- paste0("layer", seq_len(nl))
      }

      if (layers == "time") {
        dimnames(arr) <- list(x = x_lab, y = y_lab, time = layer_lab)
      }

      if (layers == "covariates") {
        dimnames(arr) <- list(x = x_lab, y = y_lab, covariate = layer_lab)
      }
    }

    if (!requireNamespace("sf", quietly = TRUE)) {
      stop("Package 'sf' is required to extract CRS units. Please install it.")
    }

    crs_str <- terra::crs(x)
    if (nzchar(crs_str)) {
      crs_sf <- try(sf::st_crs(crs_str), silent = TRUE)
      if (!inherits(crs_sf, "try-error")) {
        sref_from_raster <- list(crs = crs_sf, units = crs_sf$units_gdal, crs_scale = 1)
      }
    }

  } else if (inherits(x, "list")) {

    arr <- list_2_3Darray(x)

  } else {

    arr <- x

  }

  if (length(dim(arr)) == 2) {

    res <- array(arr, c(dim(arr), 1))
    if (!is.null(dimnames(arr))) dimnames(res) <- c(dimnames(arr), list(NULL))
    dimnams0 <- c(dimnams0, list(NULL))

  } else {

    res <- arr

  }

  dimnames(res) <- .validate_dimnames(list(x_centers,
                                           y_centers,
                                           times),
                                      dim(res),
                                      dimnames(res),
                                      dimnams0,
                                      strict = strict,
                                      verbose = verbose)

  dati <- dimnames(res)[[3]]

  ## convert date
  if (!is.null(date_origin) ||
        !is.null(date_format) ||
         isTRUE(date_decimal)) {

    dati <- .parse_dates(dati, date_format = date_format,
                         date_origin = date_origin,
                         date_decimal = date_decimal, tz = tz)

    ## convert date to numeric
    d2t <- date_2_time(dati, tref)
    dati <- d2t
    attributes(dati) <- NULL
    tref <- create_tref(attr(d2t, "tref")[["origin"]],
                        attr(d2t, "tref")[["units"]])
    if (isTRUE(attr(d2t, "tref")[["inferred"]])) {
      if (verbose) message("tref (time origin and units) was inferred from dates. Please check and adjust if needed.")
    }
  }

  dati_num <- as.numeric(as.character(dati))
  if (all(is.na(dati_num))) stop("The time information cannot be interpreted as numeric. Please provide the time information as a numeric value or use the 'date_format' and/or 'date_origin' argument to convert the date into a numeric value (see ?as.Date). If your input has multiple columns representing time steps of a single covariate (e.g. one column per date), consider setting layers = \"time\".")

  dimnames(res)[[3]] <- dati


  res <- .add_class(res, "admove_cov")
  if (!is.null(sref0) && !inherits(sref0, "try-error")) {
    res <- add_sref(res, sref0)
  }
  if (!is.null(tref0) && !inherits(tref0, "try-error")) {
    res <- add_tref(res, tref0)
  }
  if (!is.null(sref_from_raster)) res <- add_sref(res, sref_from_raster)
  res <- add_sref(res, sref)
  res <- add_tref(res, tref)

  if (plot) plot(res, plot_land = plot_land)

  return(res)
}



##' Plot covariate fields
##'
##' @description
##' `plot_cov()` plots one or more covariate fields stored in an `admove_cov`
##' object or in a higher-level \emph{admove} object that contains covariate
##' data.
##'
##' Each time slice is plotted as a separate panel, optionally with a colour
##' bar, contour lines, and land added to the plot.
##'
##' @param x An object of class `admove_cov`, or an object containing covariate
##'   data such as `admove_data`, `admove_sim`, or `admove`.
##' @param i Index (scalar or vector) used when `x` is an `admove_cov_list`.
##'   A scalar selects one covariate and plots all its time steps. A vector
##'   selects multiple covariates and produces one panel per element (showing
##'   the first time step, or the first element of `select`). Default: `1`.
##' @param select Optional vector of time-step indices to plot, i.e. indices of
##'   the covariate's time slices (not model times). Default: `NULL`, in which
##'   case all time steps are plotted. To plot the slices covering given model
##'   times, e.g. the prediction times of a fit, match them to the slice start
##'   times with [findInterval()] (see Examples).
##' @param main Main title of the plot. Default: `"Covariate fields"`.
##' @param labels Logical; currently reserved for plotting cell labels. Default:
##'   `TRUE`.
##' @param plot_land Logical; if `TRUE`, add land to the plot. Default:
##'   `FALSE`.
##' @param auto_layout Logical; if `TRUE`, plotting parameters are set
##'   automatically and restored afterwards. Default: `TRUE`.
##' @param xlab Label for the x-axis. If `NULL` (default), `"x"` with the
##'   spatial units in brackets, e.g. `"x [km]"` (`"lon [°]"` for degrees).
##' @param ylab Label for the y-axis. If `NULL` (default), `"y"` with the
##'   spatial units in brackets, e.g. `"y [km]"` (`"lat [°]"` for degrees).
##' @param bg Optional background colour for the plot. Default: `NULL`.
##' @param plot_contour Logical; if `TRUE`, add contour lines. Default:
##'   `TRUE`.
##' @param xlim xlim
##' @param ylim ylim
##' @param col Colour palette for the covariate values. Default:
##'   `hcl.colors(100, "viridis")`.
##' @param zlim Optional numeric range of values mapped onto `col`. Default:
##'   `NULL`, in which case the range over all plotted time steps of a covariate
##'   is used, so that time steps of the same covariate share one colour scale.
##'   With several covariates (`i` a vector), each covariate gets its own scale.
##' @param legend Logical; if `TRUE`, a colour bar is drawn to the right of each
##'   panel. Default: `NULL`, which means `TRUE` when `auto_layout = TRUE` and
##'   `FALSE` otherwise (when the caller controls the margins, there may be no
##'   room for the bar).
##' @param titles Optional character vector of panel titles, one per panel.
##'   Default: `NULL`, in which case covariate names (or `"Covariate i"`) are
##'   used for several covariates and the start time of each slice for several
##'   time steps, as a date at the resolution of the time units (e.g.
##'   `"Jan 2007"` for monthly units; `"t = 48"` if the time reference has no
##'   origin); a single panel gets no title.
##' @param land_col Fill colour for land. Default: `grey(0.85)` (opaque, so
##'   covariate colours do not show through).
##' @param land_border Border colour for land. Default: `grey(0.5)`.
##' @param ... Additional graphical arguments passed to [plot()].
##'
##' @return
##' Invisibly returns `NULL`.
##'
##' @details
##' The function uses the x and y dimension names, if present, as plotting
##' coordinates. Otherwise, row and column indices are used. Each selected time
##' slice is displayed with [graphics::image()], and optional contours are added
##' with [graphics::contour()].
##'
##' @examples
##' plot_cov(skjepo$cov)
##'
##' \dontrun{
##' ## covariate slices covering the prediction times of a fit
##' tt <- as.numeric(dimnames(fit$dat$cov[[1]])[[3]])
##' plot_cov(fit, 1, select = findInterval(fit$dat$pred$time, tt))
##' }
##'
##' @name plot_cov
##' @export
plot_cov <- function(x,
                     i = 1,
                     select = NULL,
                     main = "Covariate fields",
                     labels = TRUE,
                     plot_land = FALSE,
                     auto_layout = TRUE,
                     xlab = NULL,
                     ylab = NULL,
                     bg = NULL,
                     plot_contour = TRUE,
                     xlim = NULL,
                     ylim = NULL,
                     col = hcl.colors(100, "viridis"),
                     zlim = NULL,
                     legend = NULL,
                     titles = NULL,
                     land_col = grey(0.85),
                     land_border = grey(0.5),
                     ...) {

  map_labs <- .map_labs(x)
  if (is.null(xlab)) xlab <- map_labs[1]
  if (is.null(ylab)) ylab <- map_labs[2]

  xlim0 <- xlim
  ylim0 <- ylim
  if (is.null(legend)) legend <- isTRUE(auto_layout)

  if (inherits(x, "admove_sim")) {
    cov <- x$cov
  } else if(inherits(x, "admove_data")) {
    cov <- x$cov
  } else if(inherits(x, "admove")) {
    cov <- x$dat$cov
  } else{
    cov <- x
  }

  ## right margin wide enough for the colour bar and its labels
  mar_right <- if (legend) 4.5 else 1.5

  ## with several panels, draw the axes only on the outer panels (x on the
  ## lowest panel of each column, y on the first column); the axis labels go
  ## once into the outer margin
  dots <- list(...)
  xaxt <- if (is.null(dots$xaxt)) "s" else dots$xaxt
  yaxt <- if (is.null(dots$yaxt)) "s" else dots$yaxt
  dots$xaxt <- dots$yaxt <- NULL

  if (inherits(cov, "admove_cov_list") && length(i) > 1) {
    sel <- cov[i]
    n <- length(sel)
    nms <- names(sel)
    t_sel <- if (is.null(select)) 1L else select[1L]
    if (!is.null(titles) && length(titles) != n) {
      stop("'titles' must have one entry per covariate (", n, ").")
    }
    if (auto_layout) {
      opar <- par(no.readonly = TRUE)
      on.exit(par(opar))
      par(mfrow = n2mfrow(n, asp = 2),
          mar = c(0.3, 0.3, 1.4, mar_right),
          oma = c(3, 3.5, ifelse(main == "", 0, 1.5), 0),
          mgp = c(2, 0.5, 0),
          tcl = -0.3)
    }
    for (j in seq_along(sel)) {
      panel_lbl <- if (!is.null(titles)) {
        titles[j]
      } else if (!is.null(nms) && nzchar(nms[j])) {
        nms[j]
      } else {
        paste0("Covariate ", i[j])
      }
      do.call(plot_cov,
              c(list(sel[[j]], select = t_sel, main = "",
                     plot_land = plot_land, auto_layout = FALSE,
                     xlab = xlab, ylab = ylab, bg = bg,
                     plot_contour = plot_contour,
                     xlim = xlim, ylim = ylim,
                     col = col, zlim = zlim, legend = legend,
                     titles = panel_lbl,
                     land_col = land_col, land_border = land_border,
                     xaxt = if (auto_layout) "n" else xaxt,
                     yaxt = if (auto_layout) "n" else yaxt),
                dots))
      if (auto_layout) .outer_panel_axes(j, n, xaxt, yaxt)
    }
    if (auto_layout) {
      mtext(main, 3, 0, outer = TRUE)
      mtext(xlab, 1, 2, outer = TRUE)
      mtext(ylab, 2, 2.2, outer = TRUE)
    }
    return(invisible(NULL))
  }

  if (inherits(cov, "list")) {
    cov <- cov[[i]]
  }

  .check_class(cov, "admove_cov")
  sref <- sref(cov)

  if (!is.null(select)) {
    if (max(select) > dim(cov)[3]) stop("Trying to select times (select) that are outside of the covariate field!")
    cov <- cov[,,select, drop = FALSE]
  }

  nt <- dim(cov)[3]

  if (is.null(titles)) {
    titles <- if (nt > 1) {
      .time_labels(as.numeric(dimnames(cov)[[3]]), tref(cov))
    } else ""
  } else if (length(titles) != nt) {
    stop("'titles' must have one entry per plotted time step (", nt, ").")
  }

  ## one colour scale for all time steps of this covariate
  if (is.null(zlim)) {
    zlim <- suppressWarnings(range(unclass(cov), na.rm = TRUE, finite = TRUE))
  }
  zlim_ok <- all(is.finite(zlim))
  if (zlim_ok && zlim[1] == zlim[2]) {
    zlim <- zlim + c(-0.5, 0.5) * max(abs(zlim[1]), 1)
  }

  if (is.null(xlim0)) {
    if(any(names(attributes(cov)) == "dimnames")){
      xlims <- range(as.numeric(attributes(cov)$dimnames[[1]]))
    }else{
      xlims <- c(0,1)
    }
  } else {
    xlims <- xlim0
  }

  if (is.null(ylim0)) {
    if(any(names(attributes(cov)) == "dimnames")){
      ylims <- range(as.numeric(attributes(cov)$dimnames[[2]]))
    }else{
      ylims <- c(0,1)
    }
  } else {
    ylims <- ylim0
  }

  if(auto_layout){
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar))
    par(mfrow = n2mfrow(nt, asp = 2),
        mar = c(if (nt > 1) 0.3 else 1.5, if (nt > 1) 0.3 else 1.5,
                ifelse(nt > 1, 1.4, 1.5), mar_right),
        oma = c(3, 3.5, ifelse(main == "", 0, 1.5), 0),
        mgp = c(2, 0.5, 0),
        tcl = -0.3)
  }
  shared <- auto_layout && nt > 1
  for(i in 1:nt){
    x <- as.numeric(rownames(cov[,,i]))
    if(length(x) == 0) x <- 1:nrow(cov[,,i])
    y <- as.numeric(colnames(cov[,,i]))
    if(length(y) == 0) y <- 1:ncol(cov[,,i])
    if(!is.null(bg)){
      par(bg = bg)
    }
    do.call(plot, c(list(1, 1, type = "n",
                         xlim = xlims, ylim = ylims,
                         xlab = "",
                         ylab = "",
                         asp = 1,
                         xaxt = if (shared) "n" else xaxt,
                         yaxt = if (shared) "n" else yaxt),
                    dots))
    if (shared) .outer_panel_axes(i, nt, xaxt, yaxt)
    z <- cov[,,i, drop = TRUE]
    if (zlim_ok) {
      ## clamp to zlim so values outside a user-supplied range are not left blank
      z <- pmin(pmax(z, zlim[1]), zlim[2])
      image(x, y, z, col = col, zlim = zlim, add = TRUE)
    }
    if (plot_land) {
      plot_land(sref, col = land_col, border = land_border)
    }
    if(plot_contour && zlim_ok) {
      contour(x, y, cov[,,i], add = TRUE, col = grey(0.2, 0.6),
              labcex = 0.6)
    }
    if (nzchar(titles[i])) {
      title(main = titles[i], line = if (shared) 0.3 else 0.4,
            font.main = 1, cex.main = 1)
    }
    box(lwd = 1.5)
    if (legend && zlim_ok) .color_bar(col, zlim)
  }
  if(auto_layout){
    mtext(main, 3, 0, outer = TRUE)
    mtext(xlab, 1, if (shared) 2 else 1, outer = TRUE)
    mtext(ylab, 2, if (shared) 2.2 else 1.5, outer = TRUE)
  }


  return(invisible(NULL))
}


## Axes of panel k of n in a multi-panel layout: x on the lowest panel of each
## column, y on the first column.
.outer_panel_axes <- function(k, n, xaxt = "s", yaxt = "s") {
  ncol_lay <- par("mfrow")[2L]
  if (xaxt != "n" && k + ncol_lay > n) axis(1)
  if (yaxt != "n" && (k - 1L) %% ncol_lay == 0L) axis(2)
}


## Vertical colour bar just right of the current plot region, in the figure
## margin (needs a right margin of about 4 lines).
.color_bar <- function(col, zlim, width = 0.03, gap = 0.015, cex = 0.7) {

  plt <- par("plt")
  x0 <- grconvertX(plt[2] + gap, from = "nfc", to = "user")
  x1 <- grconvertX(plt[2] + gap + width, from = "nfc", to = "user")
  y0 <- grconvertY(plt[3], from = "nfc", to = "user")
  y1 <- grconvertY(plt[4], from = "nfc", to = "user")

  n <- length(col)
  yb <- seq(y0, y1, length.out = n + 1)
  op <- par(xpd = NA)
  on.exit(par(op))
  rect(x0, yb[-(n + 1)], x1, yb[-1], col = col, border = NA)
  rect(x0, y0, x1, y1, border = grey(0.3))

  at <- pretty(zlim)
  at <- at[at >= zlim[1] & at <= zlim[2]]
  yat <- y0 + (at - zlim[1]) / diff(zlim) * (y1 - y0)
  xt <- grconvertX(plt[2] + gap + width + 0.01, from = "nfc", to = "user")
  segments(x1, yat, x1 + 0.3 * (x1 - x0), yat, col = grey(0.3))
  text(xt, yat, labels = format(at), adj = c(0, 0.5), cex = cex)

  invisible(NULL)
}




##' Summarise covariate fields
##'
##' @description
##' `summarise_cov()` prints a compact summary of one or more covariate fields.
##' It works on objects of class `admove_cov`, `admove_cov_list`, or on
##' higher-level \emph{admove} objects that contain covariate data.
##'
##' The summary includes dimensions, cell size, spatial and temporal ranges,
##' covariate value range, number of missing cells, and the associated spatial
##' and temporal units where available.
##'
##' @param object An object of class `admove_cov`, `admove_cov_list`, or an object
##'   containing covariate data such as `admove_data`, `admove_sim`, or
##'   `admove`.
##' @param ... Additional arguments
##'
##' @return
##' Invisibly returns the corresponding covariate object, coerced internally to
##' a covariate list if needed.
##'
##' @examples
##' summarise_cov(skjepo$cov)
##'
##' @name summarise_cov
##' @export
summarise_cov <- function(object, ...) {
  x <- object

  if(inherits(x, "admove_sim")) {
    cov <- x$cov
  } else if(inherits(x, "admove_data")) {
    cov <- x$cov
  } else if(inherits(x, "admove")) {
    cov <- x$dat$cov
  } else if(inherits(x, "admove_cov")) {
    cov <- x
  }  else if(inherits(x, "admove_cov_list")) {
    cov <- x
  } else stop("Please provide an object of class 'admove_cov' or an object containing an such an object (e.g. admove_data, admove_sim, admove).")


  cov <- .make_cov_list(cov)
  ncov <- length(cov)

  for (i in 1:ncov) {
    covi <- cov[[i]]

    rans <- sprintf("%.2f", range(covi, na.rm = TRUE))
    dims <- dim(covi)

    xcen <- as.numeric(dimnames(covi)[[1]])
    ycen <- as.numeric(dimnames(covi)[[2]])
    tstart <- as.numeric(dimnames(covi)[[3]])
    cellsize <- c(median(diff(xcen)),
                  median(diff(ycen)))
    xrange <- range(xcen - cellsize[1]/2, xcen + cellsize[1]/2)
    yrange <- range(ycen - cellsize[2]/2,  ycen + cellsize[2]/2)

    tmp <- apply(covi, c(1,2), function(x) any(is.na(x)))
    n_na <- sum(tmp)

    units_sp <- try(units_space(covi), silent = TRUE)
    if (is.null(units_sp) || is.na(units_sp) || units_sp == "" || inherits(units_sp, "try-error")) units_sp <- "not specified"

    units_t <- try(units_time(covi), silent = TRUE)
    if (is.null(units_t) || is.na(units_t) || units_t == "" || inherits(units_t, "try-error")) units_t <- "not specified"

    rans_t <- sprintf("%.2f", range(as.numeric(dimnames(covi)[[3]])))

    rans_t_abs <- NULL
    if (!is.na(tref(covi)$origin)) {
      u <- .normalise_time_unit(tref(covi)$units)
      t_rel <- range(as.numeric(dimnames(covi)[[3]]))
      if (!is.na(u)) {
        rans_t_abs <- .add_rel_time_posix(tref(covi)$origin, t_rel, units = u,
                                          days_per_year = 365)
        rans_t_abs <- as.character(format(rans_t_abs))
      }
    }

    labw <- 10

    if (ncov == 1) {
      cat("<admove_cov>\n")
    } else {
      cat(paste0("<admove_cov> --- [",i,"]\n"))
    }
    cat(sprintf(paste0("  %-", labw, "s %s\n"), "cells:",
                prod(dims[1:2])))
    cat(sprintf(paste0("  %-", labw, "s %s\n"), "dims:",
                paste(dims[1], "x", dims[2], "x", dims[3])))
    cat(sprintf(paste0("  %-", labw, "s %s\n"), "cellsize:",
                paste(cellsize[1], "x", cellsize[1])))
    cat(sprintf(paste0("  %-", labw, "s %s\n"), "xrange:",
                paste0("[",sprintf("%.2f", xrange[1]), ", ",
                       sprintf("%.2f", xrange[2]),"]")))
    cat(sprintf(paste0("  %-", labw, "s %s\n"), "yrange:",
                paste0("[",sprintf("%.2f", yrange[1]), ", ",
                       sprintf("%.2f", yrange[2]),"]")))
    cat(sprintf(paste0("  %-", labw, "s %s\n"), "trange:",
                paste0("[",rans_t[1], ", ", rans_t[2],"]")))
    if (!is.null(rans_t_abs)) {
      cat(sprintf(paste0("  %-", labw, "s %s\n"), "",
                  paste0("[",rans_t_abs[1])))
      cat(sprintf(paste0("  %-", labw, "s %s\n"), "",
                  paste0("\t",rans_t_abs[2],"]")))
    }
    cat(sprintf(paste0("  %-", labw, "s %s\n"), "cov range:",
                paste0("[",rans[1], ", ", rans[2],"]")))
    cat(sprintf(paste0("  %-", labw, "s %s\n"), "NAs:",
                n_na))
    cat(sprintf(paste0("  %-", labw, "s %s\n"), "units:",
                paste0(units_sp, " x ", units_t)))

  }


  invisible(cov)
}





##' Check and standardise covariate data
##'
##' @description
##' `check_cov()` validates covariate data and ensures that it can be used by
##' \emph{admove}. If the input is not already of class `admove_cov` or
##' `admove_cov_list`, the function first tries to convert it using
##' [prep_cov()].
##'
##' Existing covariate objects are checked for missing values, missing dimension
##' names, and incompatible dimension-name lengths. Invalid covariate fields are
##' removed.
##'
##' @param x Covariate data to check. Can be `NULL`, an object of class
##'   `admove_cov` or `admove_cov_list`, or another object that can be converted
##'   with [prep_cov()].
##' @param verbose Logical; if `TRUE`, print informative messages about removed
##'   covariate fields. Default: `TRUE`.
##'
##' @return
##' `NULL` if no valid covariate field remains; otherwise an object of class
##' `admove_cov` or `admove_cov_list`.
##'
##' @details
##' For existing covariate objects, the function removes covariate fields that:
##' \itemize{
##'   \item contain only missing values,
##'   \item have missing or incomplete dimension names, or
##'   \item have dimension names whose lengths do not match the array dimensions.
##' }
##'
##' @export
check_cov <- function(x, verbose = TRUE) {

  if (is.null(x)) {

    return(NULL)

  } else if(!inherits(x, "admove_cov") &&
              !inherits(x, "admove_cov_list")) {

    res <- try(prep_cov(x), silent = TRUE)

    if (inherits(res, "try-error")) {

      stop("Provided object is not of class admove_cov and couldn't convert covariates to admove_cov using prep_cov. Please check your covariate object.")

    } else {

      return(res)

    }

  } else {

    x <- .make_cov_list(x)

    ## all entries NA
    idx <- which(sapply(x, function(x) all(is.na(x))))
    if (length(idx) > 0) {
      if (verbose) writeLines(paste0("All entries in covariate field(s) ",
                                 paste(idx, collapse = ","),
                                 " are NA! This covariate field is removed."))
      x <- x[-idx]
    }
    if (length(x) == 0) return(NULL)

    ## any dimnames missing
    idx <- which(sapply(x, function(x) is.null(dimnames(x)) || length(dimnames(x)) != 3))
    if (length(idx) > 0) {
      if (verbose) writeLines(paste0("Dimnames in covariate field(s) ",
                                 paste(idx, collapse = ","),
                                 " are missing! This covariate field is removed."))
      x <- x[-idx]
    }
    if (length(x) == 0) return(NULL)

    ## dimnames not expected format or length
    idx <- which(sapply(x, function(x) {
      !all(c(length(dimnames(x)[[1]]) == nrow(x),
             all(!is.na(dimnames(x)[[1]])),
             length(dimnames(x)[[2]]) == ncol(x),
             all(!is.na(dimnames(x)[[2]])),
             length(dimnames(x)[[3]]) == dim(x)[3],
             all(!is.na(dimnames(x)[[3]]))
             ))
    }))
    if (length(idx) > 0) {
      if (verbose) writeLines(paste0("Dimnames in covariate field(s) ",
                                 paste(idx, collapse = ","),
                                 " do not have the expected format/length! This covariate field is removed."))
      x <- x[-idx]
    }
    if (length(x) == 0) return(NULL)


    return(x)
  }
}



##' Fill missing covariate cells next to data
##'
##' @description
##' Fill `NA` cells of a covariate (land, cloud or ice holes) ring by ring from
##' the neighbouring non-missing cells. By default only the first ring is
##' filled: the cells that touch a non-missing cell, diagonally included, e.g.
##' the coastline but not the interior of land.
##'
##' @param x A covariate: an `admove_cov` object (see [prep_cov()]), a list of
##'   them (`admove_cov_list`), or a plain `[x, y, time]` array or `[x, y]`
##'   matrix.
##' @param n_rings Number of rings to fill: a whole number of at least zero, or
##'   `Inf` to fill every cell that can be reached from data. Default `1`. `0`
##'   returns `x` unchanged.
##' @param sd Bandwidth of the Gaussian weights in cells. Default `1`. Each
##'   cell of a ring gets the weighted mean of the non-missing cells within
##'   `max(1, ceiling(2 * sd))` cells, where cells filled in earlier rings count
##'   as data.
##' @param verbose Logical; if `TRUE` (default), report how many cells were
##'   filled.
##'
##' @details
##' Covariates are interpolated bilinearly in the likelihood, which gives `NaN`
##' where the surrounding cells are missing. Filling one ring has two effects:
##' \itemize{
##'   \item [setup_data()] keeps grid cells whose centre lies in the ring. Those
##'     are mostly coastal cells that are part water, part land, which the CTMC
##'     would otherwise lose.
##'   \item The covariate and its gradient become finite up to about one cell
##'     into the gap, so a predicted mean of the Kalman filter can reach the
##'     coast. A mean further inside a gap still gives `NaN`; fill more rings
##'     for that (e.g. for conventional tags, whose mean is not updated between
##'     release and recapture and can drift far).
##' }
##' Filled values are extrapolated from the water side, so they change the
##' covariate gradients at the coast and thus the fit. The further a ring
##' reaches into a gap, the less its values mean, which is why only one ring is
##' filled by default.
##'
##' Distances use the cell sizes in the dimension names, so the weights also
##' suit cells that are not square. Each time slice is filled separately; a
##' slice without any data stays missing.
##'
##' @return `x` with the filled cells, keeping its class and attributes. For an
##'   array or `admove_cov`, the attribute `"filled"` is a logical array marking
##'   the filled cells (accumulated over repeated calls).
##'
##' @examples
##' cov <- skjepo$sim$cov
##' cov1 <- fill_cov(cov)
##' c(sum(is.na(cov)), sum(is.na(cov1)))
##'
##' ## fill every reachable cell, with smoother values
##' cov_all <- fill_cov(cov, n_rings = Inf, sd = 2)
##'
##' @seealso [setup_data()] (argument `fill_na`)
##'
##' @export
fill_cov <- function(x, n_rings = 1, sd = 1, verbose = TRUE) {

  n_rings <- .check_n_rings(n_rings, "n_rings")
  if (!is.numeric(sd) || length(sd) != 1L || !is.finite(sd) || sd <= 0) {
    stop("'sd' must be a single positive number.", call. = FALSE)
  }
  if (n_rings == 0) return(x)

  if (is.list(x) && !is.array(x)) {
    nms <- names(x)
    for (i in seq_along(x)) {
      lab <- if (is.null(nms) || nms[i] == "") paste0("covariate ", i) else
        paste0("covariate '", nms[i], "'")
      x[[i]] <- .fill_cov_one(x[[i]], n_rings, sd, verbose, lab)
    }
    return(x)
  }

  .fill_cov_one(x, n_rings, sd, verbose, "covariate")
}



## Internal functions ---------------------------------------------------------------

## Validate a number of rings for fill_cov(): a whole number >= 0 or Inf.
.check_n_rings <- function(n, name) {

  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 0 ||
        (is.finite(n) && n != round(n))) {
    stop("'", name, "' must be a single whole number of at least 0, or Inf.",
         call. = FALSE)
  }

  n
}


## Ring-wise fill of one covariate array; see fill_cov() and
## dev/code_notes.org, "Filling covariate gaps". All time slices are shifted at
## once, so the cost is one pass over the kernel offsets per ring. The window
## must contain the 8 neighbours (w >= 1), or a ring cell could get den = 0.
.fill_cov_one <- function(x, n_rings, sd, verbose, lab) {

  d <- dim(x)
  if (is.null(d) || !(length(d) %in% 2:3)) {
    stop("Each covariate must be an [x, y, time] array or an [x, y] matrix.",
         call. = FALSE)
  }
  nt <- if (length(d) == 3) d[3] else 1L
  A <- array(as.numeric(unclass(x)), c(d[1], d[2], nt))

  ## cell aspect ratio from the cell centres, for the distances in the weights
  dn <- dimnames(x)
  ry <- 1
  if (!is.null(dn) && length(dn[[1]]) > 1 && length(dn[[2]]) > 1) {
    dx <- abs(mean(diff(suppressWarnings(as.numeric(dn[[1]])))))
    dy <- abs(mean(diff(suppressWarnings(as.numeric(dn[[2]])))))
    if (is.finite(dx) && is.finite(dy) && dx > 0 && dy > 0) ry <- dy / dx
  }

  w <- max(1L, as.integer(ceiling(2 * sd)))
  offs <- expand.grid(di = -w:w, dj = -w:w)
  offs <- offs[offs$di != 0 | offs$dj != 0, ]
  offs$wt <- exp(-0.5 * (offs$di^2 + (offs$dj * ry)^2) / sd^2)
  queen <- abs(offs$di) <= 1 & abs(offs$dj) <= 1

  ## out[i, j, ] = A[i + di, j + dj, ], NA beyond the edge
  shift <- function(A, di, dj) {
    out <- array(NA_real_, dim(A))
    si <- seq_len(dim(A)[1]) + di
    sj <- seq_len(dim(A)[2]) + dj
    vi <- si >= 1 & si <= dim(A)[1]
    vj <- sj >= 1 & sj <= dim(A)[2]
    out[which(vi), which(vj), ] <- A[si[vi], sj[vj], , drop = FALSE]
    out
  }

  n_na0 <- sum(is.na(A))
  filled <- array(FALSE, dim(A))
  k <- 0
  while (k < n_rings) {
    na <- is.na(A)
    if (!any(na)) break
    touch <- array(FALSE, dim(A))
    for (o in which(queen)) {
      touch <- touch | !is.na(shift(A, offs$di[o], offs$dj[o]))
    }
    ring <- na & touch
    if (!any(ring)) break
    ## weighted means at the ring cells only: a full-array pass per offset
    ## made n_rings = Inf with a wide kernel take minutes
    ix <- which(ring, arr.ind = TRUE)
    num <- numeric(nrow(ix))
    den <- numeric(nrow(ix))
    for (o in seq_len(nrow(offs))) {
      ii <- ix[, 1] + offs$di[o]
      jj <- ix[, 2] + offs$dj[o]
      inside <- ii >= 1 & ii <= dim(A)[1] & jj >= 1 & jj <= dim(A)[2]
      v <- rep(NA_real_, nrow(ix))
      v[inside] <- A[cbind(ii[inside], jj[inside], ix[inside, 3])]
      ok <- !is.na(v)
      num[ok] <- num[ok] + offs$wt[o] * v[ok]
      den[ok] <- den[ok] + offs$wt[o]
    }
    A[ix] <- num / den
    filled <- filled | ring
    k <- k + 1
  }

  if (verbose) {
    empty <- sum(apply(is.na(A), 3, all))
    message(lab, ": filled ", sum(filled), " of ", n_na0, " NA cells (",
            k, " ring", if (k == 1) "" else "s", "); ", sum(is.na(A)),
            " remain NA",
            if (empty > 0) paste0(", ", empty, " time slice(s) without any data"),
            ".")
  }

  x[] <- A
  dim(filled) <- d
  dimnames(filled) <- dimnames(x)
  prev <- attr(x, "filled")
  if (!is.null(prev) && identical(dim(prev), d)) filled <- filled | prev
  attr(x, "filled") <- filled
  x
}

.get_cov_trange <- function(cov) {

  cov <- .make_cov_list(cov)

  res <- sapply(cov, function(x) range(as.numeric(dimnames(x)[[3]])))

  return(res)
}


.get_cov_xyrange <- function(x) {

  ncov <- length(x)
  xr <- yr <- matrix(NA, ncov, 2)
  for (i in 1:ncov) {
    xgr <- as.numeric(as.character(dimnames(x[[i]])[[1]]))
    xr[i,] <- range(xgr)
    ygr <- as.numeric(as.character(dimnames(x[[i]])[[2]]))
    yr[i,] <- range(ygr)
  }

  return(list(xr = xr,
              yr = yr))
}


.validate_dimnames <- function(x1, dims, x2 = NULL, x3 = NULL,
                               strict = TRUE,
                               digits = 6,
                               verbose = TRUE) {

  x_list <- list(x1, x2, x3)
  dimnams_out <- vector("list", 3)
  err <- NULL
  ## x, y, t
  for (i in 1:3) {
    ## inputs
    for (j in 1:3) {
      if (length(x_list[[j]]) >= j &&
            !is.null(x_list[[j]][[i]]) &&
            !all(is.na(x_list[[j]][[i]]))) {
        if (length(x_list[[j]][[i]]) != dims[i]) {
          stop(paste0("The length of the specified dimension names in ",
                      c("'x_centers'","'y_centers'","'time'")[i],
                      " does not match the dimensions of the data! Please check!"))
        }
        d_char <- as.character(x_list[[j]][[i]])
        if (!all(is.na(d_char))) {
          dimnams_out[[i]] <- d_char
          break
        } else if (j < 3) {
          next
        } else {
          err <- c(err, i)
          if (!strict) {
            dxyt <- 1/dims[i]
            if (i == 3) {
              d_dummy <- seq(0, dims[i] - 1, length.out = dxyt)
            } else {
              d_dummy <- seq(dxyt/2, 1 - dxyt/2, dxyt)
            }
            dimnams_out[[i]] <- round(d_dummy, digits = digits)
            break
          }
        }
      } else if (j < 3) {
        next
      } else {
        err <- c(err, i)
        if (!strict) {
          dxyt <- 1/dims[i]
          if (i == 3) {
            d_dummy <- seq(0, dims[i] - 1, length.out = dxyt)
          } else {
            d_dummy <- seq(dxyt/2, 1 - dxyt/2, dxyt)
          }
          dimnams_out[[i]] <- round(d_dummy, digits = digits)
          break
        }
      }
    }
  }
  err <- unique(err)
  if (!is.null(err) && length(err) > 0) {
    if (strict) {
      stop(paste0("No valid dimension names for ",
                  paste(c("x coordinates",
                          "y coordinates",
                          "time")[err], collapse = ", "),
                  "! Use the ",
                  paste(c("'x_centers'","'y_centers'","'time'")[err], collapse = ", "),
                  " argument to specify the missing dimension!"))
    } else {
      if (verbose) {
        message(paste0("No valid dimension names for ",
                       paste(c("x coordinates",
                               "y coordinates",
                               "time")[err], collapse = ", "),
                       "! Use the ",
                       paste(c("'x_centers'","'y_centers'","'time'")[err],
                             collapse = ", "),
                       " argument to specify the missing dimension! Using some dummy defaults, which might create a mismatch with grid or tags later!"))
      }
    }
  }

  return(dimnams_out)
}



.make_cov_list <- function(x, verbose = FALSE){

  if (is.null(x) || inherits(x, "admove_cov_list")) return(x)


  if (is.list(x) && all(sapply(x, function(x) inherits(x, "admove_cov")))) {
    x <- .add_class(x, "admove_cov_list")
    return(x)
  }

  .check_class(x, "admove_cov")

  x0 <- x
  x <- list(x)
  attributes(x[[1]]) <- attributes(x0)
  x <- .add_class(x, "admove_cov_list")
  x <- add_sref(x, sref(x0), verbose = verbose)
  x <- add_tref(x, tref(x0), verbose = verbose)
  x
}



## s3 methods ----------------------------------------------------------------------

##' Subset an `admove_cov` object
##'
##' @description
##' Subsetting method for objects of class `admove_cov`.
##'
##' This method preserves the `admove_cov` class and associated attributes when
##' the result remains a three-dimensional covariate array. If subsetting
##' returns an object with fewer than three dimensions, the result is returned as
##' a regular R object without the `admove_cov` class.
##'
##' In particular, subsetting only along the third dimension (time) returns a
##' subsetted `admove_cov` object with the corresponding time-related attributes
##' updated where available.
##'
##' @param x An object of class `admove_cov`.
##' @param i Indices for the first dimension.
##' @param j Indices for the second dimension.
##' @param k Indices for the third dimension, typically corresponding to time.
##' @param ... Further indices passed to the underlying array subsetting
##'   operation.
##' @param drop Logical; should dimensions of length one be dropped? Default:
##'   `TRUE`.
##'
##' @return
##' An object of class `admove_cov` if the subset retains three dimensions;
##' otherwise a regular subsetted R object.
##'
##' @details
##' The method first performs standard array subsetting on the unclassed object.
##' If the resulting object still has three dimensions, attributes other than
##' `dim` and `dimnames` are restored and the original class is reattached.
##'
##' @name subset-admove_cov
##' @export
`[.admove_cov` <- function(x, i, j, k, ..., drop = TRUE) {

  ux <- unclass(x)

  if (!missing(i) && missing(j) && missing(k) && length(list(...)) == 0) {
    y <- ux[i]
    return(y)
  }

  if (missing(i) && missing(j) && !missing(k) && length(list(...)) == 0) {
    y <- ux[,,k, drop = drop]

    d <- dim(y)
    if (!is.null(d) && length(d) == 3) {
      ax <- attributes(x)

      if (!is.null(ax$time)) ax$time <- ax$time[k]

      keep <- setdiff(names(ax), c("dim", "dimnames"))
      for (nm in keep) attr(y, nm) <- ax[[nm]]
      class(y) <- class(x)
    }

    return(y)
  }

  mc <- match.call(expand.dots = TRUE)
  mc[[1]] <- base::`[`
  mc$x <- quote(ux)

  y <- eval(mc, envir = list(ux = ux), enclos = parent.frame())

  d <- dim(y)
  if (is.null(d) || length(d) != 3) {
    return(y)
  }

  ax <- attributes(x)

  if (!missing(k) && !is.null(ax$time)) ax$time <- ax$time[k]

  keep <- setdiff(names(ax), c("dim", "dimnames"))
  for (nm in keep) attr(y, nm) <- ax[[nm]]

  class(y) <- class(x)
  y
}


##' @method summary admove_cov
##' @rdname summarise_cov
##' @export
summary.admove_cov <- function(object, ...) {
  summarise_cov(object, ...)
}


##' @method summary admove_cov
##' @rdname summarise_cov
##' @export
summary.admove_cov_list <- function(object, ...) {
  summarise_cov(object, ...)
}


##' @rdname print-admove
##' @method print admove_cov
##' @export
print.admove_cov <- function(x, ...) {
  tmp <- x
  attributes(tmp) <- NULL
  NextMethod("print", tmp, ...)
}


##' @rdname plot_cov
##' @export
plot.admove_cov <- function(x, ...) {
  plot_cov(x, ...)
  return(invisible(NULL))
}


##' @rdname plot_cov
##' @export
plot.admove_cov_list <- function(x, ...) {
  plot_cov(x, ...)
  return(invisible(NULL))
}



##' @rdname sref
##' @export
sref.admove_cov_list <- function(x, ...) {
  sp <- attr(x, "sref")
  if (!is.null(sp)) return(sp)
  if (length(x) == 0L) stop("Object has no 'sref' attribute.")
  srefs <- lapply(x, function(co) try(sref(co), silent = TRUE))
  srefs_ok <- Filter(function(s) !inherits(s, "try-error"), srefs)
  if (length(srefs_ok) == 0L) stop("No elements have an 'sref' attribute.")
  ref <- srefs_ok[[1L]]
  differ <- any(!vapply(srefs_ok[-1L], function(s) sref_equal(ref, s), logical(1L)))
  if (differ)
    warning("Spatial references differ across covariate list elements; returning the first.")
  ref
}

##' @rdname sref-set
##' @export
`sref<-.admove_cov_list` <- function(x, value) {
  validated <- validate_sref(value)
  for (i in seq_along(x)) attr(x[[i]], "sref") <- validated
  attr(x, "sref") <- validated
  x
}

##' @rdname tref
##' @export
tref.admove_cov_list <- function(x, ...) {
  tr <- attr(x, "tref")
  if (!is.null(tr)) return(tr)
  if (length(x) == 0L) stop("Object has no 'tref' attribute.")
  trefs <- lapply(x, function(co) try(tref(co), silent = TRUE))
  trefs_ok <- Filter(function(tr) !inherits(tr, "try-error"), trefs)
  if (length(trefs_ok) == 0L) stop("No elements have a 'tref' attribute.")
  ref <- trefs_ok[[1L]]
  differ <- any(!vapply(trefs_ok[-1L], function(t)
    identical(t$units, ref$units) && identical(t$period, ref$period), logical(1L)))
  if (differ)
    warning("Time references differ across covariate list elements; returning the first.")
  ref
}

##' @rdname tref-set
##' @export
`tref<-.admove_cov_list` <- function(x, value) {
  for (i in seq_along(x)) x[[i]] <- add_tref(x[[i]], value, verbose = FALSE)
  attr(x, "tref") <- attr(x[[1L]], "tref")
  x
}
