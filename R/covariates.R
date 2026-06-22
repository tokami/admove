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
##' @param date_format Optional character string passed to [as.Date()] to parse
##'   character time labels. Default: `NULL`.
##' @param date_origin Optional origin passed to [as.Date()] when time labels are
##'   stored numerically. Default: `NULL`.
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

    if (is.null(date_format) && !is.null(date_origin)) {
      dati <- as.Date(dati, origin = date_origin, tz = tz)
    } else if (is.null(date_origin) && !is.null(date_format)) {
      dati <- as.Date(dati, format = date_format, tz = tz)
    } else if (isFALSE(date_decimal)) {
      dati <- as.Date(dati, format = date_format, origin = date_origin,
                      tz = tz)
    } else {
      dati <- .decimal_year_2_date(as.numeric(dati), tz = tz)
    }

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
##' Each time slice is plotted as a separate panel, optionally with contour
##' lines and land added to the plot.
##'
##' @param x An object of class `admove_cov`, or an object containing covariate
##'   data such as `admove_data`, `admove_sim`, or `admove`.
##' @param i Index (scalar or vector) used when `x` is an `admove_cov_list`.
##'   A scalar selects one covariate and plots all its time steps. A vector
##'   selects multiple covariates and produces one panel per element (showing
##'   the first time step, or the first element of `select`). Default: `1`.
##' @param select Optional vector of time-step indices to plot. Default:
##'   `NULL`, in which case all time steps are plotted.
##' @param main Main title of the plot. Default: `"Covariate fields"`.
##' @param labels Logical; currently reserved for plotting cell labels. Default:
##'   `TRUE`.
##' @param plot_land Logical; if `TRUE`, add land to the plot. Default:
##'   `FALSE`.
##' @param auto_layout Logical; if `TRUE`, plotting parameters are set
##'   automatically and restored afterwards. Default: `TRUE`.
##' @param xlab Label for the x-axis. Default: `"x"`.
##' @param ylab Label for the y-axis. Default: `"y"`.
##' @param bg Optional background colour for the plot. Default: `NULL`.
##' @param plot_contour Logical; if `TRUE`, add contour lines. Default:
##'   `TRUE`.
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
##' @name plot_cov
##' @export
plot_cov <- function(x,
                     i = 1,
                     select = NULL,
                     main = "Covariate fields",
                     labels = TRUE,
                     plot_land = FALSE,
                     auto_layout = TRUE,
                     xlab = "x",
                     ylab = "y",
                     bg = NULL,
                     plot_contour = TRUE,
                     ...) {

  if (inherits(x, "admove_sim")) {
    cov <- x$cov
  } else if(inherits(x, "admove_data")) {
    cov <- x$cov
  } else if(inherits(x, "admove")) {
    cov <- x$dat$cov
  } else{
    cov <- x
  }

  if (inherits(cov, "admove_cov_list") && length(i) > 1) {
    sel <- cov[i]
    n <- length(sel)
    nms <- names(sel)
    t_sel <- if (is.null(select)) 1L else select[1L]
    if (auto_layout) {
      opar <- par(no.readonly = TRUE)
      on.exit(par(opar))
      par(mfrow = n2mfrow(n, asp = 2),
          mar = c(1.5, 1.5, 1.5, 1.5),
          oma = c(3, 3, ifelse(main == "", 0, 1.5), 0))
    }
    for (j in seq_along(sel)) {
      plot_cov(sel[[j]], select = t_sel, main = "",
               plot_land = plot_land, auto_layout = FALSE,
               xlab = xlab, ylab = ylab, bg = bg,
               plot_contour = plot_contour, ...)
      panel_lbl <- if (!is.null(nms) && nzchar(nms[j])) nms[j] else paste0("Layer ", i[j])
      legend("topleft", legend = panel_lbl, bg = "white", pch = NA)
    }
    if (auto_layout) {
      mtext(main, 3, 0, outer = TRUE)
      mtext(xlab, 1, 1, outer = TRUE)
      mtext(ylab, 2, 1.5, outer = TRUE)
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

  if(any(names(attributes(cov)) == "dimnames")){
    xlims <- range(as.numeric(attributes(cov)$dimnames[[1]]))
    ylims <- range(as.numeric(attributes(cov)$dimnames[[2]]))
  }else{
    xlims <- c(0,1)
    ylims <- c(0,1)
  }

  if(auto_layout){
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar))
    par(mfrow = n2mfrow(nt, asp = 2),
        mar = c(1.5,1.5,1.5,1.5),
        oma = c(3,3,ifelse(main == "", 0, 1.5),0))
  }
  for(i in 1:nt){
    x <- as.numeric(rownames(cov[,,i]))
    if(length(x) == 0) x <- 1:nrow(cov[,,i])
    y <- as.numeric(colnames(cov[,,i]))
    if(length(y) == 0) y <- 1:ncol(cov[,,i])
    if(!is.null(bg)){
      par(bg = bg)
    }
    plot(1,1, type = "n",
         xlim = xlims, ylim = ylims,
         xlab = "",
         ylab = "",
         asp = 1,
         ...)
    ## if(!is.null(bg)){
    ##     usr <- par("usr")
    ##     rect(usr[1], usr[3], usr[2], usr[4], col = bg, border = NA)
    ## }
    image(x, y, cov[,,i, drop = TRUE],
          col = terrain.colors(100),
          add = TRUE)
    if (plot_land) {
      plot_land(sref)
    }
    if(plot_contour) contour(x, y, cov[,,i], add = TRUE)
    if(nt > 1) legend("topleft", legend = paste0("Field ", i),
                      bg = "white", pch = NA)
    box(lwd = 1.5)
  }
  if(auto_layout){
    mtext(main, 3, 0, outer = TRUE)
    mtext(xlab, 1, 1, outer = TRUE)
    mtext(ylab, 2, 1.5, outer = TRUE)
  }


  return(invisible(NULL))
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



## Internal functions ---------------------------------------------------------------

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
