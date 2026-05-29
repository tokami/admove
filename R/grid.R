
## Main functions --------------------------------------------------------------------

##' Create or modify a spatial grid for admove
##'
##' @description
##' Creates a new spatial grid or derives one from an existing object for use in
##' \code{admove}. Grids are used for spatial prediction and are required for
##' model formulations based on a discrete spatial domain, such as the
##' continuous-time Markov chain (CTMC) approach.
##'
##' @param x An optional object from which to derive the grid. Supported inputs
##'   include objects of class \code{"admove_grid"}, \code{"admove_cov"},
##'   \code{"admove_tags"}, \code{"admove_data"}, \code{"admove_sim"},
##'   \code{"admove"}, as well as \code{sf}/\code{sfc} objects, raster objects,
##'   and matrices. If \code{NULL}, a new grid is created from the supplied
##'   \code{xrange}, \code{yrange}, and \code{cellsize}.
##' @param cellsize Numeric vector giving the grid cell size in the x- and
##'   y-direction. If a single value is supplied, it is used for both
##'   directions. If \code{NULL}, a default value is derived from the spatial
##'   extent.
##' @param xrange Numeric vector of length 2 giving the range of the x
##'   dimension of the spatial domain. Ignored if extracted from \code{x}.
##' @param yrange Numeric vector of length 2 giving the range of the y
##'   dimension of the spatial domain. Ignored if extracted from \code{x}.
##' @param select Controls which cells are retained in the grid. If
##'   \code{FALSE} (default), all eligible cells are kept. If \code{TRUE} or an
##'   integer value, cells can be selected interactively. If a numeric vector of
##'   length greater than 1 is supplied, it is interpreted as cell indices to
##'   retain; if all supplied indices are negative, those cells are removed.
##' @param crs Optional coordinate reference system for the grid.
##' @param units Optional spatial units for the grid, for example
##'   \code{"degree"}, \code{"m"}, or \code{"km"}.
##' @param crs_scale Optional scaling between CRS units and the numeric units
##'   used in the grid.
##' @param plot_land Logical; if \code{TRUE}, land masses are added to plots,
##'   where supported.
##' @param auto_layout Logical; if \code{TRUE}, plotting methods may adjust
##'   graphical parameters automatically.
##' @param plot Logical; if \code{TRUE}, the resulting grid is plotted.
##' @param force Logical; if \code{TRUE}, the grid is created even when the
##'   implied number of grid breaks in x or y exceeds 1000. Use with care.
##' @param verbose Logical; if \code{TRUE}, informative messages are printed.
##'
##' @details
##' The function can construct a grid from scratch, inherit the extent and
##' missing-cell structure from an existing \code{admove} object, or derive the
##' extent from spatial objects such as \code{sf} geometries, rasters, or
##' covariate arrays.
##'
##' When \code{x} contains missing spatial cells, these are propagated to the
##' new grid. For \code{sf} objects, cells can be restricted to locations inside
##' the supplied geometry. Optional interactive selection can be used to include
##' or exclude cells manually.
##'
##' The returned object stores the grid-cell centres, cell indices, lookup
##' table, grid breaks, and spatial extent, together with the associated spatial
##' reference.
##'
##' @return
##' An object of class \code{"admove_grid"} describing the spatial grid.
##'
##' @examples
##' grid <- create_grid()
##' dim(grid)
##' bbox_grid(grid)
##'
##' @export
create_grid <- function(x = NULL,
                        cellsize = NULL,
                        xrange = NULL,
                        yrange = NULL,
                        select = FALSE,
                        crs = NULL,
                        units = NULL,
                        crs_scale = NULL,
                        plot_land = FALSE,
                        auto_layout = TRUE,
                        plot = FALSE,
                        force = FALSE,
                        verbose = TRUE) {

  cellsize0 <- cellsize
  xrange0 <- xrange
  yrange0 <- yrange
  use_sf <- FALSE
  use_raster <- FALSE
  use_grid <- FALSE
  is_na <- NULL
  crs0 <- crs
  units0 <- units
  crs_scale0 <- crs_scale
  xcen_cov <- NULL
  ycen_cov <- NULL

  if (!is.null(crs) && !all(is.na(crs))) {
    units_try <- try(crs$units_gdal, silent = TRUE)
    if (!inherits(units_try, "try-error")) {
      units <- units_try
    } ## else if (verbose) message("Both crs and units specified! They should align.")
  }

  sref <- create_sref(crs, units, 1)
  crs_scale <- sref$crs_scale

  ## helper: map points -> NA-grid cell indices
  cell_index <- function(x, y, x0, y0, cellsize) {
    ## small epsilon to reduce boundary issues from floating point noise
    eps <- 1e-12
    ix <- floor((x - x0 + eps) / cellsize[1])
    iy <- floor((y - y0 + eps) / cellsize[2])
    paste(ix, iy, sep = "_")
  }

  ## default if not provided
  if (is.null(xrange0)) xrange <- c(0,1)
  if (is.null(yrange0)) yrange <- c(0,1)

  ## Checks
  if (is.null(x) && length(xrange) == 1) stop("Only a single value provided for xrange. Please provide upper and lower range!")
  if (is.null(x) && length(yrange) == 1) stop("Only a single value provided for yrange. Please provide upper and lower range!")

  xrange <- sort(xrange)
  if (abs(xrange[2] - xrange[1]) < 1e-10) stop("Are the xrange values correct? They are basically identical!")
  yrange <- sort(yrange)
  if (abs(yrange[2] - yrange[1]) < 1e-10) stop("Are the yrange values correct? They are basically identical!")

  ## default if not provided
  if (is.null(cellsize0)) cellsize <- c(diff(xrange) / 10, diff(yrange) / 10)

  if (length(cellsize) == 1) {
    cellsize <- rep(cellsize, 2)
  } else if(cellsize[1] != cellsize[2]) {
    if (verbose) message(paste0("Non regular square cells are not yet implemented in admove. ", min(cellsize), " is used for the cell size in both x and y direction."))
    cellsize <- rep(min(cellsize), 2)
  }


  if (!is.null(x) && (inherits(x, "admove") ||
                        inherits(x, "admove_sim") ||
                        inherits(x, "admove_data"))) {

    if (inherits(x, "admove") || inherits(x, "admove_sim")) {
      dat <- x$dat
    } else {
      dat <- x
    }

    cellsize_na <- dat$grid$cellsize
    xgr <- dat$grid$xgr
    xcen <- xgr[-1] - 0.5 * cellsize_na[1]
    ygr <- dat$grid$ygr
    ycen <- ygr[-1] - 0.5 * cellsize_na[2]

    tmp <- is.na(dat$grid$celltable)
    tmp <- which(tmp, arr.ind = TRUE)

    is_na <- data.frame(x = xcen[tmp[,1]],
                        y = ycen[tmp[,2]])

    if (is.null(cellsize0)) cellsize <- cellsize_na
    if (is.null(xrange0)) xrange <- dat$grid$xrange
    if (is.null(yrange0)) yrange <- dat$grid$yrange

    sref <- sref(dat$grid)

  } else if (!is.null(x) && (inherits(x, "admove_cov") ||
                               inherits(x, "admove_cov_list") ||
                               inherits(x, "matrix"))) {

    if (inherits(x, "admove_cov_list")) {
      xcens <- lapply(x, function(co) as.numeric(dimnames(co)[[1]]))
      ycens <- lapply(x, function(co) as.numeric(dimnames(co)[[2]]))
      xcen <- sort(Reduce(intersect, xcens))
      ycen <- sort(Reduce(intersect, ycens))
      if (length(xcen) == 0L || length(ycen) == 0L)
        stop("Covariate fields have no common spatial grid cells.")
      if (verbose && (length(xcen) < length(xcens[[1L]]) ||
                        length(ycen) < length(ycens[[1L]])))
        message("Covariate fields have different spatial extents; using the intersection.")
      na_union <- matrix(FALSE, length(xcen), length(ycen))
      for (co in x) {
        xi <- match(xcen, as.numeric(dimnames(co)[[1]]))
        yi <- match(ycen, as.numeric(dimnames(co)[[2]]))
        na_co <- apply(is.na(unclass(co)[xi, yi, , drop = FALSE]), c(1L, 2L), any)
        na_union <- na_union | na_co
      }
      tmp <- which(na_union, arr.ind = TRUE)
    } else {
      xcen <- as.numeric(dimnames(x)[[1]])
      ycen <- as.numeric(dimnames(x)[[2]])
      tmp <- apply(unclass(x), c(1L, 2L), function(co) any(is.na(co)))
      tmp <- which(tmp, arr.ind = TRUE)
    }

    cellsize_na <- c(min(diff(xcen)), min(diff(ycen)))
    if (is.null(cellsize0)) {
      cellsize <- cellsize_na
      xcen_cov <- xcen
      ycen_cov <- ycen
    }
    xrange_cov <- range(xcen - cellsize_na[1]/2, xcen + cellsize_na[1]/2)
    yrange_cov <- range(ycen - cellsize_na[2]/2, ycen + cellsize_na[2]/2)
    if (is.null(xrange0)) xrange <- xrange_cov
    if (is.null(yrange0)) yrange <- yrange_cov

    is_na <- data.frame(x = xcen[tmp[, 1]],
                        y = ycen[tmp[, 2]])

    if (inherits(x, "admove_cov") || inherits(x, "admove_cov_list")) {
      sref <- sref(x)
    }

  } else if (!is.null(x) && inherits(x, "admove_grid")) {

    if (is.null(xrange0)) xrange <- x$xrange
    if (is.null(yrange0)) yrange <- x$yrange
    idx <- which(is.na(x$celltable), arr.ind = TRUE)

    if (is.null(cellsize0)) cellsize <- x$cellsize

    use_grid <- TRUE

    sref <- sref(x)

  } else if (!is.null(x) && (inherits(x, "admove_tags"))) {

    xrange <- range(x[,"x"])
    yrange <- range(x[,"y"])

    sref <- sref(x)

    if (is.null(cellsize0)) cellsize <- c(diff(xrange) / 10, diff(yrange) / 10)

  } else if (!is.null(x) && (inherits(x, "sf") || inherits(x, "sfc"))) {

    if (!requireNamespace("sf", quietly = TRUE)) {
      stop("Package 'sf' is required to extract information from sf objects. Please install it or provide information to create_grid and set x to NULL.")
    }

    ## dimensions
    bb <- sf::st_bbox(x)
    xrange <- unname(c(bb["xmin"], bb["xmax"]))
    yrange <- unname(c(bb["ymin"], bb["ymax"]))

    if (is.null(cellsize0)) cellsize <- c(diff(xrange) / 10, diff(yrange) / 10)

    use_sf <- TRUE

    crs <- sf::st_crs(x)
    sref <- create_sref(crs = crs,
                          units = crs$units_gdal,
                          crs_scale = 1)

  } else if (!is.null(x) && inherits(x, "RasterLayer")) {

    if (!requireNamespace("raster", quietly = TRUE)) {
      stop("Package 'raster' is required to extract information from sf objects. Please install it or provide information to create_grid and set x to NULL.")
    }

    ## dimensions
    bb <- raster::bbox(x)
    xrange <- c(bb[1,1], bb[1,2])
    yrange <- c(bb[2,1], bb[2,2])

    if (is.null(cellsize0)) cellsize <- c(diff(xrange) / 10, diff(yrange) / 10)

    use_raster <- TRUE

    p4 <- as.character(raster::crs(x))

    if (!requireNamespace("sf", quietly = TRUE)) {
      stop("Package 'sf' is required to extract information from sf objects. Please install it or provide information to create_grid and set x to NULL.")
    }

    sref <- list()
    sref$crs <- sf::st_crs(pf)
    sref$units <- sref$crs$units_gdal
    sref$crs_scale <- 1

  }

  tmp <- round(diff(xrange) / cellsize[1])
  if (tmp > 1e3 && !force) stop(paste0("Are you sure your settings are correct? The provided cell size (cellsize) implies ", tmp, " x breaks. That is too large. Check the cellsize values (units) or use force = TRUE if you are sure about what you are doing."))

  if (!is.null(xcen_cov)) {
    xcen <- xcen_cov
    xgr <- c(xcen - cellsize[1] / 2, xcen[length(xcen)] + cellsize[1] / 2)
  } else {
    xgr <- seq(xrange[1], xrange[2], by = cellsize[1])
    xcen <- xgr[-1] - 0.5 * cellsize[1]
  }

  tmp <- round(diff(yrange) / cellsize[2])
  if (tmp > 1e3 && !force) stop(paste0("Are you sure your settings are correct? The provided cell size (cellsize) implies ", tmp, " y breaks. That is too large. Check the cellsize values (units) or use force = TRUE if you are sure about what you are doing."))

  if (!is.null(ycen_cov)) {
    ycen <- ycen_cov
    ygr <- c(ycen - cellsize[2] / 2, ycen[length(ycen)] + cellsize[2] / 2)
  } else {
    ygr <- seq(yrange[1], yrange[2], by = cellsize[2])
    ycen <- ygr[-1] - 0.5 * cellsize[2]
  }

  xygrid <- expand.grid(x = xcen, y = ycen)
  igrid <- expand.grid(idx = seq_along(xcen), idy = seq_along(ycen))
  celltable <- matrix(NA,
                      nrow = length(xcen),
                      ncol = length(ycen))

  if (length(select) > 1) {

    idx <- sort(unique(select))
    xygrid <- xygrid[idx,]
    igrid <- igrid[idx,]

  } else if(as.integer(select) > 0) {
    opts <- options()
    on.exit(options(opts))
    options(locatorBell = FALSE)

    if(auto_layout){
      opar <- par(no.readonly = TRUE)
      on.exit(par(opar))
      par(mfrow = c(1,1))
    }

    if(verbose) message("Please refer to the graphical device and select the cells to include/exclude in the grid by clicking onto the map. Right-click to exit the selection function_")

    if(as.integer(select) == 1){
      idx <- NULL
    }else if(as.integer(select) == 2){
      idx <- 1:(length(xcen) * length(ycen))
    }

    if (use_sf) {
      pts <- sf::st_as_sf(xygrid, coords = c("x", "y"), crs = sf::st_crs(x))
      inside <- sf::st_within(pts, x, sparse = FALSE)[, 1]
      if (as.integer(select) == 1) {
        idx <- which(inside)
      } else if(as.integer(select) == 2){
        idx <- which(!inside)
      }
    }

    if (!is.null(is_na) && nrow(is_na) > 0) {
      x0 <- min(is_na$x) - cellsize_na[1]/2
      y0 <- min(is_na$y) - cellsize_na[2]/2
      na_key <- cell_index(is_na$x, is_na$y, x0, y0, cellsize_na)
      xy_key <- cell_index(xygrid$x, xygrid$y, x0, y0, cellsize_na)
      keep <- !(xy_key %in% na_key)

      if (as.integer(select) == 1) {
        idx <- which(keep)
      } else if(as.integer(select) == 2){
        idx <- which(!keep)
      }
    }


    repeat{
      plot(xygrid[,1], xygrid[,2],
           type = "n",
           main = "Click on points to select/de-select,\nright-click to exit",
           xlim = xrange,
           ylim = yrange,
           xaxs = "i", yaxs = "i",
           xlab = "x", ylab = "y")
      c0 <- matrix(NA,
                   nrow = length(xcen),
                   ncol = length(ycen))
      if (length(idx) > 0) {
        sel <- rep(NA,nrow(igrid))
        sel[idx] <- 1
        c0[cbind(igrid$idx, igrid$idy)] <- sel
      }
      image(xgr, ygr,
            c0,
            col = adjustcolor("dodgerblue2", 0.2),
            xlim = xrange,
            ylim = yrange,
            add = TRUE)
      if (plot_land) {
        plot_land(sref = sref)
      }
      if (use_sf) {
        plot(sf::st_geometry(x), add = TRUE, border = "goldenrod2", lwd = 2)
      }
      abline(v = xgr)
      abline(h = ygr)
      points(xygrid[,1], xygrid[,2], col = "dodgerblue2")
      points(xygrid[idx,1], xygrid[idx,2], pch = 16,
             col = "dodgerblue2")
      idx_sel <- identify(xygrid[,1], xygrid[,2], n = 1)
      idx <- if(any(idx == idx_sel)) idx[idx != idx_sel] else c(idx, idx_sel)

      if (length(idx_sel) == 0) {
        cat("No more points selected, exiting...\n")
        break
      }
    }

    idx <- sort(unique(idx))
    xygrid <- xygrid[idx,]
    igrid <- igrid[idx,]

  }

  if (length(select) < 2) {
    if (use_grid) {

      idx1 <- x$celltable[cbind(cut(xygrid$x, x$xgr),
                                cut(xygrid$y, x$ygr))]
      idx2 <- x$celltable[cbind(cut(xygrid$x, x$xgr,
                                    right = FALSE),
                                cut(xygrid$y, x$ygr,
                                    right = FALSE))]
      idx3 <- x$celltable[cbind(cut(xygrid$x, x$xgr),
                                cut(xygrid$y, x$ygr,
                                    right = FALSE))]
      idx4 <- x$celltable[cbind(cut(xygrid$x, x$xgr,
                                    right = FALSE),
                                cut(xygrid$y, x$ygr))]
      idx <- 1:nrow(xygrid)
      idx[apply(cbind(idx1,idx2,idx3,idx4),1,function(x) sum(is.na(x)) > 2)] <- NA
      idx <- which(!is.na(idx))
      if(length(idx) > 0){
        xygrid <- xygrid[idx,]
        igrid <- igrid[idx,]
      }

    } else if (!is.null(is_na)  && nrow(is_na) > 0) {

      x0 <- min(is_na$x) - cellsize_na[1]/2
      y0 <- min(is_na$y) - cellsize_na[2]/2
      na_key <- cell_index(is_na$x, is_na$y, x0, y0, cellsize_na)
      xy_key <- cell_index(xygrid$x, xygrid$y, x0, y0, cellsize_na)
      keep <- !(xy_key %in% na_key)

      idx <- which(keep)
      if(length(idx) > 0){
        xygrid <- xygrid[idx,]
        igrid <- igrid[idx,]
      }

    } else if (use_sf && as.integer(select) < 1) {

      pts <- sf::st_as_sf(xygrid, coords = c("x", "y"), crs = sf::st_crs(x))
      inside <- sf::st_within(pts, x, sparse = FALSE)[, 1]

      if (as.integer(select) == 0) {
        idx <- which(inside)
      } else {
        idx <- which(!inside)
      }

      if(length(idx) > 0){
        xygrid <- xygrid[idx,]
        igrid <- igrid[idx,]
      }
    }
  }

  celltable[cbind(igrid$idx, igrid$idy)] <- 1:nrow(igrid)

  res <- list(
    xygrid = xygrid,
    igrid = igrid,
    celltable = celltable,
    cellsize = cellsize,
    xgr = xgr,
    ygr = ygr,
    xrange = xrange,
    yrange = yrange
  )

  ## Return
  res <- .add_class(res, "admove_grid")
  res <- add_sref(res, sref)

  if (!is.null(units0) && units0 != units) {
    res <- scale_sref(res, units = units0)
  } else if (!is.null(crs_scale0) && crs_scale0 != crs_scale) {
    res <- scale_sref(res, scale = crs_scale0)
  }


  ## if (!is.null(crs_scale0) && crs_scale0 != crs_scale(res)) {
  ##     res <- change_units(res, crs_scale)
  ## } else if (!is.null(units0) && units0 != units(res)) {
  ##   crs_scale <- try(guess_crs_scale(crs(res), units0), silent = TRUE)
  ##   if(!inherits(crs_scale, "try-error")) {
  ##     res <- change_units(res, crs_scale)
  ##   } else if (verbose) message("Couldn't change the units.")
  ## }

  if (plot) plot(res, plot_land = plot_land, auto_layout = auto_layout)

  return(res)
}



##' Return the bounding box of an admove grid
##'
##' @description
##' Returns the spatial extent of an \code{admove_grid} object as a named vector
##' containing the minimum and maximum x and y coordinates.
##'
##' @param grid A grid object of class \code{"admove_grid"}, as returned by
##'   [create_grid()].
##'
##' @return
##' A named numeric vector with elements \code{xmin}, \code{xmax},
##' \code{ymin}, and \code{ymax}.
##'
##' @export
bbox_grid <- function(grid) {
  xr <- range(grid$xgr)
  yr <- range(grid$ygr)
  c(xmin = xr[1], xmax = xr[2], ymin = yr[1], ymax = yr[2])
}



##' Summarise an admove grid
##'
##' @description
##' Prints a short summary of an \code{admove_grid} object, including the number
##' of active cells, grid dimensions, cell size, spatial extent, number of
##' missing cells, and spatial units.
##'
##' @param object An object of class \code{"admove_grid"} or an object containing a
##'   grid, such as \code{"admove_data"}, \code{"admove_sim"}, or
##'   \code{"admove"}.
##' @param ... Additional arguments
##'
##' @details
##' If \code{x} is not itself an \code{admove_grid}, the grid is extracted from
##' the corresponding component of the supplied object.
##'
##' @return
##' Invisibly returns the extracted \code{admove_grid} object.
##'
##' @examples
##' summarise_grid(skjepo$grid)
##'
##' @name summarise_grid
##' @export
summarise_grid <- function(object, ...) {
  x <- object

  if(inherits(x, "admove_sim")) {
    grid <- x$grid
  } else if(inherits(x, "admove_data")) {
    grid <- x$grid
  } else if(inherits(x, "admove")) {
    grid <- x$dat$grid
  } else if(inherits(x, "admove_grid")) {
    grid <- x
  } else stop("Please provide an object of class 'admove_grid' or an object containing an such an object (e.g. admove_data, admove_sim, admove).")

  bb <- sprintf("%.2f", bbox_grid(grid))
  dims <- dim.admove_grid(grid)
  n_na <- sum(is.na(grid$celltable))

  units <- try(units_space(grid), silent = TRUE)
  if (is.null(units) || is.na(units) || units == "" || inherits(units, "try-error")) units <- "not specified"

  labw <- 10

  cat("<admove_grid>\n")
  cat(sprintf(paste0("  %-", labw, "s %s\n"), "cells:",
              nrow(grid$xygrid)))
  cat(sprintf(paste0("  %-", labw, "s %s\n"), "dims:",
              paste(dims[1], "x", dims[2])))
  cat(sprintf(paste0("  %-", labw, "s %s\n"), "cellsize:",
              paste(sprintf("%.2f", grid$cellsize[1]), "x",
                    sprintf("%.2f", grid$cellsize[1]))))
  cat(sprintf(paste0("  %-", labw, "s %s\n"), "xrange:",
              paste0("[",bb[1], ", ", bb[2],"]")))
  cat(sprintf(paste0("  %-", labw, "s %s\n"), "yrange:",
              paste0("[",bb[3], ", ", bb[4],"]")))
  cat(sprintf(paste0("  %-", labw, "s %s\n"), "NAs:",
              n_na))
  cat(sprintf(paste0("  %-", labw, "s %s\n"), "units:",
              units))

  invisible(grid)
}




##' Plot an admove grid
##'
##' @description
##' Plots the spatial grid of an \code{admove_grid} object or of an
##' \code{admove} object containing a grid.
##'
##' @param x A grid object of class \code{"admove_grid"} or an object containing
##'   a grid, such as \code{"admove_data"}, \code{"admove_sim"}, or
##'   \code{"admove"}.
##' @param main Main title for the plot. Default is \code{"Grid"}.
##' @param labels Logical; if \code{TRUE}, cell numbers are plotted at the cell
##'   centres.
##' @param plot_grid Logical; if \code{TRUE}, grid lines are added.
##' @param plot_land Logical; if \code{TRUE}, land masses are added to the plot.
##' @param plot_bg Logical; if \code{TRUE}, active cells are shaded in the
##'   background.
##' @param auto_layout Logical; if \code{TRUE}, graphical parameters are set
##'   automatically and restored on exit.
##' @param xlab Label for the x-axis.
##' @param ylab Label for the y-axis.
##' @param bg Optional background colour for the plotting device.
##' @param ... Additional arguments passed to [plot()].
##'
##' @return
##' Invisibly returns \code{NULL}. The function is called for its plotting side
##' effects.
##'
##' @examples
##' plot_grid(skjepo$grid)
##'
##' @name plot_grid
##' @export
plot_grid <- function(x,
                      main = "Grid",
                      labels = TRUE,
                      plot_grid = TRUE,
                      plot_land = FALSE,
                      plot_bg = TRUE,
                      auto_layout = TRUE,
                      xlab = "x",
                      ylab = "y",
                      bg = NULL,
                      ...) {

  if(inherits(x, "admove_sim") || inherits(x, "admove_data")){
    grid <- x$grid
  }else if(inherits(x, "admove")){
    grid <- x$dat$grid
  }else{
    grid <- x
  }

  bb <- bbox_grid(grid)
  xlims <- bb[1:2]
  ylims <- bb[3:4]

  if(auto_layout){
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar))
    par(mfrow = c(1,1))
  }
  if(!is.null(bg)){
    par(bg = bg)
  }
  plot(xlims, ylims,
       xlim = xlims,
       ylim = ylims,
       main = main,
       ty = "n",
       asp = 1,
       xaxs = "i", yaxs = "i",
       xlab = xlab, ylab = ylab,
       ...)
  ## if(!is.null(bg)){
  ##     usr <- par("usr")
  ##     rect(usr[1], usr[3], usr[2], usr[4], col = bg, border = NA)
  ## }
  c0 <- grid$celltable
  c0[c0 > 0] <- 1
  if (isTRUE(plot_bg)) {
    image(grid$xgr,
          grid$ygr,
          c0,
          col = adjustcolor("dodgerblue2",0.2),
          xlim = xlims,
          ylim = ylims,
          add = TRUE)
  }
  if (isTRUE(plot_land)) {
    plot_land(sref = sref(grid))
  }
  labs <- as.numeric(grid$celltable)
  labs <- labs[!is.na(labs)]
  if(labels) text(grid$xygrid[,1], grid$xygrid[,2], labs)
  if (isTRUE(plot_grid)) {
    abline(v = grid$xgr, lty = 3)
    abline(h = grid$ygr, lty = 3)
  }
  box(lwd = 1.5)
}




##' Get neighbouring cells of an admove grid
##'
##' @description
##' Returns the four-neighbour structure of an \code{admove_grid}, listing for
##' each active cell the indices of its neighbouring cells above, below, left,
##' and right.
##'
##' @param grid A grid object of class \code{"admove_grid"}, as returned by
##'   [create_grid()].
##'
##' @return
##' A matrix with one row per active cell and columns giving the cell index and
##' the indices of its neighbouring cells.
##'
##' @examples
##' grid <- create_grid()
##' neighbours <- get_neighbours(grid)
##'
##' @export
get_neighbours <- function(grid) {

  .check_class(grid, "admove_grid")

  celltable <- grid$celltable

  AT <- array(c(celltable,
                ## north/top
                cbind(celltable[,-1], NA),
                ## south/down
                cbind(NA, celltable[,-ncol(celltable)]),
                ## west/right
                rbind(NA, celltable[-nrow(celltable),]),
                ## east/left
                rbind(celltable[-1,], NA)
                ),
              dim = c(nrow(celltable), ncol(celltable), 5))
  labs <- c("c", "top", "down", "left", "right")
  nextTo <- apply(AT, 3, function(x) x)
  nextTo <- nextTo[!is.na(nextTo[,1]),]
  nextTo <- nextTo[order(nextTo[,1]),]
  colnames(nextTo) <- labs

  return(nextTo)
}


##' Get x-coordinates of grid-cell centres
##'
##' @description
##' Returns the x-coordinates of the cell centres of an \code{admove_grid}.
##'
##' @param x A grid object of class \code{"admove_grid"}, as returned by
##'   [create_grid()].
##'
##' @return
##' A numeric vector of x-coordinates for the grid-cell centres.
##'
##' @examples
##' grid <- create_grid()
##' xcens <- x_centers(grid)
##'
##' @export
x_centers <- function(x) {
    .check_class(x, "admove_grid")
  x$xgr[-1] - 0.5 * x$cellsize[1]
}

##' Get y-coordinates of grid-cell centres
##'
##' @description
##' Returns the y-coordinates of the cell centres of an \code{admove_grid}.
##'
##' @param x A grid object of class \code{"admove_grid"}, as returned by
##'   [create_grid()].
##'
##' @return
##' A numeric vector of y-coordinates for the grid-cell centres.
##'
##' @examples
##' grid <- create_grid()
##' ycens <- y_centers(grid)
##'
##' @export
y_centers <- function(x) {
  .check_class(x, "admove_grid")
  x$ygr[-1] - 0.5 * x$cellsize[2]
}



##' Add a one-cell buffer around an admove grid
##'
##' @description
##' Expands an \code{admove_grid} by adding a one-cell buffer around its outer
##' boundary. Missing cells in the original grid are propagated to the buffered
##' grid, and boundary cells adjacent to missing regions may also be excluded.
##'
##' @param grid A grid object of class \code{"admove_grid"}, as returned by
##'   [create_grid()].
##' @param plot Logical; if \code{TRUE}, the buffered grid is plotted.
##'
##' @return
##' An object of class \code{"admove_grid"} with an additional one-cell buffer
##' around the original spatial domain.
##'
##' @examples
##' grid <- create_grid()
##' grid_with_buffer <- add_buffer(grid)
##'
##' @export
add_buffer <- function(grid, plot = FALSE) {

  .check_class(grid, "admove_grid")

  gridin <- grid

  cellsize <- gridin$cellsize
  xrange <- gridin$xrange + c(-cellsize[1], cellsize[1])
  yrange <- gridin$yrange + c(-cellsize[2], cellsize[2])

  xgr <- c(grid$xgr[1] - cellsize[1],
           grid$xgr,
           grid$xgr[length(grid$xgr)] + cellsize[1])
  ygr <- c(grid$ygr[1] - cellsize[2],
           grid$ygr,
           grid$ygr[length(grid$ygr)] + cellsize[2])

  xcen <- x_centers(gridin)
  xcen <- c(xcen[1] - cellsize[1],
            xcen,
            xcen[length(xcen)] + cellsize[1])
  ycen <- y_centers(gridin)
  ycen <- c(ycen[1] - cellsize[2],
            ycen,
            ycen[length(ycen)] + cellsize[2])

  dims <- dim(grid)
  nx <- dims[1] + 2
  ny <- dims[2] + 2

  xygrid <- expand.grid(x = xcen, y = ycen)
  igrid <- expand.grid(idx = 1:length(xcen), idy = 1:length(ycen))
  celltable <- matrix(rep(NA, ((length(xgr)-1)*(length(ygr)-1))),
                      nrow = (length(xgr)-1))
  celltable[cbind(igrid$idx,igrid$idy)] <- 1:nrow(igrid)

  ## account for NAs
  mask <- is.na(gridin$celltable)
  mask.buffer <- matrix(FALSE,
                        nrow = nrow(celltable),
                        ncol = ncol(celltable))
  is.buffer <- mask.buffer
  is.buffer[c(1,nrow(is.buffer)),] <- TRUE
  is.buffer[,c(1,ncol(is.buffer))] <- TRUE
  mask.buffer[2:(nrow(mask)+1), 2:(ncol(mask)+1)] <- mask

  up <- rbind(FALSE, mask.buffer[-nrow(mask.buffer), ])
  down <- rbind(mask.buffer[-1, ], FALSE)
  left <- cbind(FALSE, mask.buffer[, -ncol(mask.buffer)])
  right <- cbind(mask.buffer[, -1], FALSE)

  mask2 <- ((mask.buffer | up | down | left | right) & is.buffer) | mask.buffer

  ## correction for corners
  mask2[1,1] <- ifelse(mask2[1,2] & mask2[2,1], TRUE, FALSE)
  mask2[1,ncol(mask2)] <- ifelse(mask2[1,ncol(mask2)-1] & mask2[2,ncol(mask2)],
                                 TRUE, FALSE)
  mask2[nrow(mask2),ncol(mask2)] <- ifelse(mask2[nrow(mask2)-1,ncol(mask2)] &
                                             mask2[nrow(mask2),ncol(mask2)-1],
                                           TRUE, FALSE)
  mask2[nrow(mask2),1] <- ifelse(mask2[nrow(mask2)-1,1] &
                                   mask2[nrow(mask2),2],
                                 TRUE, FALSE)

  celltable[mask2] <- NA

  igrid <- as.data.frame(which(!is.na(celltable), arr.ind = TRUE))
  colnames(igrid) <- c("idx","idy")

  xcen <- xgr[-1] - 0.5 * cellsize[1]
  ycen <- ygr[-1] - 0.5 * cellsize[2]
  xygrid <- data.frame(x = xcen[igrid[,1]],
                       y = ycen[igrid[,2]])
  celltable[cbind(igrid$idx,igrid$idy)] <- 1:nrow(igrid)

  grid <- list(xygrid = xygrid,
               igrid = igrid,
               celltable = celltable,
               cellsize = cellsize,
               xgr = xgr,
               ygr = ygr,
               xrange = xrange,
               yrange = yrange)
  attributes(grid) <- attributes(gridin)

  if (plot) plot(grid)

  return(grid)
}





##' Check or coerce an admove grid
##'
##' @description
##' Checks whether an object is a valid \code{admove_grid}. If the object is not
##' already of class \code{"admove_grid"}, the function attempts to convert it
##' using [create_grid()].
##'
##' @param x An object to be checked or converted.
##'
##' @details
##' If \code{x} is \code{NULL}, it is returned unchanged. If \code{x} is not an
##' \code{admove_grid}, the function tries to construct one from \code{x} using
##' [create_grid()]. If this fails, an error is thrown.
##'
##' @return
##' An object of class \code{"admove_grid"}, or \code{NULL} if \code{x} is
##' \code{NULL}.
##'
##' @export
check_grid <- function(x) {

  if (is.null(x)) {
    return(x)
  } else if(!inherits(x, "admove_grid")) {
    res <- try(create_grid(x), silent = TRUE)
    if (inherits(res, "try-error")) {
      stop("Provided grid is not of class admove_grid and couldn't convert grid to admove_grid using create_grid. Please check your grid.")
    } else {
      return(res)
    }
  } else {

    return(x)

  }
}



## s3 methods ----------------------------------------------------------------------

##' Return dimensions
##'
##' @description Return dimensions of a grid of class `admove_grid`.
##'
##' @param x Grid of class `admove_grid` as returned by the function
##'     [create_grid()].
##'
##' @return A vector with x and y dimensions (number of cells).
##'
##' @export
dim.admove_grid <- function(x) {
  c(nx = length(x$xgr) - 1,
    ny = length(x$ygr) - 1)
}

##' @method summary admove_grid
##' @rdname summarise_grid
##' @export
summary.admove_grid <- function(object, ...) {
  summarise_grid(object, ...)
}


##' @rdname print-admove
##' @method print admove_grid
##' @export
print.admove_grid <- function(x, ...) {
  tmp <- x
  attributes(tmp) <- NULL
  NextMethod("print", tmp, ...)
}


##' @rdname plot_grid
##' @export
plot.admove_grid <- function(x, ...) {
  plot_grid(x, ...)
}
