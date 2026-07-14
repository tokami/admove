

##' Plot land masses in the current plotting region
##'
##' @description
##' Add land polygons to an existing plot, using the spatial reference stored in
##' `sref`. Land is transformed to the requested coordinate reference system,
##' optionally rescaled to match the plotting units, cropped to the current plot
##' extent, and then added to the active graphics device.
##'
##' @param sref Optional spatial reference object, typically as returned by
##'   [sref()]. It should contain at least a valid CRS in `sref$crs`, and may
##'   also contain plotting units in `sref$units` and a scaling factor in
##'   `sref$crs_scale`.
##' @param col Fill colour for land polygons. Default is
##'   `grDevices::adjustcolor(grey(0.7), 0.5)`.
##' @param border Border colour for land polygons. Default is `grey(0.5)`.
##' @param download_map Logical; if `TRUE`, the land map is downloaded if needed.
##'   Otherwise, a locally available map is used when possible. Default is
##'   `FALSE`.
##' @param scale Numeric map scale passed to the internal land-data loader.
##'   Default is `110`.
##' @param verbose Logical; if `TRUE`, informative messages are printed.
##'   Default is `TRUE`.
##' @param warn_once Logical; if `TRUE`, warnings about missing CRS information
##'   are printed only once per session. Default is `TRUE`.
##'
##' @details
##' The function requires package \pkg{sf}. If `sref` does not contain a valid
##' CRS, no land is plotted. The land polygons are transformed to the requested
##' CRS, optionally multiplied by `sref$crs_scale`, and cropped to the current
##' plotting region defined by `par("usr")`.
##'
##' For geographic coordinate systems, the plotting extent is truncated to valid
##' longitude and latitude ranges before cropping.
##'
##' @return
##' Invisibly returns `NULL`. Called for its side effect of adding land masses
##' to an existing plot.
##'
##' @export
plot_land <- local({
  warned_no_crs <- FALSE

  function(sref = NULL,
           col = grDevices::adjustcolor(grey(0.7), 0.5),
           border = grey(0.5),
           download_map = FALSE,
           scale = 110,
           verbose = TRUE,
           warn_once = TRUE) {

    if (!requireNamespace("sf", quietly = TRUE)) {
      stop("Package 'sf' is required to plot land with CRS support. Please install it.")
    }

    crs <- if (!is.null(sref)) sref$crs else NULL
    units <- if (!is.null(sref)) sref$units else NULL
    crs_scale <- if (!is.null(sref)) sref$crs_scale else NULL

    if (is.null(crs) || all(is.na(crs)) || all(!nzchar(as.character(crs)))) {
      if (isTRUE(verbose) && (!isTRUE(warn_once) || !isTRUE(warned_no_crs))) {
        message("No valid crs provided in sref. Don't know how to plot land.")
        warned_no_crs <<- TRUE
      }
      return(invisible(NULL))
    } else {
      crs <- sf::st_crs(crs)
      if (is.na(crs)) stop("Invalid 'crs' provided.")
    }

    if (is.null(units)) units <- "degree"
    if (is.null(crs_scale)) crs_scale <- 1

    land_ll <- .get_land(download_map, scale = scale)
    land_crs <- sf::st_transform(land_ll, crs)
    land_fix <- sf::st_make_valid(land_crs)

    ## crs_scale maps: CRS-units -> grid-units
    ## For CRS in meters and grid in km: crs_scale = 0.001
    if (!isTRUE(all.equal(crs_scale, 1))) {
      land_fix <- sf::st_set_geometry(
        land_fix,
        sf::st_geometry(land_fix) * crs_scale
        )
    }

    usr <- par("usr")  ## c(x1, x2, y1, y2)

    is_longlat <- isTRUE(sf::st_is_longlat(sf::st_crs(land_crs)))

    if (is_longlat) {
      usr[3] <- max(-90, usr[3])
      usr[4] <- min( 90, usr[4])

      ## Split any polygon crossing the antimeridian (e.g. Fiji) at +/-180 so
      ## none spans the whole globe; otherwise such a polygon is drawn as a
      ## full-width horizontal sliver.
      land_fix <- suppressWarnings(sf::st_wrap_dateline(land_fix))
    }

    ## Crop land to a longitude window [x1, x2] (in the [-180, 180] frame for
    ## geographic coordinates) and draw it, shifting the result east by `offset`
    ## degrees. The offset lets the region east of +180 be taken from its
    ## [-180, 180] equivalent and drawn continuously across the dateline.
    draw_window <- function(x1, x2, offset) {
      if (x2 <= x1) return(invisible(NULL))
      bb <- sf::st_bbox(c(
        xmin = x1, xmax = x2,
        ymin = usr[3], ymax = usr[4]
      ), crs = sf::st_crs(land_fix))
      cr <- try(suppressWarnings(sf::st_crop(land_fix, sf::st_as_sfc(bb))),
                silent = TRUE)
      if (inherits(cr, "try-error")) {
        warning("Counldn't plot land masses. Check the spatial reference info: sref(x).")
        return(invisible(NULL))
      }
      if (nrow(cr) > 0) {
        g <- sf::st_geometry(cr)
        if (offset != 0) g <- g + c(offset, 0)
        plot(g, add = TRUE, col = col, border = border)
      }
      invisible(NULL)
    }

    if (is_longlat && usr[2] > 180) {
      ## Window crosses the antimeridian: draw the part up to +180, then the
      ## part beyond it from its negative-longitude equivalent, shifted east.
      draw_window(max(-180, usr[1]), 180, 0)
      draw_window(max(-180, usr[1] - 360), usr[2] - 360, 360)
    } else if (is_longlat) {
      draw_window(max(-180, usr[1]), min(180, usr[2]), 0)
    } else {
      draw_window(usr[1], usr[2], 0)
    }

    invisible(NULL)
  }
})



##' Plot taxis on a spatial grid
##'
##' @description
##' Plot the taxis component as arrows over the spatial prediction grid for a
##' fitted or simulated `admove` object. The function can display taxis at
##' selected time steps or the average taxis across multiple time steps.
##'
##' @param x An object of class `admove` or `admove_sim`.
##' @param select Optional index vector specifying which prediction time steps to
##'   plot. If `NULL`, all available prediction time steps are used.
##' @param select_sea Optional index vector specifying which seasons to plot. If
##'   `NULL` (default), all seasons are used.
##' @param average Logical; if `TRUE` (default), the taxis vectors are averaged
##'   over the selected time steps. If `FALSE`, taxis is plotted separately for
##'   each selected time step.
##' @param cor Optional scaling factor for arrow lengths. If `NULL`, the longest
##'   arrow is automatically scaled to one grid cell width.
##' @param col Colour of the arrows. Default is `"black"`.
##' @param alpha Transparency value. Currently not used directly in the plotting
##'   call. Default is `0.5`.
##' @param lwd Line width of the arrows. Default is `1`.
##' @param main Main title of the plot. Default is `"Taxis"`.
##' @param plot_land Logical; if `TRUE`, land masses are added using
##'   [plot_land()]. Default is `FALSE`.
##' @param image_bg Logical; if `TRUE` (default), a colour image of taxis
##'   magnitude is drawn underneath the arrows.
##' @param auto_layout Logical; if `TRUE`, the plotting layout is set
##'   automatically. If multiple time steps are plotted and `average = FALSE`,
##'   panels are arranged using [n2mfrow()]. Default is `TRUE`.
##' @param add Logical; if `TRUE`, taxis arrows are added to an existing plot.
##'   If `FALSE` (default), a new plot is created.
##' @param xlab Label for the x-axis. Default is `"x"`.
##' @param ylab Label for the y-axis. Default is `"y"`.
##' @param xaxt A character specifying the x-axis type, passed to [plot()].
##'   Default is `"s"`.
##' @param yaxt A character specifying the y-axis type, passed to [plot()].
##'   Default is `"s"`.
##' @param bg Optional background colour for the plotting device. If `NULL`
##'   (default), the current background setting is used.
##' @param ... Additional arguments passed to [plot()] when a new plot is
##'   created.
##'
##' @details
##' For objects of class `admove`, the function plots predicted taxis from
##' `x$pred$hTdx` and `x$pred$hTdy`. For objects of class `admove_sim`, taxis is
##' recomputed from the simulated covariates and parameter values.
##'
##' If `average = TRUE`, the mean taxis over the selected time steps is plotted.
##' Otherwise, one panel per selected time step is produced unless `add = TRUE`.
##'
##' @return
##' Invisibly returns `NULL`. Called for its side effect of producing a plot.
##'
##' @export
plot_taxis <- function(x,
                       select = NULL,
                       select_sea = NULL,
                       average = TRUE,
                       cor = NULL,
                       col = "black",
                       alpha = 0.5,
                       lwd = 1,
                       main = "Taxis",
                       plot_land = FALSE,
                       image_bg = TRUE,
                       auto_layout = TRUE,
                       add = FALSE,
                       xlab = "x",
                       ylab = "y",
                       xaxt = "s",
                       yaxt = "s",
                       bg = NULL,
                       ...) {

  if (inherits(x, "admove")) {
    if (is.null(select)) select <- 1:length(x$dat$pred$time)
  } else if(inherits(x, "admove_sim")) {
    if (is.null(select)) select <- 1:length(x$dat$pred$time)
  }  else stop("Don't know how to plot taxis for this object. Only implemented yet for objects of class `admove` or `admove_sim`.")

  ## detect seasonal setup (admove only)
  nsea <- 1L
  is_seasonal <- FALSE
  if (inherits(x, "admove") && !is.null(x$par$alpha)) {
    nsea <- dim(x$par$alpha)[3L]
    is_seasonal <- nsea > 1L
  }
  if (is_seasonal) {
    if (is.null(select_sea)) select_sea <- seq_len(nsea)
    nsea_plot <- length(select_sea)
  } else {
    nsea_plot <- 1L
  }

  if(auto_layout){
    opar <- par(no.readonly = TRUE)
    on.exit(suppressWarnings(graphics::par(opar)))
    n_panels <- if (is_seasonal) nsea_plot
                else if (average || length(select) == 1L) 1L
                else length(select)
    par(mfrow = if (n_panels == 1L) c(1, 1) else n2mfrow(n_panels, asp = 2))
  }

  if (inherits(x, "admove")) {

    if (is_seasonal) {
      ## one taxis field per seasonal component: find break points
      ts_len <- vapply(x$dat$time_spline, length, integer(1L))
      i_sea <- which(ts_len == nsea)
      i_sea <- if (length(i_sea) > 0L) i_sea[1L] else 1L
      ts_breaks <- x$dat$time_spline[[i_sea]]
      per <- x$dat$period
      ts_upper <- c(ts_breaks[-1L], ts_breaks[1L] + per)

      ## representative absolute times: mid-season, shifted into dat$trange
      t_mid <- (ts_breaks + ts_upper) / 2
      t_ref <- x$dat$trange[1L]
      t_sea_abs <- t_ref + (t_mid - t_ref %% per + per) %% per

      kappa <- exp(x$pl$logKappa)
      ncp <- nrow(x$dat$pred$grid$xygrid)
      tax.x <- tax.y <- matrix(NA_real_, ncp, nsea)
      for (s in seq_len(nsea)) {
        tmp <- x$pred$habi$tax$grad(x$dat$pred$grid$xygrid, t_sea_abs[s])
        tax.x[, s] <- kappa * tmp[, 1]
        tax.y[, s] <- kappa * tmp[, 2]
      }
      tax.x <- tax.x[, select_sea, drop = FALSE]
      tax.y <- tax.y[, select_sea, drop = FALSE]

      ## per-panel titles showing season time interval
      mains <- if (length(main) == nsea_plot) {
        main
      } else {
        paste0(main[1L], " (Season ", select_sea, " [",
               round(ts_breaks[select_sea], 2), ", ",
               round(ts_upper[select_sea], 2), "))")
      }

    } else {
      if (average) {
        if (length(select) > 1) {
          tax.x <- apply(x$pred$hTdx[,select], 1, mean, na.rm = TRUE)
          tax.y <- apply(x$pred$hTdy[,select], 1, mean, na.rm = TRUE)
        }else{
          tax.x <- x$pred$hTdx[,select]
          tax.y <- x$pred$hTdy[,select]
        }
      }else{
        tax.x <- x$pred$hTdx[,select]
        tax.y <- x$pred$hTdy[,select]
      }
      mains <- rep(main[1L], ncol(as.matrix(tax.x)))
    }

    if(!inherits(tax.x, "matrix")){
      tax.x <- as.matrix(tax.x)
      tax.y <- as.matrix(tax.y)
    }

  if (is.null(cor)) {
    max_mag <- max(sqrt(tax.x^2 + tax.y^2), na.rm = TRUE)
    cor <- if (is.finite(max_mag) && max_mag > 0)
      x$dat$grid$cellsize[1] / max_mag else 1
  }

    for(i in 1:ncol(tax.x)){

      if (is_seasonal && add && i > 1L) {
        ## advance to the next panel in the caller's layout so each seasonal
        ## component overlays its own panel (not all on the first panel)
        mfg <- par("mfg")
        nc <- mfg[4L]; nr <- mfg[3L]
        r <- mfg[1L]; co <- mfg[2L] + 1L
        if (co > nc) { co <- 1L; r <- r + 1L }
        if (r <= nr) par(mfg = c(r, co, nr, nc))
      }

      if(!add){
        if(!is.null(bg)){
          graphics::par(bg = bg)
        }
        plot(NA,
             xlim = x$dat$pred$grid$xrange,
             ylim = x$dat$pred$grid$yrange,
             xlab = xlab,
             ylab = ylab,
             xaxt = xaxt,
             yaxt = yaxt,
             main = mains[i],
             asp = 1,
             ...)
        if (image_bg) {
          ig <- x$dat$pred$grid$igrid
          mag <- sqrt(tax.x[, i]^2 + tax.y[, i]^2)
          z <- matrix(NA_real_, length(x$dat$pred$grid$xgr) - 1L,
                      length(x$dat$pred$grid$ygr) - 1L)
          z[cbind(ig$idx, ig$idy)] <- mag
          image(x$dat$pred$grid$xgr, x$dat$pred$grid$ygr, z,
            col = adjustcolor(rev(hcl.colors(100, "YlOrRd")), 0.4),
                add = TRUE)
        }
      }
      if(plot_land){
        plot_land(sref = sref(x$dat))
      }

      arrows(x$dat$pred$grid$xygrid[,1],
             x$dat$pred$grid$xygrid[,2],
             x$dat$pred$grid$xygrid[,1] + tax.x[,i] * cor,
             x$dat$pred$grid$xygrid[,2] + tax.y[,i] * cor,
             col = col,
             lwd = lwd,
             length = .1)

      if(!add) box(lwd = 1.5)

    }

  } else if(inherits(x, "admove_sim")) {

    grid <- x$grid
    cov <- x$cov
    par <- x$par_sim
    dat <- x$dat
    funcs <- NULL

    if(is.null(par)) stop("No parameters provided! Use par = list() to specify parameters for taxis.")

    par <- default_sim_par(par)
    cov <- .make_cov_list(cov)

    trange <- range(as.numeric(attributes(cov[[1]])$dimnames[[3]]))
    if(diff(trange) == 0) trange[2] <- trange[1] + 1

    dat <- setup_data(cov = cov,
                      grid = grid,
                      trange = trange,
                      knots_tax = dat$knots_tax,
                      knots_dif = dat$knots_dif,
                      verbose = FALSE)

    dat$pred$grid$xygrid <- x$dat$pred$grid$xygrid
    dat$pred$grid$igrid <- x$dat$pred$grid$igrid

    conf <- default_conf(dat)
    funcs <- default_sim_funcs(dat, conf, par, funcs)
    hTdx.true <- sapply(dat$pred$time,
                        function(t) apply(dat$pred$grid$xygrid, 1,
                                          function(x) exp(par$logKappa) * funcs$tax(t(x),t)[1]))
    hTdy.true <- sapply(dat$pred$time,
                        function(t) apply(dat$pred$grid$xygrid, 1,
                                          function(x) exp(par$logKappa) * funcs$tax(t(x),t)[2]))

    if(average){
      if(length(select) > 1){
        tax.x <- apply(hTdx.true[,select], 1, mean, na.rm = TRUE)
        tax.y <- apply(hTdy.true[,select], 1, mean, na.rm = TRUE)
      }else{
        tax.x <- hTdx.true[,select]
        tax.y <- hTdy.true[,select]
      }
    }else{
      tax.x <- hTdx.true[,select]
      tax.y <- hTdy.true[,select]
    }

    if(!inherits(tax.x, "matrix")){
      tax.x <- as.matrix(tax.x)
      tax.y <- as.matrix(tax.y)
    }

    if (is.null(cor)) {
      max_mag <- max(sqrt(tax.x^2 + tax.y^2), na.rm = TRUE)
      cor <- if (is.finite(max_mag) && max_mag > 0)
        grid$cellsize[1] / max_mag else 1
    }

    if(!add){
      if(!is.null(bg)){
        graphics::par(bg = bg)
      }
      plot(NA,
           xlim = grid$xrange,
           ylim = grid$yrange,
           xlab = xlab,
           ylab = ylab,
           xaxt = xaxt,
           yaxt = yaxt,
           main = main,
           asp = 1,
           ...)
      if (image_bg) {
        ig <- dat$pred$grid$igrid
        mag <- rowMeans(sqrt(tax.x^2 + tax.y^2))
        z <- matrix(NA_real_, length(dat$pred$grid$xgr) - 1L,
                    length(dat$pred$grid$ygr) - 1L)
        z[cbind(ig$idx, ig$idy)] <- mag
        image(dat$pred$grid$xgr, dat$pred$grid$ygr, z,
            col = adjustcolor(rev(hcl.colors(100, "YlOrRd")), 0.4),
              add = TRUE)
      }
    }

    if(plot_land){
      plot_land(sref = sref(x$dat))
    }

    for(i in 1:ncol(tax.x)){

      arrows(dat$pred$grid$xygrid[,1],
             dat$pred$grid$xygrid[,2],
             dat$pred$grid$xygrid[,1] + tax.x[,i] * cor,
             dat$pred$grid$xygrid[,2] + tax.y[,i] * cor,
             col = col,
             lwd = lwd,
             length = .1)

    }

    if(!add) box(lwd = 1.5)

  }
}


##' Plot advection on a spatial grid
##'
##' @description
##' Plot the advection component as arrows over the spatial prediction grid for a
##' fitted or simulated `admove` object. Advection is the directed transport of
##' the animal by the environmental flow (e.g. ocean currents), scaled by the
##' estimated entrainment coefficients `gamma`; the arrows therefore show
##' `gamma * current` in coordinate units per time step. The function can display
##' advection at selected time steps or the average across multiple time steps.
##'
##' @param x An object of class `admove` or `admove_sim`.
##' @param select Optional index vector specifying which prediction time steps to
##'   plot. If `NULL`, all available prediction time steps are used.
##' @param select_sea Optional index vector selecting which seasonal components to
##'   plot when the advection coefficients are seasonal. If `NULL`, all seasons
##'   are shown.
##' @param average Logical; if `TRUE` (default), the advection vectors are
##'   averaged over the selected time steps. If `FALSE`, advection is plotted
##'   separately for each selected time step.
##' @param cor Optional scaling factor for arrow lengths. If `NULL`, the longest
##'   arrow is automatically scaled to one grid cell width.
##' @param col Colour of the arrows. Default is `"black"`.
##' @param alpha Transparency value. Currently not used directly in the plotting
##'   call. Default is `0.5`.
##' @param lwd Line width of the arrows. Default is `1`.
##' @param main Main title of the plot. Default is `"Advection"`.
##' @param plot_land Logical; if `TRUE`, land masses are added using
##'   [plot_land()]. Default is `FALSE`.
##' @param image_bg Logical; if `TRUE` (default), a colour image of advection
##'   magnitude is drawn underneath the arrows.
##' @param auto_layout Logical; if `TRUE`, the plotting layout is set
##'   automatically. If multiple time steps are plotted and `average = FALSE`,
##'   panels are arranged using [n2mfrow()]. Default is `TRUE`.
##' @param add Logical; if `TRUE`, advection arrows are added to an existing plot.
##'   If `FALSE` (default), a new plot is created.
##' @param xlab Label for the x-axis. Default is `"x"`.
##' @param ylab Label for the y-axis. Default is `"y"`.
##' @param xaxt A character specifying the x-axis type, passed to [plot()].
##'   Default is `"s"`.
##' @param yaxt A character specifying the y-axis type, passed to [plot()].
##'   Default is `"s"`.
##' @param bg Optional background colour for the plotting device. If `NULL`
##'   (default), the current background setting is used.
##' @param ... Additional arguments passed to [plot()] when a new plot is
##'   created.
##'
##' @details
##' For objects of class `admove`, the function plots predicted advection from
##' `x$pred$hAx` and `x$pred$hAy` (which already incorporate the estimated
##' `gamma`). For objects of class `admove_sim`, advection is recomputed from the
##' simulated covariates and parameter values.
##'
##' If `average = TRUE`, the mean advection over the selected time steps is
##' plotted. Otherwise, one panel per selected time step is produced unless
##' `add = TRUE`.
##'
##' @return
##' Invisibly returns `NULL`. Called for its side effect of producing a plot.
##'
##' @seealso [plot_taxis()], [plot_diffusion()]
##'
##' @export
plot_advection <- function(x,
                           select = NULL,
                           select_sea = NULL,
                           average = TRUE,
                           cor = NULL,
                           col = "black",
                           alpha = 0.5,
                           lwd = 1,
                           main = "Advection",
                           plot_land = FALSE,
                           image_bg = TRUE,
                           auto_layout = TRUE,
                           add = FALSE,
                           xlab = "x",
                           ylab = "y",
                           xaxt = "s",
                           yaxt = "s",
                           bg = NULL,
                           ...) {

  if (inherits(x, "admove")) {
    if (is.null(select)) select <- 1:length(x$dat$pred$time)
    if (isFALSE(x$conf$use_advection))
      warning("This model was fitted without advection (conf$use_advection = FALSE); all advection arrows will be zero.")
  } else if(inherits(x, "admove_sim")) {
    if (is.null(select)) select <- 1:length(x$dat$pred$time)
  }  else stop("Don't know how to plot advection for this object. Only implemented yet for objects of class `admove` or `admove_sim`.")

  ## detect seasonal setup from the advection coefficients (admove only)
  nsea <- 1L
  is_seasonal <- FALSE
  if (inherits(x, "admove") && !is.null(x$par$gamma)) {
    nsea <- dim(x$par$gamma)[3L]
    is_seasonal <- nsea > 1L
  }
  if (is_seasonal) {
    if (is.null(select_sea)) select_sea <- seq_len(nsea)
    nsea_plot <- length(select_sea)
  } else {
    nsea_plot <- 1L
  }

  if(auto_layout){
    opar <- par(no.readonly = TRUE)
    on.exit(suppressWarnings(graphics::par(opar)))
    n_panels <- if (is_seasonal) nsea_plot
                else if (average || length(select) == 1L) 1L
                else length(select)
    par(mfrow = if (n_panels == 1L) c(1, 1) else n2mfrow(n_panels, asp = 2))
  }

  if (inherits(x, "admove")) {

    if (is_seasonal) {
      ## one advection field per seasonal component: find break points
      ts_len <- vapply(x$dat$time_spline, length, integer(1L))
      i_sea <- which(ts_len == nsea)
      i_sea <- if (length(i_sea) > 0L) i_sea[1L] else 1L
      ts_breaks <- x$dat$time_spline[[i_sea]]
      per <- x$dat$period
      ts_upper <- c(ts_breaks[-1L], ts_breaks[1L] + per)

      ## representative absolute times: mid-season, shifted into dat$trange
      t_mid <- (ts_breaks + ts_upper) / 2
      t_ref <- x$dat$trange[1L]
      t_sea_abs <- t_ref + (t_mid - t_ref %% per + per) %% per

      ncp <- nrow(x$dat$pred$grid$xygrid)
      adv.x <- adv.y <- matrix(NA_real_, ncp, nsea)
      for (s in seq_len(nsea)) {
        adv.x[, s] <- x$pred$habi$adv_x$val(x$dat$pred$grid$xygrid, t_sea_abs[s])
        adv.y[, s] <- x$pred$habi$adv_y$val(x$dat$pred$grid$xygrid, t_sea_abs[s])
      }
      adv.x <- adv.x[, select_sea, drop = FALSE]
      adv.y <- adv.y[, select_sea, drop = FALSE]

      ## per-panel titles showing season time interval
      mains <- if (length(main) == nsea_plot) {
        main
      } else {
        paste0(main[1L], " (Season ", select_sea, " [",
               round(ts_breaks[select_sea], 2), ", ",
               round(ts_upper[select_sea], 2), "))")
      }

    } else {
      if (average) {
        if (length(select) > 1) {
          adv.x <- apply(x$pred$hAx[,select], 1, mean, na.rm = TRUE)
          adv.y <- apply(x$pred$hAy[,select], 1, mean, na.rm = TRUE)
        }else{
          adv.x <- x$pred$hAx[,select]
          adv.y <- x$pred$hAy[,select]
        }
      }else{
        adv.x <- x$pred$hAx[,select]
        adv.y <- x$pred$hAy[,select]
      }
      mains <- rep(main[1L], ncol(as.matrix(adv.x)))
    }

    if(!inherits(adv.x, "matrix")){
      adv.x <- as.matrix(adv.x)
      adv.y <- as.matrix(adv.y)
    }

  if (is.null(cor)) {
    max_mag <- max(sqrt(adv.x^2 + adv.y^2), na.rm = TRUE)
    cor <- if (is.finite(max_mag) && max_mag > 0)
      x$dat$grid$cellsize[1] / max_mag else 1
  }

    for(i in 1:ncol(adv.x)){

      if (is_seasonal && add && i > 1L) {
        ## advance to the next panel in the caller's layout so each seasonal
        ## component overlays its own panel (not all on the first panel)
        mfg <- par("mfg")
        nc <- mfg[4L]; nr <- mfg[3L]
        r <- mfg[1L]; co <- mfg[2L] + 1L
        if (co > nc) { co <- 1L; r <- r + 1L }
        if (r <= nr) par(mfg = c(r, co, nr, nc))
      }

      if(!add){
        if(!is.null(bg)){
          graphics::par(bg = bg)
        }
        plot(NA,
             xlim = x$dat$pred$grid$xrange,
             ylim = x$dat$pred$grid$yrange,
             xlab = xlab,
             ylab = ylab,
             xaxt = xaxt,
             yaxt = yaxt,
             main = mains[i],
             asp = 1,
             ...)
        if (image_bg) {
          ig <- x$dat$pred$grid$igrid
          mag <- sqrt(adv.x[, i]^2 + adv.y[, i]^2)
          z <- matrix(NA_real_, length(x$dat$pred$grid$xgr) - 1L,
                      length(x$dat$pred$grid$ygr) - 1L)
          z[cbind(ig$idx, ig$idy)] <- mag
          image(x$dat$pred$grid$xgr, x$dat$pred$grid$ygr, z,
            col = adjustcolor(rev(hcl.colors(100, "YlOrRd")), 0.4),
                add = TRUE)
        }
      }
      if(plot_land){
        plot_land(sref = sref(x$dat))
      }

      arrows(x$dat$pred$grid$xygrid[,1],
             x$dat$pred$grid$xygrid[,2],
             x$dat$pred$grid$xygrid[,1] + adv.x[,i] * cor,
             x$dat$pred$grid$xygrid[,2] + adv.y[,i] * cor,
             col = col,
             lwd = lwd,
             length = .1)

      if(!add) box(lwd = 1.5)

    }

  } else if(inherits(x, "admove_sim")) {

    grid <- x$grid
    cov <- x$cov
    par <- x$par_sim
    dat <- x$dat
    funcs <- NULL

    if(is.null(par)) stop("No parameters provided! Use par = list() to specify parameters for advection.")

    par <- default_sim_par(par)
    cov <- .make_cov_list(cov)

    trange <- range(as.numeric(attributes(cov[[1]])$dimnames[[3]]))
    if(diff(trange) == 0) trange[2] <- trange[1] + 1

    dat <- setup_data(cov = cov,
                      grid = grid,
                      trange = trange,
                      knots_tax = dat$knots_tax,
                      knots_dif = dat$knots_dif,
                      verbose = FALSE)

    dat$pred$grid$xygrid <- x$dat$pred$grid$xygrid
    dat$pred$grid$igrid <- x$dat$pred$grid$igrid

    conf <- default_conf(dat)
    conf$use_advection <- TRUE
    funcs <- default_sim_funcs(dat, conf, par, funcs)
    if (is.null(funcs$adv))
      stop("No advection function available for this simulated object (no gamma / currents).")
    hAx.true <- sapply(dat$pred$time,
                       function(t) apply(dat$pred$grid$xygrid, 1,
                                         function(x) funcs$adv(t(x),t)[1]))
    hAy.true <- sapply(dat$pred$time,
                       function(t) apply(dat$pred$grid$xygrid, 1,
                                         function(x) funcs$adv(t(x),t)[2]))

    if(average){
      if(length(select) > 1){
        adv.x <- apply(hAx.true[,select], 1, mean, na.rm = TRUE)
        adv.y <- apply(hAy.true[,select], 1, mean, na.rm = TRUE)
      }else{
        adv.x <- hAx.true[,select]
        adv.y <- hAy.true[,select]
      }
    }else{
      adv.x <- hAx.true[,select]
      adv.y <- hAy.true[,select]
    }

    if(!inherits(adv.x, "matrix")){
      adv.x <- as.matrix(adv.x)
      adv.y <- as.matrix(adv.y)
    }

    if (is.null(cor)) {
      max_mag <- max(sqrt(adv.x^2 + adv.y^2), na.rm = TRUE)
      cor <- if (is.finite(max_mag) && max_mag > 0)
        grid$cellsize[1] / max_mag else 1
    }

    if(!add){
      if(!is.null(bg)){
        graphics::par(bg = bg)
      }
      plot(NA,
           xlim = grid$xrange,
           ylim = grid$yrange,
           xlab = xlab,
           ylab = ylab,
           xaxt = xaxt,
           yaxt = yaxt,
           main = main,
           asp = 1,
           ...)
      if (image_bg) {
        ig <- dat$pred$grid$igrid
        mag <- rowMeans(sqrt(adv.x^2 + adv.y^2))
        z <- matrix(NA_real_, length(dat$pred$grid$xgr) - 1L,
                    length(dat$pred$grid$ygr) - 1L)
        z[cbind(ig$idx, ig$idy)] <- mag
        image(dat$pred$grid$xgr, dat$pred$grid$ygr, z,
            col = adjustcolor(rev(hcl.colors(100, "YlOrRd")), 0.4),
              add = TRUE)
      }
    }

    if(plot_land){
      plot_land(sref = sref(x$dat))
    }

    for(i in 1:ncol(adv.x)){

      arrows(dat$pred$grid$xygrid[,1],
             dat$pred$grid$xygrid[,2],
             dat$pred$grid$xygrid[,1] + adv.x[,i] * cor,
             dat$pred$grid$xygrid[,2] + adv.y[,i] * cor,
             col = col,
             lwd = lwd,
             length = .1)

    }

    if(!add) box(lwd = 1.5)

  }
}


##' Plot diffusion on a spatial grid
##'
##' @description
##' `plot_diffusion()` plots the diffusion component over the spatial prediction
##' grid for a fitted or simulated `admove` object.
##'
##' For fitted `admove` objects, the function plots the mean predicted diffusion
##' across prediction times. For `admove_sim` objects, diffusion is reconstructed
##' from the simulated covariates and parameter values and then averaged over
##' prediction times.
##'
##' @param x An object of class `admove` or `admove_sim`.
##' @param cor Optional scaling factor controlling the size of the diffusion
##'   symbols. If `NULL`, the largest circle is automatically scaled so that its
##'   diameter equals one grid cell width.
##' @param col Colour of the plotted diffusion symbols. Default: `"black"`.
##' @param alpha Transparency value. Currently not used directly in the plotting
##'   call. Default: `0.5`.
##' @param lwd Line width used for the plotted symbols. Default: `1`.
##' @param main Main title of the plot. Default: `"Diffusion"`.
##' @param plot_land Logical; if `TRUE`, land masses are added using
##'   [plot_land()]. Default: `FALSE`.
##' @param image_bg Logical; if `TRUE` (default), a colour image of diffusion
##'   intensity is drawn underneath the circles.
##' @param auto_layout Logical; if `TRUE`, graphical parameters are set and
##'   restored automatically. Default: `TRUE`.
##' @param add Logical; if `TRUE`, diffusion is added to an existing plot. If
##'   `FALSE` (default), a new plot is created.
##' @param xlab Label for the x-axis. Default: `"x"`.
##' @param ylab Label for the y-axis. Default: `"y"`.
##' @param xaxt A character specifying the x-axis type, passed to [plot()].
##'   Default: `"s"`.
##' @param yaxt A character specifying the y-axis type, passed to [plot()].
##'   Default: `"s"`.
##' @param bg Optional background colour for the plotting device. If `NULL`
##'   (default), the current background setting is used.
##' @param ... Additional graphical arguments passed to [plot()] when a new plot
##'   is created.
##'
##' @details
##' Diffusion is represented by point sizes on the prediction grid, with larger
##' symbols indicating higher diffusion. For fitted objects, diffusion is based
##' on `x$pred$hD`. For simulated objects, diffusion is reconstructed from the
##' simulation setup using [default_sim_funcs()].
##'
##' @return
##' Invisibly returns `NULL`. Called for its side effect of producing a plot.
##'
##' @export
plot_diffusion <- function(x,
                           cor = NULL,
                           col = "black",
                           alpha = 0.5,
                           lwd = 1,
                           main = "Diffusion",
                           plot_land = FALSE,
                           image_bg = TRUE,
                           auto_layout = TRUE,
                           add = FALSE,
                           xlab = "x",
                           ylab = "y",
                           xaxt = "s",
                           yaxt = "s",
                           bg = NULL,
                           ...) {

  if(auto_layout){
    opar <- par(no.readonly = TRUE)
    on.exit(suppressWarnings(graphics::par(opar)))
    par(mfrow = c(1,1))
  }

  if (inherits(x, "admove")) {

    if (!add) {
      plot(NA,
           xlim = x$dat$pred$grid$xrange,
           ylim = x$dat$pred$grid$yrange,
           xlab = xlab,
           ylab = ylab,
           xaxt = xaxt,
           yaxt = yaxt,
           main = main,
           asp = 1,
           ...)
    }

    dif.est <- exp(apply(x$pred$hD, 1, mean))

    if (is.null(cor)) {
      max_size <- max(sqrt(dif.est), na.rm = TRUE)
      char_u <- graphics::par("cxy")[1]
      cor <- if (is.finite(max_size) && max_size > 0 && is.finite(char_u) && char_u > 0)
        (x$dat$pred$grid$cellsize[1] / char_u) / max_size else 1
    }

    if (image_bg && !add) {
      ig <- x$dat$pred$grid$igrid
      z <- matrix(NA_real_, length(x$dat$pred$grid$xgr) - 1L,
                  length(x$dat$pred$grid$ygr) - 1L)
      z[cbind(ig$idx, ig$idy)] <- dif.est
      image(x$dat$pred$grid$xgr, x$dat$pred$grid$ygr, z,
            col = adjustcolor(rev(hcl.colors(100, "YlOrRd")), 0.4),
            add = TRUE)
    }

    if (isTRUE(plot_land)) {
      plot_land(sref = sref(x$dat$grid))
    }

    points(x$dat$pred$grid$xygrid[,1],
           x$dat$pred$grid$xygrid[,2],
           col = col,
           lwd = lwd,
           cex = sqrt(dif.est) * cor)

  } else if(inherits(x, "admove_sim")) {

    grid <- x$grid
    cov <- x$cov
    par <- x$par_sim
    dat <- x$dat
    funcs <- NULL

    if(is.null(par)) stop("No parameters provided! Use par = list() to specify parameters for diffusion.")

    if (!add) {
      if(!is.null(bg)){
        graphics::par(bg = bg)
      }
      plot(NA,
           xlim = grid$xrange,
           ylim = grid$yrange,
           xlab = xlab,
           ylab = ylab,
           xaxt = xaxt,
           yaxt = yaxt,
           main = main,
           asp = 1,
           ...)
      ## if(!is.null(bg)){
      ##     usr <- par("usr")
      ##     rect(usr[1], usr[3], usr[2], usr[4], col = bg, border = NA)
      ## }
    }

    par <- default_sim_par(par)
    cov <- .make_cov_list(cov)

    trange <- range(as.numeric(attributes(cov[[1]])$dimnames[[3]]))
    if(diff(trange) == 0) trange[2] <- trange[1] + 1

    dat <- setup_data(cov = cov,
                      grid = grid,
                      trange = trange,
                      ## trange = c(0,
                      ##            max(sapply(cov,
                      ##                       function(x) dim(x)[3]))),
                      knots_tax = dat$knots_tax,
                      knots_dif = dat$knots_dif,
                      verbose = FALSE)

    dat$pred$grid$xygrid <- x$dat$pred$grid$xygrid
    dat$pred$grid$igrid <- x$dat$pred$grid$igrid

    conf <- default_conf(dat)
    funcs <- default_sim_funcs(dat, conf, par, funcs)

    D.true <- sapply(dat$pred$time,
                     function(t) apply(dat$pred$grid$xygrid, 1,
                                       function(x) exp(funcs$dif(as.matrix(x),t)[1])))
    dif_avg <- rowMeans(D.true)

    if (is.null(cor)) {
      max_size <- max(sqrt(dif_avg), na.rm = TRUE)
      char_u <- graphics::par("cxy")[1]
      cor <- if (is.finite(max_size) && max_size > 0 && is.finite(char_u) && char_u > 0)
        (dat$pred$grid$cellsize[1] / char_u) / max_size else 1
    }

    if (image_bg && !add) {
      ig <- dat$pred$grid$igrid
      z <- matrix(NA_real_, length(dat$pred$grid$xgr) - 1L,
                  length(dat$pred$grid$ygr) - 1L)
      z[cbind(ig$idx, ig$idy)] <- dif_avg
      image(dat$pred$grid$xgr, dat$pred$grid$ygr, z,
            col = adjustcolor(rev(hcl.colors(100, "YlOrRd")), 0.4),
            add = TRUE)
    }

    if (isTRUE(plot_land)) {
      plot_land(sref = sref(x$dat$grid))
    }

    points(dat$pred$grid$xygrid[,1],
           dat$pred$grid$xygrid[,2],
           col = col,
           lwd = lwd,
           cex = sqrt(dif_avg) * cor)

  }
  if(!add) box(lwd = 1.5)
}




##' Compare a single quantity across fitted or simulated admove objects
##'
##' @description
##' `plot_compare_one()` provides a compact plotting interface for comparing a
##' single quantity across multiple fitted or simulated \emph{admove} objects.
##' Depending on the selected `quantity`, the function compares habitat
##' preference functions, taxis, diffusion, or parameter estimates.
##'
##' The input can be supplied either as individual `admove` / `admove_sim`
##' objects or as a list of such objects.
##'
##' @param fit An object of class `admove` or `admove_sim`, or a list of such
##'   objects.
##' @param ... Additional `admove` or `admove_sim` objects to compare.
##' @param quantity Character string specifying which quantity to compare.
##'   Currently implemented options are:
##'   \itemize{
##'     \item `"pref"` for habitat preference functions,
##'     \item `"taxis"` for taxis,
##'     \item `"dif"` for diffusion, and
##'     \item `"par"` for parameter estimates.
##'   }
##' @param plot_land Logical; if `TRUE`, add land masses to spatial plots using
##'   [plot_land()]. Default: `FALSE`.
##' @param auto_layout Logical; if `TRUE`, restore graphical parameters on exit.
##'   Default: `TRUE`.
##' @param asp Positive numeric value giving the target aspect ratio
##'   (columns / rows) of the plot arrangement. Default: `2`.
##' @param col Vector of colours used for the different objects being compared.
##'   By default, colours are taken from `.admove_cols(10)`.
##' @param lty Vector of line types used for the different objects being
##'   compared. Default: `1:10`.
##' @param cor_tax Optional scaling factor for taxis arrows. If `NULL`
##'   (default), the longest arrow is automatically scaled to one grid cell
##'   width.
##' @param cor_dif Optional scaling factor for diffusion symbols. If `NULL`
##'   (default), the largest circle is automatically scaled to one grid cell
##'   width.
##' @param cor_adv Optional scaling factor for advection arrows. If `NULL`
##'   (default), the longest arrow is automatically scaled to one grid cell
##'   width.
##' @param plot.legend Logical or integer indicating whether, or which, legend
##'   should be plotted. Default: `1`.
##' @param bg Optional background colour for the plotting device. Default:
##'   `NULL`.
##' @param panel_lab Optional character string added as a panel label in the
##'   top-left corner of the plot. Default: `NULL`.
##' @param select Optional index vector specifying which prediction time steps to
##'   use for spatial quantities. If `NULL` (default), all available prediction
##'   time steps are used.
##'
##' @return
##' No return value. Called for its side effect of producing a comparison plot.
##'
##' @details
##' For `quantity = "pref"`, habitat preference curves are overlaid in a single
##' panel. For `quantity = "taxis"` and `quantity = "dif"`, spatial movement
##' components are compared on the prediction grid. For `quantity = "par"`,
##' fitted parameter estimates are shown, with confidence intervals for fitted
##' `admove` objects where available.
##'
##' Simulated objects of class `admove_sim` can be included in the comparison,
##' allowing direct visual comparison between simulated and fitted quantities.
##'
##' @export
plot_compare_one <- function(fit, ...,
                             quantity = c("pref","taxis","advection",
                                          "pref_dif","dif",
                                          "par"),
                             plot_land = FALSE,
                             auto_layout = TRUE,
                             asp = 2,
                             col = .admove_cols(10),
                             lty = 1:10,
                             cor_tax = NULL,
                             cor_dif = NULL,
                             cor_adv = NULL,
                             plot.legend = 1,
                             panel_lab = NULL,
                             select = NULL,
                             bg = NULL) {

  if("admove" %in% class(fit) || "admove_sim" %in% class(fit)){
    fitlist <- list(fit, ...)
  }else if(inherits(fit, "list")){
    fitlist <- c(fit, ...)
  }else stop("Please provide fitted admove objects either individually or as list.")

  sim_ind <- lapply(fitlist, function(x) inherits(x, "admove_sim"))

  quantity <- match.arg(quantity)
  n <- length(fitlist)

  if (quantity == "pref") {
    ylims <- range(sapply(fitlist, function(x)
      plot_pref_func(x, select = select, return_limits = TRUE)$ylim))

    plot_pref_func(fitlist[[1]],
                    select = select,
                    cols = if (n == 1L) col else col[1],
                    main = "",
                    auto_layout = FALSE,
                    panel_lab = panel_lab,
                    bg = bg,
                    ylim = ylims)
    if (n > 1) {
      for(i in 2:n){
        plot_pref_func(fitlist[[i]], add = TRUE,
                        select = select,
                        cols = col[i], lty = lty[i],
                        auto_layout = FALSE,
                        bg = bg)
      }
    }
  }

  if (quantity == "pref_dif") {
    ## diffusion as a function of the covariate(s) (the diffusion spline),
    ## mirroring the "pref" (taxis) panel but with type = "diffusion"
    ylims <- range(sapply(fitlist, function(x)
      plot_pref_func(x, type = "diffusion", select = select, return_limits = TRUE)$ylim))

    plot_pref_func(fitlist[[1]],
                    type = "diffusion",
                    select = select,
                    cols = if (n == 1L) col else col[1],
                    main = "",
                    auto_layout = FALSE,
                    panel_lab = panel_lab,
                    bg = bg,
                    ylim = ylims)
    if (n > 1) {
      for(i in 2:n){
        plot_pref_func(fitlist[[i]], type = "diffusion", add = TRUE,
                        select = select,
                        cols = col[i], lty = lty[i],
                        auto_layout = FALSE,
                        bg = bg)
      }
    }
  }

  if (quantity == "taxis") {
    ## detect number of seasonal panels from the first fit
    nsea_cmp <- if (!is.null(fitlist[[1L]]$par$alpha)) dim(fitlist[[1L]]$par$alpha)[3L] else 1L
    is_sea_cmp <- nsea_cmp > 1L

    ## draw one panel per season: all fits overlaid before advancing to the
    ## next panel — mirrors the pattern used by plot_pref_func so no
    ## par(mfg=...) back-navigation is needed
    for (s in seq_len(nsea_cmp)) {
      sel_s <- if (is_sea_cmp) s else NULL
      lab_s <- if (!is.null(panel_lab) && s <= length(panel_lab)) panel_lab[s] else NULL
      plot_taxis(fitlist[[1L]], col = col[1L],
                 main = "",
                 cor = cor_tax,
                 select_sea = sel_s,
                 auto_layout = FALSE,
                 plot_land = plot_land,
                 bg = bg)
      if (!is.null(lab_s)) add_lab(lab_s)
      if (n > 1L) {
        for (i in 2L:n) {
          plot_taxis(fitlist[[i]], add = TRUE,
                     col = col[i], lty = lty[i],
                     cor = cor_tax,
                     select_sea = sel_s,
                     auto_layout = FALSE,
                     plot_land = plot_land,
                     bg = bg)
        }
      }
    }
  }

  if (quantity == "advection") {
    ## one panel per season (advection coefficients gamma may be seasonal); all
    ## fits overlaid before advancing to the next panel — mirrors the "taxis" case
    nsea_cmp <- if (!is.null(fitlist[[1L]]$par$gamma)) dim(fitlist[[1L]]$par$gamma)[3L] else 1L
    is_sea_cmp <- nsea_cmp > 1L

    for (s in seq_len(nsea_cmp)) {
      sel_s <- if (is_sea_cmp) s else NULL
      lab_s <- if (!is.null(panel_lab) && s <= length(panel_lab)) panel_lab[s] else NULL
      plot_advection(fitlist[[1L]], col = col[1L],
                     main = "",
                     cor = cor_adv,
                     select_sea = sel_s,
                     auto_layout = FALSE,
                     plot_land = plot_land,
                     bg = bg)
      if (!is.null(lab_s)) add_lab(lab_s)
      if (n > 1L) {
        for (i in 2L:n) {
          plot_advection(fitlist[[i]], add = TRUE,
                         col = col[i],
                     cor = cor_adv,
                         select_sea = sel_s,
                         auto_layout = FALSE,
                         plot_land = plot_land,
                         bg = bg)
        }
      }
    }
  }

  if(quantity == "dif"){

    plot_diffusion(fitlist[[1]],
                   col = col[1], lty = lty[1],
                   cor = cor_dif,
                   main = "",
                   auto_layout = FALSE,
                   plot_land = plot_land,
                   bg = bg)
    if (n > 1) {
      for (i in 2:n) {
        plot_diffusion(fitlist[[i]], add = TRUE,
                       col = col[i], lty = lty[i],
                       cor = cor_dif,
                       auto_layout = FALSE,
                       plot_land = plot_land,
                       bg = bg)
      }
    }
    if (!is.null(panel_lab)) add_lab(panel_lab)
  }

  if(quantity == "par"){

    idx <- which(sapply(fitlist, function(x) inherits(x, "admove")))
    if(length(idx) > 0){
      pars <- unique(unlist(lapply(fitlist[idx], .get_par_names)))
    }


    tmp <- lapply(fitlist, function(x) {
      if (inherits(x, "admove_sim")) {
        nam <- names(x$par_sim)
        map <- names(x$map)[match(nam,names(x$map))]
        map <- map[!is.na(map)]
        notMapped <- unlist(lapply(x$map[map], function(x) !is.na(x) & !duplicated(x)))
        ## mapped <- unlist(x$map[map])
        ## mapped <- is.na(mapped)
        pars <- unlist(x$par_sim)
        pars <- pars[names(pars) %in% names(notMapped)[notMapped]]
        ## ind <- unlist(sapply(c("beta","logSdO"),
        ##                      function(x) grep(x, names(pars))))
        ## if (length(ind) > 0) {
        ##   pars[ind] <- exp(pars[ind])
        ## }
        lo <- pars
        hi <- pars
      } else if(inherits(x, "admove")) {
        nam <- unique(names(x$opt$par))
        map <- names(x$map)[match(nam,names(x$map))]
        map <- map[!is.na(map)]
        notMapped <- unlist(lapply(x$map[map], function(x) !is.na(x) & !duplicated(x)))
        if (is.null(x$pl)) {
          pars <- unlist(x$opt$par)
        } else {
          pars <- unlist(x$pl[nam])
        }
        pars <- pars[names(pars) %in% names(notMapped)[notMapped]]
        sds <- unlist(x$plsd[nam])
        sds <- sds[names(sds) %in% names(notMapped)[notMapped]]
        lo <- pars - 1.96 * sds
        hi <- pars + 1.96 * sds
        ## ind <- unlist(sapply(c("beta","logSdO"),
        ##                      function(x) grep(x, names(pars))))
        ## if(length(ind) > 0){
        ##   lo[ind] <- exp(pars[ind] - 1.96 * sds[ind])
        ##   hi[ind] <- exp(pars[ind] + 1.96 * sds[ind])
        ##   pars[ind] <- exp(pars[ind])
        ## }
      }
      return(c(pars,lo,hi))
    })

    r <- range(unlist(tmp), na.rm = TRUE)
    pad <- 0.1 * diff(r)
    ylim <- c(r[1] - pad, r[2] + pad)
    xlim <- c(1, max(unique(unlist(lapply(tmp, function(x)
      length(unique(names(x)))))))) + 0.5 * c(-1,1)

    i = 1
    if (inherits(fitlist[[i]], "admove_sim")){

      nam <- names(fitlist[[i]]$par_sim)

      map <- names(fitlist[[i]]$map)[match(nam,names(fitlist[[i]]$map))]
      map <- map[!is.na(map)]
      notMapped <- unlist(lapply(fitlist[[i]]$map[map], function(x) !is.na(x) & !duplicated(x)))
      ## mapped <- unlist(fitlist[[i]]$map[map])
      ## mapped <- is.na(mapped)

      pars <- unlist(fitlist[[i]]$par_sim[nam])
      pars <- pars[names(pars) %in% names(notMapped)[notMapped]]


        ## ind <- unlist(sapply(c("beta","logSdO"),
        ##                      function(x) grep(x, names(pars))))
        ## if (length(ind) > 0) {
        ##   pars[ind] <- exp(pars[ind])
        ## }

      ## ind <- which(names(pars) %in% c("beta","logSdO"))
      ## if(length(ind) > 0){
      ##   pars[ind] <- exp(pars[ind])
      ## }

    }else if(inherits(fitlist[[i]], "admove")){

      nam <- unique(names(fitlist[[i]]$opt$par))

      map <- names(fitlist[[i]]$map)[match(nam,names(fitlist[[i]]$map))]
      map <- map[!is.na(map)]
      ## mapped <- unlist(fitlist[[i]]$map[map])
      ## mapped <- is.na(mapped)
      notMapped <- unlist(lapply(fitlist[[i]]$map[map], function(x) !is.na(x) & !duplicated(x)))

      if(is.null(fitlist[[i]]$pl)){
        pars <- unlist(fitlist[[i]]$opt$par)
      }else{
        pars <- unlist(fitlist[[i]]$pl[nam])
      }

      pars <- pars[names(pars) %in% names(notMapped)[notMapped]]
      sds <- unlist(fitlist[[i]]$plsd[nam])
      sds <- sds[names(sds) %in% names(notMapped)[notMapped]]
      lo <- pars - 1.96 * sds
      hi <- pars + 1.96 * sds
      ##   ind <- unlist(sapply(c("beta","logSdO"),
      ##                        function(x) grep(x, names(pars))))
      ## if(length(ind) > 0){
      ##   lo[ind] <- exp(pars[ind] - 1.96 * sds[ind])
      ##   hi[ind] <- exp(pars[ind] + 1.96 * sds[ind])
      ##   pars[ind] <- exp(pars[ind])
      ## }

    }

    labs <- names(pars)
    names(labs) <- names(notMapped)[notMapped]
    if(length(grep("logSdO",labs)) > 0){
      names(labs)[grep("logSdO",labs)] <- "sdO"
    }

    if(!is.null(bg)){
      graphics::par(bg = bg)
    }
    plot(seq(pars), pars,
         ty = "n",
         xlim = xlim,
         xaxt = "n",
         ylim = ylim,
         xlab = "Parameter",
         ylab = "Value")
    ## if(!is.null(bg)){
    ##     usr <- par("usr")
    ##     rect(usr[1], usr[3], usr[2], usr[4], col = bg, border = NA)
    ## }
    axis(1, at = seq(pars), labels = names(labs))

    addi <- seq(-0.1, 0.1, length.out = n)

    if(inherits(fitlist[[i]], "admove") && length(lo) > 0){
      arrows(seq(pars) + addi[i], lo,
             seq(pars) + addi[i], hi,
             length = 0.1,
             angle = 90,
             code = 3,
             col = col[i])
    }
    points(seq(pars) + addi[i], pars, col = col[i])

    if(n > 1){
      for(i in 2:n){
        if(inherits(fitlist[[i]], "admove_sim")){
          pars <- unlist(fitlist[[i]]$par_sim)
        ## ind <- unlist(sapply(c("beta","logSdO"),
        ##                      function(x) grep(x, names(pars))))

        ##   if(length(ind) > 0){
        ##     pars[ind] <- exp(pars[ind])
        ##   }
          pars <- pars[match(names(labs), names(pars))]
        }else if(inherits(fitlist[[i]], "admove")){

          nam <- unique(names(fitlist[[i]]$opt$par))

          map <- names(fitlist[[i]]$map)[match(nam,names(fitlist[[i]]$map))]
          map <- map[!is.na(map)]
          mapped <- unlist(fitlist[[i]]$map[map])
          mapped <- is.na(mapped)

          pars <- unlist(fitlist[[i]]$pl[nam])
          pars <- pars[!names(pars) %in% names(mapped)[mapped]]
          sds <- unlist(fitlist[[i]]$plsd[nam])
          sds <- sds[!names(sds) %in% names(mapped)[mapped]]
          lo <- pars - 1.96 * sds
          hi <- pars + 1.96 * sds
        ## ind <- unlist(sapply(c("beta","logSdO"),
        ##                      function(x) grep(x, names(pars))))

        ##   if(length(ind) > 0){
        ##     lo[ind] <- exp(pars[ind] - 1.96 * sds[ind])
        ##     hi[ind] <- exp(pars[ind] + 1.96 * sds[ind])
        ##     pars[ind] <- exp(pars[ind])
        ##   }
          arrows(seq(pars) + addi[i], lo,
                 seq(pars) + addi[i], hi,
                 length = 0.1,
                 angle = 90,
                 code = 3,
                 col = col[i])
        }
        points(seq(pars) + addi[i], pars, col = col[i])
      }
    }
    box(lwd = 1.5)
    if (!is.null(panel_lab)) add_lab(panel_lab)
  }
}


##' Compare fitted and simulated admove objects
##'
##' @description
##' Create comparison plots for one or more fitted or simulated `admove`
##' objects. Depending on `quantity`, the function can compare habitat
##' preference, taxis, diffusion, or parameter estimates across objects.
##' Multiple requested quantities are arranged automatically in a multi-panel
##' layout.
##'
##' @param fit Either a single object of class `admove` or `admove_sim`, or a
##'   list of such objects. If a named list is supplied, names are used in the
##'   legend.
##' @param ... Additional `admove` or `admove_sim` objects to compare.
##' @param quantity Character vector specifying which quantities to compare.
##'   Implemented options are:
##'   \describe{
##'     \item{`"pref"`}{Habitat preference.}
##'     \item{`"taxis"`}{Taxis.}
##'     \item{`"dif"`}{Diffusion.}
##'     \item{`"par"`}{Parameter estimates.}
##'   }
##'   Multiple quantities can be selected.
##' @param plot_land Logical; if `TRUE`, land masses are added to spatial plots
##'   using [plot_land()]. Default is `FALSE`.
##' @param auto_layout Logical; if `TRUE`, the plot layout and graphical
##'   parameters are set automatically. Default is `TRUE`.
##' @param col Colours used for the different objects being compared. Defaults to
##'   `admove:::.admove_cols(10)`.
##' @param cor_dif Optional scaling factor for diffusion symbols. If `NULL`,
##'   the default internal scaling is used.
##' @param cor_tax Optional scaling factor for taxis arrows. If `NULL`,
##'   the default internal scaling is used.
##' @param asp Positive numeric value giving the target aspect ratio
##'   (columns / rows) for the plot arrangement. Default is `2`.
##' @param plot.legend Logical or integer controlling legend placement. If set
##'   to `1` (default), a shared legend is drawn in a separate layout panel. If
##'   set to `2`, the legend is added within the final plot panel.
##' @param bg Optional background colour for the plotting device. If `NULL`
##'   (default), the current background setting is used.
##'
##' @details
##' If `auto_layout = TRUE`, the function arranges the requested comparison plots
##' automatically. For spatial quantities, land can optionally be added via
##' [plot_land()]. Simulated objects of class `admove_sim` can be compared
##' directly with fitted objects of class `admove`.
##'
##' When `plot.legend = 1`, a shared legend is drawn below the plots. If the
##' input objects are unnamed, fitted objects are labelled sequentially and
##' simulated objects are labelled `"Sim"`.
##'
##' @return
##' Invisibly returns `NULL`. Called for its side effect of producing plots.
##'
##' @export
plot_compare <- function(fit, ...,
                         quantity = c("pref","taxis","advection",
                                      "dif","par"),
                         plot_land = FALSE,
                         auto_layout = TRUE,
                         col = .admove_cols(10),
                         cor_dif = NULL,
                         cor_tax = NULL,
                         asp = 2,
                         plot.legend = 1,
                         bg = NULL) {

  if("admove" %in% class(fit) || "admove_sim" %in% class(fit)){
    fitlist <- list(fit = fit, ...)
  }else if(inherits(fit, "list")){
    fitlist <- c(fit, ...)
  }else stop("Please provide fitted admove objects either individually or as list.")

  sim_ind <- lapply(fitlist, function(x) inherits(x, "admove_sim"))

  quantity <- match.arg(quantity, several.ok = TRUE)

  ## drop "advection" when no fit was configured with advection — otherwise it
  ## would produce an empty/misleading panel
  if ("advection" %in% quantity) {
    any_adv <- any(vapply(fitlist, function(x)
      inherits(x, c("admove", "admove_sim")) && isTRUE(x$conf$use_advection),
      logical(1L)))
    if (!any_adv) quantity <- quantity[quantity != "advection"]
  }

  if (length(quantity) == 0L)
    stop("Nothing to plot: the only requested quantity (\"advection\") is not ",
         "available because no fit was configured with advection ",
         "(conf$use_advection = FALSE).")

  nq <- length(quantity)

  if(!is.null(bg)){
    graphics::par(bg = bg)
  }

  ref_fit <- Filter(function(x) inherits(x, c("admove", "admove_sim")), fitlist)[[1L]]
  ncov <- if (!is.null(ref_fit$dat$cov)) length(ref_fit$dat$cov) else 1L
  nsea_ref <- if (!is.null(ref_fit$par$alpha)) dim(ref_fit$par$alpha)[3L] else 1L
  nsea_adv <- if (!is.null(ref_fit$par$gamma)) dim(ref_fit$par$gamma)[3L] else 1L
  panels_per_q <- vapply(quantity, function(q) {
    if (q == "pref") ncov
    else if (q == "taxis") nsea_ref
    else if (q == "advection") nsea_adv
    else 1L
  }, integer(1L))
  total_panels <- sum(panels_per_q)

  if(auto_layout){
    opar <- par(no.readonly = TRUE)
    on.exit(suppressWarnings(graphics::par(opar)))
    mfrow <- n2mfrow(total_panels, asp = asp)
    par(mar = c(4.5,4,1,1)+0.1, oma = c(1,1,1,1))
    if(as.integer(plot.legend) == 1){
      cells <- rep(0L, prod(mfrow))
      cells[seq_len(total_panels)] <- seq_len(total_panels)
      layout(rbind(matrix(cells,
                          nrow = mfrow[1],
                          ncol = mfrow[2],
                          byrow = TRUE),
                   rep(total_panels + 1L, mfrow[2])),
             heights = c(rep(1, mfrow[1]), 0.15))
    }else{
      layout(matrix(seq_len(total_panels),
                    nrow = mfrow[1],
                    ncol = mfrow[2],
                    byrow = TRUE))
    }
  }

  panel_start <- cumsum(c(0L, panels_per_q[-nq]))
  all_labs <- if (nq > 1L) LETTERS[seq_len(total_panels)] else NULL
  for(i in 1:nq){
    q_labs <- if (!is.null(all_labs)) all_labs[panel_start[i] + seq_len(panels_per_q[i])] else NULL
    plot_compare_one(fitlist,
                     quantity = quantity[i],
                     col = col,
                     plot.legend = as.integer(plot.legend) == 2 && i == nq,
                     plot_land = plot_land,
                     auto_layout = FALSE,
                     panel_lab = q_labs,
                     cor_dif = cor_dif,
                     cor_tax = cor_tax,
                     bg = bg)
  }

  if(as.integer(plot.legend) == 1){
    nfit <- sum(!unlist(sim_ind))
    if(is.null(names(fitlist))){
      leg.text <- ifelse(unlist(sim_ind), "Sim",
                         paste0("Fit ",
                                cumsum(!unlist(sim_ind))))
    }else{
      leg.text <- names(fitlist)
    }

    par(mar = c(1,5,0,0))
    plot.new()
    legend("center", legend = leg.text,
           lwd = 2,
           horiz = TRUE,
           bty = "s",
           box.lwd = 1.5,
           col = col[1:length(fitlist)],
           bg = "white")
  }
}




##' Plot predicted location distributions for a single tag
##'
##' @description
##' Displays the model's predicted spatial location distribution alongside the
##' observed track for a single archival tag. The predicted distributions must
##' be precomputed by [add_tag_dist()] before calling this function; an
##' informative error is raised otherwise.
##'
##' The plot is arranged as a grid: one row per tag, one column per time step.
##' Time steps are evenly spaced across each tag's observation sequence. The
##' release observation (first time step) is shown without a density background
##' because the location is known exactly at release.
##'
##' @param x A fitted object of class `admove` with `$tag_dist` added by
##'   [add_tag_dist()].
##' @param select Integer vector of tag indices to display. `NULL` (default)
##'   shows all computed tags (subject to `n_tags`).
##' @param n_tags Maximum number of tags to display (rows). `NULL` shows all
##'   selected tags.
##' @param n_time_steps Number of time steps (columns) per tag. Default is `6`.
##' @param plot_land Logical; if `TRUE`, land masses are added. Default is
##'   `FALSE`.
##' @param plot_contour Logical; if `TRUE`, contour lines are added on top of
##'   the predicted density image. Default is `FALSE`.
##' @param xlab Label for the x-axis. Default is `"x"`.
##' @param ylab Label for the y-axis. Default is `"y"`.
##' @param asp Target aspect ratio passed to [grDevices::n2mfrow()] when a
##'   single tag is displayed, controlling the shape of the time-step panel
##'   grid (larger values favour more columns). Ignored when several tags are
##'   shown (each tag keeps its own row). Default is `1`.
##'
##' @return
##' Invisibly returns `NULL`. Called for its side effect of producing plots.
##'
##' @seealso [add_tag_dist()]
##'
##' @importFrom mvtnorm dmvnorm
##'
##' @export
plot_tag_dist <- function(x,
                          select = NULL,
                          n_tags = NULL,
                          n_time_steps = 6L,
                          plot_land = FALSE,
                          plot_contour = FALSE,
                          xlab = "x",
                          ylab = "y",
                          asp = 1) {

  if (is.null(x$tag_dist))
    stop("No precomputed tag distributions found in this object.\n",
         "  Run add_tag_dist() first, e.g.:\n",
         "    fit <- add_tag_dist(fit, i = 1)")

  td_store <- x$tag_dist

  ## handle old single-entry format (list with $engine at top level)
  if (!is.null(td_store$engine)) {
    td_store <- setNames(list(td_store), as.character(td_store$i))
  }

  available <- names(td_store)

  ## resolve which tags to show
  if (is.null(select)) {
    sel_keys <- available
  } else {
    sel_keys <- as.character(select)
    missing_keys <- sel_keys[!sel_keys %in% available]
    if (length(missing_keys) > 0L)
      stop("Tag(s) ", paste(missing_keys, collapse = ", "),
           " have no precomputed distribution. Available: ",
           paste(available, collapse = ", "), ".")
  }
  if (!is.null(n_tags)) sel_keys <- head(sel_keys, n_tags)

  n_row <- length(sel_keys)
  ## columns are capped at the most observations any selected tag has, so
  ## n_time_steps acts as an upper bound rather than a fixed width (avoids
  ## trailing blank panels when tags have fewer time steps than n_time_steps)
  max_nobs <- max(vapply(sel_keys,
                         function(k) nrow(td_store[[k]]$tag), integer(1L)))
  n_col <- min(as.integer(n_time_steps), max_nobs)

  ## draw a single time-step panel for tag `td` at observation index `j`,
  ## honouring the requested axis styles (shared by both layouts below)
  draw_panel <- function(td, ind.tag, j, xaxt, yaxt) {
    tag <- td$tag
    engine <- td$engine

    xrange <- td$xrange + c(-0.1, 0.1) * diff(td$xrange)
    yrange <- td$yrange + c(-0.1, 0.1) * diff(td$yrange)

    is_release <- (j == 1L)

    plot(NA, NA,
         xlim = xrange, ylim = yrange,
         xaxt = xaxt, yaxt = yaxt,
         asp = 1, xlab = "", ylab = "")

    if (plot_land) plot_land(sref = sref(x$dat))

    ## no density at release: location is known exactly
    if (!is_release) {
      if (engine == 2L) {

        image(td$xg, td$yg, td$dens_list[[j]],
              xlim = xrange, ylim = yrange,
              col = adjustcolor(terrain.colors(100), 0.4),
              asp = 1, add = TRUE)

      } else {

        tagi <- td$traj[td$ind.track[j], ]
        mu <- as.numeric(tagi[1:2])
        Sigma <- matrix(c(tagi[3], 0, 0, tagi[4]), 2L, 2L)

        xg <- seq(xrange[1L], xrange[2L], length.out = 150L)
        yg <- seq(yrange[1L], yrange[2L], length.out = 150L)
        Z <- matrix(mvtnorm::dmvnorm(as.matrix(expand.grid(xg, yg)), mu, Sigma),
                    length(xg), length(yg))

        image(xg, yg, Z,
              xlim = xrange, ylim = yrange,
              col = adjustcolor(terrain.colors(100), 0.4),
              asp = 1, add = TRUE)

        if (plot_contour && all(!is.na(Z)))
          contour(xg, yg, Z, nlevels = 4, add = TRUE)
      }
    }

    points(tag$x[ind.tag], tag$y[ind.tag],
           type = "b", col = adjustcolor("grey20", 0.2))
    points(tag$x[j], tag$y[j],
           col = "dodgerblue3", pch = 16, cex = 1.2)

    legend("topright", legend = round(tag[j, 1L], 3),
           title.font = 2, cex = 0.8, pch = NA, x.intersp = -0.5,
           bg = "white")

    box(lwd = 1.5)
  }

  opar <- par(no.readonly = TRUE)
  on.exit(suppressWarnings(graphics::par(opar)))

  if (n_row == 1L) {
    ## single tag: spread the time steps over a compact grid rather than a
    ## single long row of panels
    td <- td_store[[sel_keys[1L]]]
    nobs <- nrow(td$tag)
    ind.tag <- if (n_col >= nobs) {
      seq_len(nobs)
    } else {
      round(seq(1L, nobs, length.out = n_col))
    }
    n_panel <- length(ind.tag)

    mfrow <- n2mfrow(n_panel, asp = asp)
    nr <- mfrow[1L]; nc <- mfrow[2L]
    par(mfrow = mfrow, mar = c(0.1, 0.1, 0.1, 0.1), oma = c(4, 4, 1, 1))

    for (p in seq_len(nr * nc)) {
      if (p > n_panel) {
        plot.new()
        next
      }
      ## mfrow fills by rows; x-axis on the bottom-most panel of each column,
      ## y-axis on the first column
      col_p <- ((p - 1L) %% nc) + 1L
      xaxt <- if (p + nc > n_panel) "s" else "n"
      yaxt <- if (col_p == 1L) "s" else "n"
      draw_panel(td, ind.tag, ind.tag[p], xaxt, yaxt)
    }
  } else {
    ## multiple tags: one row per tag, one column per time step
    par(mfrow = c(n_row, n_col), mar = c(0.1, 0.1, 0.1, 0.1), oma = c(4, 4, 1, 1))

    for (r in seq_len(n_row)) {

      td <- td_store[[sel_keys[r]]]
      nobs <- nrow(td$tag)
      ind.tag <- if (n_col >= nobs) {
        seq_len(nobs)
      } else {
        round(seq(1L, nobs, length.out = n_col))
      }

      for (c in seq_len(n_col)) {

        xaxt <- if (r == n_row) "s" else "n"
        yaxt <- if (c == 1L)   "s" else "n"

        if (c > length(ind.tag)) {
          ## pad with a blank panel to keep the grid regular
          plot.new()
          next
        }

        draw_panel(td, ind.tag, ind.tag[c], xaxt, yaxt)
      }
    }
  }

  mtext(xlab, 1, 2, outer = TRUE)
  mtext(ylab, 2, 2, outer = TRUE)

  invisible(NULL)
}


##' Add a label to a plot
##'
##' @param lab label to be added
##'
##' @export
add_lab <- function(lab){
  legend("topleft", legend = lab,
         bg = "white", x.intersp = -0.4,
         cex = 1.8, text.font = 2)
}





##' Plot habitat preference functions
##'
##' @description
##' Plot estimated or simulated habitat preference functions against covariate
##' values for taxis or diffusion. The function can display one or several
##' covariate-specific preference functions, optionally with confidence bands for
##' fitted `admove` objects.
##'
##' @param x An object of class `admove` or `admove_sim`.
##' @param type Character string specifying which preference function to plot:
##'   `"taxis"` (default) or `"diffusion"`.
##' @param select Optional index vector specifying which covariates to plot. If
##'   `NULL`, all available covariates for the selected `type` are shown.
##' @param main Optional main title. Can be a single character string or a
##'   character vector with one title per panel. If `NULL`, covariate names are
##'   used where available.
##' @param cols Colours used for the plotted preference functions. Defaults to
##'   `admove:::.admove_cols(10)`.
##' @param lwd Line width. Default is `1`.
##' @param ci Confidence level for pointwise confidence intervals. Default is
##'   `0.95`.
##' @param auto_layout Logical; if `TRUE`, the plotting layout is set
##'   automatically. Default is `TRUE`.
##' @param add Logical; if `TRUE`, the preference function is added to an
##'   existing plot. If `FALSE` (default), a new plot is created.
##' @param xlab Label for the x-axis. If `NULL` (default), the covariate name
##'   from `dat$cov` is used for each panel; falls back to `"Covariate"` when no
##'   name is available.
##' @param ylab Label for the y-axis. Default is `"Preference"`.
##' @param bg Optional background colour for the plotting device. If `NULL`
##'   (default), the current background setting is used.
##' @param ylim Optional y-axis limits.
##' @param xlim Optional x-axis limits.
##' @param return_limits Logical; if `TRUE`, no plot is produced and a list with
##'   `xlim` and `ylim` is returned instead. Default is `FALSE`.
##' @param data.range Logical; if `TRUE`, x-axis limits are based on the range of
##'   the observed covariate data rather than the prediction grid. Default is
##'   `FALSE`.
##' @param asp Positive numeric value giving the target aspect ratio
##'   (columns / rows) for multi-panel plot arrangements. Default is `2`.
##' @param leg_ncol Number of columns in the legend when seasonal curves are
##'   shown. Default is `1`.
##' @param panel_lab Optional character string added as a panel label in the
##'   top-left corner of the plot. Default is `NULL`.
##' @param ... Additional arguments passed to [plot()] when a new plot is
##'   created.
##'
##' @details
##' For fitted objects of class `admove`, the function plots estimated
##' preference functions based on predicted values stored in the fitted object.
##' If standard deviations are available, pointwise confidence bands are drawn.
##'
##' For simulated objects of class `admove_sim`, the function reconstructs the
##' corresponding preference function directly from the simulated spline knots
##' and coefficients.
##'
##' If multiple seasonal curves are available for a covariate, they are shown as
##' separate line types, and a legend is added.
##'
##' @return
##' Invisibly returns `NULL` when plotting. If `return_limits = TRUE`, returns a
##' list with components `xlim` and `ylim`.
##'
##' @export
plot_pref_func <- function(x,
                           type = "taxis",
                           select = NULL,
                           main = NULL,
                           cols = .admove_cols(10),
                           lwd = 1,
                           ci = 0.95,
                           auto_layout = TRUE,
                           add = FALSE,
                           panel_lab = NULL,
                           xlab = NULL,
                           ylab = "Preference",
                           bg = NULL,
                           ylim = NULL,
                           xlim = NULL,
                           return_limits = FALSE,
                           data.range = FALSE,
                           asp = 2,
                           leg_ncol = 1,
                           ...) {

  main0 <- main
  ylim0 <- ylim

  if (inherits(x, "admove") || inherits(x, "admove_sim")) {
    if(auto_layout && !return_limits){
      opar <- par(no.readonly = TRUE)
      on.exit(suppressWarnings(graphics::par(opar)))
      par(mfrow = c(1,1))
    }
  }

  if (inherits(x, "admove")) {

    sdr <- x$sdr
    cov_pred <- x$dat$pred$cov

    if (type == "taxis") {

      if (is.null(select)) {
        select <- 1:dim(x$par$alpha)[2]
      }

      if (!is.null(sdr)) {
        ind <- which(names(sdr$value) == "pref_taxis_pred")
        par_est <- x$pl$alpha[,select,, drop = FALSE]
      }else{
        ind <- which(names(x$rep) == "pref_taxis_pred")
        tmp <- array(x$opt$par[names(x$opt$par) == "alpha"],
                     dim = c(nrow(x$par$alpha)-1,
                             dim(x$par$alpha)[2:3]))
        par_est <- abind::abind(array(0,c(1,dim(tmp)[2:3])), tmp, along = 1)[,select, , drop = FALSE]
      }
      knots <- x$dat$knots_tax[,select]

    } else if(type == "diffusion") {

      if(is.null(select)){
        select <- 1:ncol(x$par$beta)
      }

      if(!is.null(sdr)){
        ind <- which(names(sdr$value) == "pref_dif_pred")
        par_est <- x$pl$beta[,select,, drop = FALSE]
      }else{
        ind <- which(names(x$rep) == "pref_dif_pred")
        tmp <- array(x$opt$par[names(x$opt$par) == "beta"],
                     dim = c(nrow(x$par$beta)-1,
                             dim(x$par$beta)[2:3]))
        par_est <- abind::abind(array(0,c(1,dim(tmp)[2:3])), tmp, along = 1)[,select, , drop = FALSE]
      }
      knots <- x$dat$knots_dif[,select]

    } else stop("only taxis and diffusion implemented yet.")

    if (!is.null(sdr)) {
      pref <- sdr$value[ind]
      prefsd <- sdr$sd[ind]
      if (all(is.na(prefsd) | is.nan(prefsd)))
        warning("Standard deviations are not available (NA/NaN); confidence intervals will not be shown.")
      preflow <- pref - qnorm(ci + (1 - ci)/2) * prefsd
      prefup <- pref + qnorm(ci + (1 - ci)/2) * prefsd
    } else {
      fit_rep <- if (!is.null(x$rep)) x$rep else x$obj$report()
      if (type == "taxis") {
        pref <- fit_rep[["pref_taxis_pred"]]
      } else {
        pref <- fit_rep[["pref_dif_pred"]]
      }
      prefsd <- preflow <- prefup <- rep(NA, length(pref))
    }

    pref <- array(pref, dim = c(nrow(cov_pred), ncol(cov_pred), dim(par_est)[3]))
    preflow <- array(preflow, dim = c(nrow(cov_pred), ncol(cov_pred), dim(par_est)[3]))
    prefup <- array(prefup, dim = c(nrow(cov_pred), ncol(cov_pred), dim(par_est)[3]))

    ## restrict the covariate dimension to the requested covariates (par_est and
    ## knots were already subset above), so the per-panel loop indexes 1:length(select)
    ## consistently across all arrays regardless of which covariates are selected
    pref <- pref[, select, , drop = FALSE]
    preflow <- preflow[, select, , drop = FALSE]
    prefup <- prefup[, select, , drop = FALSE]
    cov_pred <- cov_pred[, select, drop = FALSE]

    if(is.null(xlim)) xlim <- apply(cov_pred, 2, range)

    if(data.range){
      xlim <- sapply(get_cov(x$dat, x$conf)[select], range, na.rm = TRUE)
    }

    if(is.null(ylim)) ylim <- apply(rbind(apply(pref, 2, range),
                                          apply(preflow, 2, range),
                                          apply(prefup, 2, range)),2,range,
                                    na.rm = TRUE)
    alpha <- 0.3
    if(is.null(cols)) cols <- .admove_cols(length(select))
    cols <- rep_len(cols, length(select))

    if(return_limits) return(list(xlim = xlim, ylim = ylim))

    if(auto_layout && !return_limits){
      par(mfrow = n2mfrow(length(select), asp))
    }

    if (!inherits(ylim, "matrix")) {
      ylim <- matrix(rep(as.numeric(ylim), length(select)), nrow = 2L)
    }


    cov_nms <- names(x$dat$cov)

    for(i in 1:length(select)){

      if (is.null(main0)) {
        main <- ""
      } else if(length(main0) > 1) {
        main <- main0[i]
      }

      xlab_i <- if (is.null(xlab)) {
        nm <- cov_nms[select[i]]
        if (!is.null(nm) && nzchar(nm)) nm else "Covariate"
      } else xlab

      if (!add) {

        if (!is.null(bg)) {
          graphics::par(bg = bg)
        }
        plot(NA, ty = 'n',
             xlim = xlim[,i],
             ylim = ylim[,i],
             xlab = xlab_i,
             ylab = ylab,
             main = main,
             ...)
        if (!is.null(panel_lab) && length(panel_lab) >= i) add_lab(panel_lab[i])
      }

      if (!is.null(sdr)) {
        for (j in 1:dim(par_est)[3]) {
          if (dim(par_est)[3] == 1) {
            polygon(c(cov_pred[,i], rev(cov_pred[,i])),
                    c(preflow[,i,1], rev(prefup[,i,1])),
                    border = NA,
                    col = rgb(t(col2rgb(cols[i]))/255, alpha=alpha))
          } else {
            polygon(c(cov_pred[,i], rev(cov_pred[,i])),
                    c(preflow[,i,j], rev(prefup[,i,j])),
                    border = NA,
                    col = rgb(t(col2rgb(cols[i]))/255, alpha=alpha))
          }

        }
        ## rug(x$dat$cov$cov_obs[,inp$cov$var[i]])
      }

      for (j in 1:dim(par_est)[3]) {
        if (length(select) > 1) {
          knoti <- knots[,i]
          esti <- par_est[,i,j]
        } else {
          knoti <- knots
          esti <- par_est[,,j]
        }
        points(knoti, esti,
               pch = 15 + j, cex = 1.2) ##, col = cols[i])
      }
      for (j in 1:dim(par_est)[3]) {
        lines(cov_pred[,i], pref[,i,j], col = cols[i], lwd = lwd,
              lty = j)
      }

      if (dim(par_est)[3] > 1) {
        tmpi <- c(x$dat$time_spline[[i]], x$dat$period)
        leg_text <- sapply(1:dim(par_est)[3],
                           function(x) paste0("[",tmpi[x],",",
                                              ifelse(x < dim(par_est)[3],
                                                     tmpi[x+1]-1, tmpi[x+1]),"]"))
        legend("topright", legend = leg_text,
               col = cols[i],
               lty = 1:dim(par_est)[3],
               lwd = lwd,
               ncol = leg_ncol)

      }

      if(!add) box(lwd = 1.5)

    }

  } else if(inherits(x, "admove_sim")) {

    i = 1

    grid <- x$grid
    cov <- x$cov
    par <- x$par_sim
    dat <- x$dat
    funcs <- NULL

    if(is.null(par)) stop("No parameters provided! Use par = list() to specify parameters for taxis.")

    par <- default_sim_par(par)
    cov <- .make_cov_list(cov)

    trange <- range(as.numeric(attributes(cov[[1]])$dimnames[[3]]))
    if(diff(trange) == 0) trange[2] <- trange[1] + 1

    dat <- setup_data(cov = cov,
                       grid = grid,
                             trange = trange,
                             ## trange = c(0,
                             ##            max(sapply(cov,
                             ##                       function(x) dim(x)[3]))),
                             knots_tax = dat$knots_tax,
                       knots_dif = dat$knots_dif,
                       verbose = FALSE)
    conf <- default_conf(dat)
    funcs <- default_sim_funcs(dat, conf, par, funcs)

    cov_pred <- dat$pred$cov

    if (is.null(xlim)) xlim <- apply(dat$pred$cov, 2, range)

    if (type == "taxis") {
      knots <- dat$knots_tax[,i]
      par <- par$alpha[,i,]
    } else if(type == "diffusion") {
      knots <- dat$knots_dif[,i]
      par <- par$beta[,i]
    }

    get_true.pref <- .poly_fun(as.numeric(knots),
                                       as.numeric(par))

    pref <- get_true.pref(dat$pred$cov[,i])

    if(is.null(ylim)) ylim <- range(pref, par)

    if(return_limits) return(list(xlim = xlim, ylim = ylim))

    alpha <- 0.3


    sim_cov_nm <- names(dat$cov)[i]
    if (is.null(main0)) main <- ""
    xlab_i <- if (is.null(xlab)) {
      if (!is.null(sim_cov_nm) && nzchar(sim_cov_nm)) sim_cov_nm else "Covariate"
    } else xlab

    if (!add) {
      if(!is.null(bg)){
        graphics::par(bg = bg)
      }
      plot(NA, ty = 'n',
           xlim = xlim,
           ylim = ylim,
           xlab = xlab_i,
           ylab = ylab,
           main = main,
           ...)
    }
    lines(cov_pred[,i], pref,
          col = cols[i],
          lwd = lwd)
    points(knots, par,
           ## col = col,
           pch = 16, cex = 1.5)
    if(!add) box(lwd = 1.5)
  } else message("This function is only implemented for objects of class 'admove' or 'admove_sim'. Did you provide the correct object? Consider running 'sim_admove()' or 'admove()'.")
}



##' Plot spatial habitat preference surfaces
##'
##' @description
##' Plot habitat preference in space for taxis or diffusion as a raster-like
##' surface over the model grid. Preference surfaces can be shown separately for
##' selected covariates and seasons, or combined across covariates and/or
##' seasons.
##'
##' @param x An object of class `admove` or `admove_sim`.
##' @param type Character string specifying which preference surface to plot:
##'   `"taxis"` (default) or `"diffusion"`.
##' @param select_cov Optional index vector specifying which covariates to plot.
##'   If `NULL`, all available covariates for the selected `type` are used.
##' @param select.y Optional time slice to use from the covariate field. If
##'   `NULL`, the first available layer is used.
##' @param select.sea Optional index vector specifying which seasonal components
##'   to plot. If `NULL`, all available seasons are used.
##' @param combine_cov Logical; if `TRUE`, preference surfaces are summed across
##'   selected covariates before plotting. Default is `FALSE`.
##' @param combine.sea Logical; if `TRUE`, preference surfaces are summed across
##'   selected seasonal components before plotting. Default is `FALSE`.
##' @param main Optional main title for the plot panels.
##' @param col Colours for the preference surface. Defaults to
##'   `hcl.colors(14, "YlOrRd", rev = TRUE)`.
##' @param ci Confidence level for pointwise confidence intervals. Default is
##'   `0.95`.
##' @param plot_land Logical; if `TRUE`, land masses are added to the plot.
##'   Default is `FALSE`.
##' @param auto_layout Logical; if `TRUE`, the plotting layout is set
##'   automatically. Default is `TRUE`.
##' @param add Logical; if `TRUE`, the preference surface is added to an
##'   existing plot. If `FALSE` (default), a new plot is created.
##' @param xlab Label for the x-axis. Default is `"x"`.
##' @param ylab Label for the y-axis. Default is `"y"`.
##' @param bg Optional background colour for the plotting device. If `NULL`
##'   (default), the current background setting is used.
##' @param asp Positive numeric value giving the target aspect ratio
##'   (columns / rows) for multi-panel plot arrangements. Default is `2`.
##' @param ... Additional arguments passed to [plot()] when a new plot is
##'   created.
##'
##' @details
##' For fitted objects of class `admove`, the function evaluates the estimated
##' preference functions on the spatial covariate fields and displays the
##' resulting preference surface for each selected covariate and season.
##'
##' Preference surfaces can be plotted separately or combined across covariates
##' (`combine_cov = TRUE`) and/or seasons (`combine.sea = TRUE`).
##'
##' For simulated objects of class `admove_sim`, the preference surface is
##' reconstructed from the simulated covariates and parameter values.
##'
##' @return
##' Invisibly returns `NULL`. Called for its side effect of producing plots.
##'
##' @export
plot_pref_grid <- function(x,
                           type = "taxis",
                           select_cov = NULL,
                           select.y = NULL,
                           select.sea = NULL,
                           combine_cov = FALSE,
                           combine.sea = FALSE,
                           main = NULL,
                           col = hcl.colors(14, "YlOrRd", rev = TRUE),
                           ci = 0.95,
                           plot_land = FALSE,
                           auto_layout = TRUE,
                           add = FALSE,
                           xlab = "x",
                           ylab = "y",
                           bg = NULL,
                           asp = 2,
                           ...) {

  main0 <- main
  select.y0 <- select.y

  if(auto_layout){
    opar <- par(no.readonly = TRUE)
    on.exit(suppressWarnings(graphics::par(opar)))
    par(mfrow = c(1,1))
  }

  if (inherits(x, "admove")) {

    sdr <- x$sdr
    cov_pred <- x$dat$pred$cov

    if (type == "taxis") {

      if(is.null(select_cov)){
        select_cov <- 1:dim(x$par$alpha)[2]
      }

      if(is.null(select.sea)){
        select.sea <- 1:dim(x$par$alpha)[3]
      }

      if (!is.null(sdr)) {
        ind <- which(names(sdr$value) == "pref_taxis_pred")
        par_est <- x$pl$alpha[,select_cov,select.sea, drop = FALSE]
      } else {
        ind <- which(names(x$rep) == "pref_taxis_pred")
        tmp <- array(x$opt$par[names(x$opt$par) == "alpha"],
                     dim = c(nrow(x$par$alpha)-1,
                             dim(x$par$alpha)[2:3]))
        par_est <- abind::abind(array(0,c(1,dim(tmp)[2:3])), tmp, along = 1)[,select_cov, select.sea , drop = FALSE]
      }
      knots <- x$dat$knots_tax[,select_cov]
      ## if (!is.null(par)) par_true <- par$alpha[,select_cov,select.sea]


    } else if(type == "diffusion") {

      if(is.null(select_cov)){
        select_cov <- 1:dim(x$par$beta)[2]
      }
      if(is.null(select.sea)){
        select.sea <- 1:dim(x$par$beta)[3]
      }

      if(!is.null(sdr)){
        ind <- which(names(sdr$value) == "pref_dif_pred")
        par_est <- x$pl$beta[,select_cov]
      }else{
        ind <- which(names(x$rep) == "pref_dif_pred")
        par_est <- x$opt$par[names(x$opt$par) == "beta"]
        tmp <- matrix(x$opt$par[names(x$opt$par) == "beta"],
                      nrow = nrow(x$par$beta),
                      ncol = ncol(x$par$beta))
        par_est <- tmp[,select_cov]
      }
      knots <- x$dat$knots_dif[,select_cov]
      ## if(!is.null(par)) par_true <- par$beta[,select_cov]

    } else stop("only taxis and diffusion implemented yet.")

    if (!is.null(sdr)) {
      pref <- sdr$value[ind]
      prefsd <- sdr$sd[ind]
      preflow <- pref - qnorm(ci + (1 - ci)/2) * prefsd
      prefup <- pref + qnorm(ci + (1 - ci)/2) * prefsd
    } else {
      fit_rep <- if (!is.null(x$rep)) x$rep else x$obj$report()
      if(type == "taxis"){
        pref <- fit_rep[["pref_taxis_pred"]]
      }else if(type == "diffusion"){
        pref <- fit_rep[["pref_dif_pred"]]
      }
      prefsd <- preflow <- prefup <- pref
    }

    pref <- array(pref, dim = c(nrow(cov_pred), ncol(cov_pred), 4))
    preflow <- array(preflow, dim = c(nrow(cov_pred), ncol(cov_pred), 4))
    prefup <- array(prefup, dim = c(nrow(cov_pred), ncol(cov_pred), 4))

    nsea <- length(select.sea)
    ncov <- length(select_cov)
    if((nsea == 1 || combine.sea) && (ncov == 1 || combine_cov)) {
      mfrow <- c(1,1)
    } else if((nsea == 1 || combine.sea) && (ncov != 1 || !combine_cov)) {
      mfrow <- n2mfrow(ncov, asp)
    } else if(ncov == 1 || combine_cov) {
      mfrow <- n2mfrow(nsea, asp)
    } else {
      mfrow <- c(nsea, ncov)
    }
    if(auto_layout){
      par(mfrow = mfrow,
          mar = c(4.5,4,1,1)+0.1, oma = c(1,1,1,1))
    }

    ## plotting data
    mat_list <- vector("list", nsea)
    for(j in 1:nsea){
      mat_list[[j]] <- vector("list", ncov)
      for(i in 1:ncov){

        years <- as.numeric(attributes(x$dat$cov[[i]])$dimnames[[3]])
        if(!is.null(select.y0)){
          indi <- which.min(abs(years - as.numeric(select.y)))
        }else{
          indi <- 1
        }

        ## if(is.null(main0)) main <- paste0(names(x$dat$cov)[i], " (",
        ##                                   years[indi],")")

        xcov <- as.numeric(rownames(x$dat$cov[[i]][,,indi]))
        ycov <- as.numeric(colnames(x$dat$cov[[i]][,,indi]))
        xycov <- expand.grid(xcov, ycov)
        xgrid <- x$dat$grid$xgr
        ygrid <- x$dat$grid$ygr

        indix <- as.integer(cut(xycov[,1], xgrid, include.lowest = TRUE))
        indiy <- as.integer(cut(xycov[,2], ygrid, include.lowest = TRUE))

        covi <- x$dat$cov[[i]][,,indi]

        isna <- x$dat$grid$celltable[cbind(indix,indiy)]
        covi[is.na(isna)] <- NA

        if (inherits(knots, "matrix")) {
          get_true.pref <- .poly_fun(as.numeric(knots[,i]),
                                             as.numeric(par_est[,i,j]))
        } else {
          get_true.pref <- .poly_fun(knots, par_est)
        }


        pref_pred <- get_true.pref(as.numeric(covi))

        mat <- x$dat$cov[[i]][,,1]
        mat[] <- pref_pred

        mat_list[[j]][[i]] <- mat
      }
    }

    if (combine.sea) {
      nsea <- 1
      res_list <- vector("list", nsea)
      res_list[[1]] <- vector("list", ncov)
      for (i in 1:ncov) {
        res_list[[i]] <- do.call("+", lapply(mat_list, "[[", i))
      }
      mat_list <- res_list
    }

    if (combine_cov) {
      ncov <- 1
      for (i in 1:length(mat_list)) mat_list[[i]][[1]] <- do.call("+", mat_list[[i]])
    }

    ## plot
    for(j in 1:nsea){
      for(i in 1:ncov){

        if(!add){
          if(!is.null(bg)){
            graphics::par(bg = bg)
          }
          plot(NA,
               xlim = x$dat$grid$xrange,
               ylim = x$dat$grid$yrange,
               xlab = xlab,
               ylab = ylab,
               main = main,
               ...)
        }

        mat <- mat_list[[j]][[i]]

        image(as.numeric(rownames(mat)),
              as.numeric(colnames(mat)),
              mat,
              xlim = x$dat$grid$xrange,
              ylim = x$dat$grid$yrange,
              add = TRUE)

        if(plot_land){
          plot_land(x$dat$grid$xrange, x$dat$grid$yrange,
                    shift = ifelse(max(x$dat$grid$xrange) > 180, TRUE, FALSE))
        }

        {
          parts <- character(0)
          if (ncov > 1) {
            cov_nm <- names(x$dat$cov)[select_cov[i]]
            if (!is.null(cov_nm) && nzchar(cov_nm)) parts <- c(parts, cov_nm)
          }
          if (nsea > 1) parts <- c(parts, paste0("t = ", years[indi]))
          leg_lab <- paste(parts, collapse = ", ")
          if (nzchar(leg_lab))
            legend("topleft", legend = leg_lab, pch = NA, bg = "white", x.intersp = 0.1)
        }

        box(lwd=1.5)

      }
    }

  } else if (inherits(x, "admove_sim")){

    grid <- x$grid
    cov <- x$cov
    par <- x$par
    dat <- x$dat
    knots_tax <- dat$knots_tax
    funcs <- NULL

    if(is.null(par)) stop("No parameters provided! Use par = list() to specify parameters for taxis.")

    if(!add){
      if(!is.null(bg)){
        graphics::par(bg = bg)
      }
      plot(NA,
           xlim = grid$xrange,
           ylim = grid$yrange,
           xlab = xlab,
           ylab = ylab,
           ...)
      ## if(!is.null(bg)){
      ##     usr <- par("usr")
      ##     rect(usr[1], usr[3], usr[2], usr[4], col = bg, border = NA)
      ## }
    }

    par <- default_sim_par(par)
    cov <- .make_cov_list(cov)
    dat <- setup_data(cov = cov,
                       grid = grid,
                       trange = c(0,
                                  max(sapply(cov,
                                             function(x) dim(x)[3]))),
                       verbose = FALSE)
    conf <- default_conf(dat)
    funcs <- default_sim_funcs(dat, conf, par, funcs)

    ## uv.true <- t(apply(x$dat$xygrid, 1, function(xy)
    ##     funcs$tax(xy,NA)))

    get_true.pref <- .poly_fun(as.numeric(knots_tax),
                                       as.numeric(x$par_sim$alpha))

    i = 1
    pref_pred <- get_true.pref(as.numeric(cov[[i]]))

    if(plot_land){
      plot_land(grid$xrange, grid$yrange,
                shift = ifelse(max(grid$xrange) > 180, TRUE, FALSE))
    }

    dims <- dim(grid)
    image(
      matrix(pref_pred, dims$nx, dims$ny),
      xlim = grid$xrange,
      ylim = grid$yrange,
      add = TRUE)

    if(!add) box(lwd = 1.5)


  } else {

    grid <- x$grid
    cov <- x$cov
    par <- x$par
    dat <- x$dat
    knots_tax <- x$knots
    funcs <- NULL

    if(is.null(par)) stop("No parameters provided! Use par = list() to specify parameters for taxis.")

    if(!add){
      if(!is.null(bg)){
        graphics::par(bg = bg)
      }
      plot(NA,
           xlim = grid$xrange,
           ylim = grid$yrange,
           xlab = xlab,
           ylab = ylab,
           ...)
      ## if(!is.null(bg)){
      ##     usr <- par("usr")
      ##     rect(usr[1], usr[3], usr[2], usr[4], col = bg, border = NA)
      ## }
    }

    par <- default_sim_par(par)
    cov <- .make_cov_list(cov)
    dat <- setup_data(cov = cov,
                       grid = grid,
                       trange = c(0,
                                                   max(sapply(cov,
                                                              function(x) dim(x)[3]))),
                       verbose = FALSE)
    conf <- default_conf(dat)
    funcs <- default_sim_funcs(dat, conf, par, funcs)

    ## uv.true <- t(apply(x$dat$xygrid, 1, function(xy)
    ##     funcs$tax(xy,NA)))

          browser()



    get_true.pref <- .poly_fun(as.numeric(x$knots),
                                       as.numeric(x$par$alpha))

    i = 1
    pref_pred <- get_true.pref(as.numeric(cov[[i]]))

    if(plot_land){
      plot_land(grid$xrange, grid$yrange,
                shift = ifelse(max(grid$xrange) > 180, TRUE, FALSE))
    }

    dims <- dim(grid)

    image(
      matrix(pref_pred, dims$nx, dims$ny),
      xlim = grid$xrange,
      ylim = grid$yrange,
      add = TRUE)

    if(!add) box(lwd = 1.5)

  }
}
