
## Main functions ---------------------------------------------------------------------

##' Simulate a complete admove data set
##'
##' @description
##' Simulates a complete data set for \code{admove}, including a spatial grid,
##' optional covariate fields, simulated tags, and the corresponding
##' \code{admove_data} object required for model fitting.
##'
##' @param grid Optional spatial grid used for simulation. If not already an
##'   \code{admove_grid}, it is converted using [create_grid()].
##' @param cov Optional covariate fields. If \code{NULL}, or if
##'   \code{simulate_cov = TRUE}, covariates are simulated internally using
##'   [sim_cov()].
##' @param par Optional named list of simulation parameters. Missing parameters
##'   are filled using [default_sim_par()].
##' @param dat Optional \code{admove_data} object used as a template for the
##'   simulation. If supplied, missing inputs such as grid, covariates, time
##'   range, and spline knots are extracted from it.
##' @param conf Optional configuration list for the movement model.
##' @param fit Optional fitted \code{admove} model. If supplied, the simulation
##'   can reuse information from the fitted model, such as the original data
##'   object.
##' @param trange Numeric vector of length 2 giving the simulation time range.
##'   If \code{NULL}, a default range is used or extracted from \code{dat}.
##' @param dt Optional model time-step size.
##' @param simulate_cov Logical; if \code{TRUE}, covariate fields are simulated
##'   even when covariates are already available from \code{dat} or \code{fit}.
##' @param simple Logical; if \code{TRUE}, simple deterministic covariate fields
##'   are simulated instead of random fields.
##' @param nt Number of covariate fields (time steps) to simulate.
##' @param rho_t Temporal autocorrelation coefficient for simulated covariate
##'   fields.
##' @param sd Standard deviation used in the simulation of covariate fields.
##' @param h Parameter controlling the covariate precision or covariance
##'   structure.
##' @param nu Smoothness parameter of the Matérn covariance structure.
##' @param rho_s Spatial range parameter of the Matérn covariance structure.
##' @param delta Small positive value added for numerical stability in the
##'   precision matrix.
##' @param zrange Numeric vector of length 2 giving the target range of the
##'   simulated covariate values.
##' @param matern Logical; if \code{TRUE}, a Matérn-based covariance structure is
##'   used for covariate simulation.
##' @param sim_buffer Logical; if \code{TRUE}, the simulation grid is extended by
##'   a one-cell buffer before covariates are simulated.
##' @param knots_tax Optional knot locations for the taxis component.
##' @param knots_dif Optional knot locations for the diffusion component.
##' @param release_events Optional data frame or matrix of release events. If
##'   \code{NULL}, release events are simulated internally using
##'   [sim_release_events()].
##' @param n_release_events Number of release events to simulate if
##'   \code{release_events} is not supplied.
##' @param trange_rel Optional time range within which releases occur.
##' @param xrange_rel Optional x-range within which release locations are drawn.
##' @param yrange_rel Optional y-range within which release locations are drawn.
##' @param trange_rec Optional time range within which recapture or final
##'   observation times are drawn.
##' @param use_dtags Logical; if \code{TRUE}, data-storage tags are simulated.
##' @param use_ctags Logical; if \code{TRUE}, conventional mark-recapture tags
##'   are simulated.
##' @param use_stags Logical; if \code{TRUE}, mark-resight tags are simulated.
##' @param n_dtags Number of data-storage tags to simulate.
##' @param n_stags Number of mark-resight tags to simulate.
##' @param n_ctags Number of conventional mark-recapture tags to simulate.
##' @param n_resightings Integer vector giving the minimum and maximum number of
##'   resightings for mark-resight tags.
##' @param sim_engine Integer specifying the simulation engine: \code{1} for
##'   continuous-space simulation and \code{2} for CTMC-based grid simulation.
##' @param use_reject Logical; if \code{TRUE}, invalid simulated locations are
##'   rejected and redrawn where relevant.
##' @param n_reject Maximum number of rejection attempts used in rejection-based
##'   simulation steps.
##' @param target_dif_frac Target diffusion strength as a fraction of the
##'   characteristic spatial scale squared per unit time, used by
##'   [default_sim_par()].
##' @param target_tax_frac Target taxis strength as a fraction of the
##'   characteristic spatial scale per unit time, used by [default_sim_par()].
##' @param target_sdO_frac Target observation error as a fraction of the
##'   characteristic spatial scale, used by [default_sim_par()].
##' @param plot Logical; if \code{TRUE}, a summary plot of the simulated data is
##'   produced.
##' @param verbose Logical; if \code{TRUE}, informative messages are printed.
##'
##' @details
##' The function provides a convenient way to generate a complete simulation
##' setup for \code{admove}. Depending on the supplied inputs, it can reuse an
##' existing grid, covariates, or fitted model, or simulate these components
##' from scratch.
##'
##' Covariates are simulated with [sim_cov()] if needed. Tag release events are
##' either provided directly or generated with [sim_release_events()]. Tag data
##' for the requested tag types are then simulated using [sim_tags()]. Finally,
##' the simulated tags and covariates are combined into an \code{admove_data}
##' object suitable for model fitting.
##'
##' The returned object also includes default configuration, parameter, and map
##' objects for downstream fitting with [admove()].
##'
##' @return
##' An object of class \code{"admove_sim"} containing the simulated grid,
##' covariates, simulated tags, simulation parameters, \code{admove_data}
##' object, and default fitting components.
##'
##' @examples
##' sim <- sim_data()
##'
##' @export
sim_data <- function(grid = NULL,
                     cov = NULL,
                     par = NULL,
                     dat = NULL,
                     conf = NULL,
                     fit = NULL,
                     ## time
                     trange = NULL,
                     dt = NULL,
                     ## cov
                     simulate_cov = FALSE,
                     simple = FALSE,
                     nt = NULL,
                     rho_t = 0.85, sd = 2, h = 0.2, nu = 2,
                     rho_s = 0.8, delta = 0.1,
                     zrange = c(20, 28),
                     matern = TRUE,
                     sim_buffer = TRUE,
                     knots_tax = NULL,
                     knots_dif = NULL,
                     ## release events
                     release_events = NULL,
                     n_release_events = 10,
                     trange_rel = NULL,
                     xrange_rel = NULL,
                     yrange_rel = NULL,
                     ## tags
                     trange_rec = NULL,
                     use_dtags = TRUE,
                     use_ctags = TRUE,
                     use_stags = FALSE,
                     n_dtags = 10,
                     n_stags = 0,
                     n_ctags = 100,
                     n_resightings = c(1,5),
                     ## other
                     sim_engine = 1,
                     use_reject = FALSE,
                     n_reject = 100,
                     target_dif_frac = 1/300,
                     target_tax_frac = 1/10,
                     target_sdO_frac = 1/30,
                     plot = FALSE,
                     verbose = TRUE) {

  if (!is.null(fit) && .check_class(fit, "admove")) {
    dat0 <- fit$dat
    par0 <- get_par_est(fit$par, fit$map, fit$opt)
    conf0 <- fit$conf
    if (is.null(dat)) dat <- dat0
  }

  if (!is.null(dat) && .check_class(dat, "admove_data")) {
    grid0 <- dat$grid
    cov0 <- dat$cov
    trange0 <- dat$trange
    knots_tax0 <- dat$knots_tax
    knots_dif0 <- dat$knots_dif
    if (is.null(trange)) trange <- trange0
    if (is.null(grid)) grid <- grid0
    if (is.null(cov)) cov <- cov0
    if (is.null(knots_tax)) knots_tax <- knots_tax0
    if (is.null(knots_dif)) knots_dif <- knots_dif0
  }

  ## time
  if (is.null(trange)) trange <- c(0,1)
  if (is.null(trange_rel)) {
    trange_rel <- trange[1] + c(0, ifelse(trange[2] > 0.1, 0.1, 0))
  }
  if (is.null(trange_rec)) {
    trange_rec <- trange[2] - c(ifelse(trange[1] < (trange[2] - 0.1), 0.1, 0), 0)
  }

  ## space
  grid <- create_grid(grid)
  if (is.null(xrange_rel)) {
    xrange_rel <- grid$xrange + c(1,-1) * grid$cellsize[1]
  }

  if (is.null(yrange_rel)) {
    yrange_rel <- grid$yrange + c(1,-1) * grid$cellsize[2]
  }

  ## covariates
  cov <- .make_cov_list(cov)
  if (!is.null(cov) && is.null(nt)) nt <- dim(cov[[1]])[3]
  if (is.null(nt)) nt <- 1
  if (is.null(cov) || simulate_cov) {

    cov <- sim_cov(grid, nt = nt,
                   simple = simple,
                   rho_t = rho_t, sd = sd, h = h, nu = nu,
                   rho_s = rho_s, delta = delta,
                   zrange = zrange,
                   matern = matern,
                   sim_buffer = sim_buffer,
                   tref = list(origin = as.Date("2025-01-01"),
                               units = "year"))

  }

  if (is.null(dat)) {

    dat <- setup_data(cov = cov,
                      grid = grid,
                      knots_tax = knots_tax,
                      knots_dif = knots_dif,
                      trange = trange,
                      verbose = FALSE)

  }

  ## par
  par <- default_sim_par(par, dat,
                         target_dif_frac = target_dif_frac,
                         target_tax_frac =  target_tax_frac,
                         target_sdO_frac = target_sdO_frac)


  ## release events
  if (is.null(release_events)) {

    release_events <- sim_release_events(grid = grid,
                                         trange_rel = trange_rel,
                                         xrange_rel = xrange_rel,
                                         yrange_rel = yrange_rel,
                                         n_release_events = n_release_events,
                                         use_reject = use_reject,
                                         n_reject = n_reject)

  }


  if(use_ctags && n_ctags > 0){

    ctags_list <- sim_tags("c",
                           grid = grid,
                           par = par,
                           dat = dat,
                           n_tags = n_ctags,
                           trange = trange,
                           trange_rel = trange_rel,
                           trange_rec = trange_rec,
                           xrange_rel = xrange_rel,
                           yrange_rel = yrange_rel,
                           release_events = release_events,
                           sim_engine = sim_engine,
                           target_dif_frac = target_dif_frac,
                           target_tax_frac = target_tax_frac,
                           target_sdO_frac = target_sdO_frac,
                           sref = sref(grid),
                           tref = tref(cov)
                           )
    ctags <- ctags_list$tags

  } else {

    ctags <- NULL

  }


  if (use_dtags && n_dtags) {

    dtags_list <- sim_tags("d",
                           grid = grid,
                           par = par,
                           dat = dat,
                           n_tags = n_dtags,
                           trange = trange,
                           trange_rel = trange_rel,
                           trange_rec = trange_rec,
                           xrange_rel = xrange_rel,
                           yrange_rel = yrange_rel,
                           release_events = release_events,
                           sim_engine = sim_engine,
                           target_dif_frac = target_dif_frac,
                           target_tax_frac = target_tax_frac,
                           target_sdO_frac = target_sdO_frac,
                           sref = sref(grid),
                           tref = tref(cov)
                           )
    dtags <- dtags_list$tags

  } else {

    dtags <- NULL

  }

  if (use_stags && n_stags) {

    stags_list <- sim_tags("s",
                           grid = grid,
                           par = par,
                           dat = dat,
                           n_tags = n_stags,
                           trange = trange,
                           trange_rel = trange_rel,
                           trange_rec = trange_rec,
                           xrange_rel = xrange_rel,
                           yrange_rel = yrange_rel,
                           release_events = release_events,
                           n_resightings = n_resightings,
                           sim_engine = sim_engine,
                           target_dif_frac = target_dif_frac,
                           target_tax_frac = target_tax_frac,
                           target_sdO_frac = target_sdO_frac,
                           sref = sref(grid),
                           tref = tref(cov)
                           )
    stags <- stags_list$tags

  } else {

    stags <- NULL

  }

  ## combine tags
  tags <- as.data.frame(c(dtags, stags, ctags))
  sref(tags) <- sref(grid)
  tref(tags) <- tref(cov)

  dat <- setup_data(cov = cov,
                    grid = grid,
                    tags = tags,
                    knots_tax = knots_tax,
                    knots_dif = knots_dif,
                    trange = trange)

  res <- list()
  res$grid <- grid
  res$cov <- cov
  res$par_sim <- par
  res$tags <- tags
  res$dat <- dat

  res$conf <- default_conf(dat)
  res$par <- default_par(dat, res$conf)
  ## copy kappa as it is fixed
  res$par$logKappa <- res$par_sim$logKappa
  res$map <- default_map(dat, res$conf, res$par)

  res <- .add_class(res, "admove_sim")

  if (plot) plot(res)

  return(res)
}


##' Simulate release events
##'
##' @description
##' Simulates release locations and release times for tagging experiments on an
##' \code{admove_grid}.
##'
##' @param grid A grid object of class \code{"admove_grid"}, as returned by
##'   [create_grid()].
##' @param trange_rel Optional numeric vector of length 2 giving the time range
##'   within which release times are generated. If \code{NULL}, the default is
##'   \code{c(0, 1)}.
##' @param xrange_rel Optional numeric vector of length 2 giving the x-range
##'   within which release locations are generated. If \code{NULL}, the full
##'   x-range of \code{grid} is used.
##' @param yrange_rel Optional numeric vector of length 2 giving the y-range
##'   within which release locations are generated. If \code{NULL}, the full
##'   y-range of \code{grid} is used.
##' @param n_release_events Number of release events to simulate.
##' @param use_reject Logical; if \code{TRUE}, candidate release locations that
##'   fall into invalid grid cells are rejected and redrawn.
##' @param n_reject Maximum number of rejection attempts for each release event.
##'
##' @details
##' Release positions are drawn uniformly from the specified x- and y-ranges,
##' and release times are drawn uniformly from \code{trange_rel}. If
##' \code{use_reject = TRUE}, locations falling into \code{NA} cells of the grid
##' are rejected and resampled.
##'
##' @return
##' A numeric matrix with columns \code{x0}, \code{y0}, and \code{t0}, giving
##' the simulated release positions and release times.
##'
##' @examples
##' rel_events <- sim_release_events(create_grid())
##'
##' @export
sim_release_events <- function(grid,
                               trange_rel = NULL,
                               xrange_rel = NULL,
                               yrange_rel = NULL,
                               n_release_events = 10,
                               use_reject = TRUE,
                               n_reject = 100) {

  if (is.null(n_release_events) || is.na(n_release_events[1])) stop("Please provide a valid number of release events (n_release_events)!")
  n_release_events <- floor(n_release_events[1])

  if (!use_reject) n_reject <- 1

  .check_class(grid, "admove_grid")

  if (is.null(trange_rel)) {
    trange_rel <- c(0,1)
  }
  trange_rel <- sort(trange_rel)

  if (length(trange_rel) == 1) trange_rel <- rep(trange_rel, 2)

  if (is.null(xrange_rel)) {
    xrange_rel <- grid$xrange
  }
  xrange_rel <- sort(xrange_rel)

  if (length(xrange_rel) == 1) xrange_rel <- rep(xrange_rel, 2)

  if (is.null(yrange_rel)) {
    yrange_rel <- grid$yrange
  }
  yrange_rel <- sort(yrange_rel)

  if (length(yrange_rel) == 1) yrange_rel <- rep(yrange_rel, 2)

  release_events <- matrix(NA, n_release_events, 3)
  for (i in 1:n_release_events) {
    t0 <- runif(1, trange_rel[1], trange_rel[2])
    x0 <- runif(1, xrange_rel[1], xrange_rel[2])
    y0 <- runif(1, yrange_rel[1], yrange_rel[2])
    cnt <- 0
    is_invalid <- TRUE
    while (cnt < n_reject && is_invalid) {
      x0 <- runif(1, xrange_rel[1], xrange_rel[2])
      y0 <- runif(1, yrange_rel[1], yrange_rel[2])
      is_invalid <- is.na(grid$celltable[cbind(cut(x0, grid$xgr,
                                                   include.lowest = TRUE),
                                               cut(y0, grid$ygr,
                                                   include.lowest = TRUE))])
      cnt <- cnt + 1
    }
    release_events[i,] <- c(x0, y0, t0)
  }
  colnames(release_events) <- c("x0","y0","t0")

  return(release_events)
}


##' Simulate covariate fields
##'
##' @description
##' Simulates one or more spatial covariate fields on an \code{admove_grid},
##' optionally with temporal autocorrelation.
##'
##' @param grid Optional spatial grid. If not already an \code{admove_grid},
##'   it is converted using [create_grid()].
##' @param nt Number of covariate fields (time steps) to simulate.
##' @param simple Logical; if \code{TRUE}, a simple deterministic spatial field
##'   is simulated instead of a random field.
##' @param rho_t Temporal autocorrelation coefficient between successive
##'   covariate fields.
##' @param sd Standard deviation used in the simulation of the spatial random
##'   field.
##' @param h Parameter controlling the precision matrix or covariance structure.
##' @param nu Smoothness parameter of the Matérn covariance structure.
##' @param rho_s Spatial range parameter of the Matérn covariance structure. If
##'   \code{NULL}, a default is chosen from the grid resolution.
##' @param delta Small positive value added to the precision matrix diagonal to
##'   improve numerical stability.
##' @param zrange Numeric vector of length 2 giving the target range to which the
##'   simulated fields are rescaled.
##' @param trange Optional numeric vector of length 2 giving the time range
##'   covered by the simulated covariate series. If \code{NULL}, the range
##'   \code{c(0, nt - 1)} is used.
##' @param matern Logical; if \code{TRUE}, a Matérn-based covariance structure is
##'   used. Otherwise, a simpler neighbour-based precision matrix is used.
##' @param sim_buffer Logical; if \code{TRUE}, the grid is extended by a
##'   one-cell buffer before simulation.
##' @param tref Optional temporal reference attached to the returned covariate
##'   object.
##' @param verbose Logical; if \code{TRUE}, informative messages may be printed.
##'
##' @details
##' The first covariate field is simulated independently. If \code{nt > 1},
##' subsequent fields are generated using an AR(1)-type temporal dependence
##' structure with coefficient \code{rho_t}. All fields are then rescaled to the
##' interval given by \code{zrange}.
##'
##' The returned object is processed with [prep_cov()] and includes spatial and
##' temporal metadata.
##'
##' @return
##' An \code{admove_cov} object containing the simulated covariate fields.
##'
##' @examples
##' cov <- sim_cov()
##'
##' @export
sim_cov <- function(grid = NULL,
                    nt = 1,
                    simple = FALSE,
                    rho_t = 0.85,
                    sd = 2,
                    h = 0.2,
                    nu = 2,
                    rho_s = NULL, ## 0.8,
                    delta = 0.1,
                    zrange = c(20,28),
                    trange = NULL,
                    matern = TRUE,
                    sim_buffer = FALSE,
                    tref = NULL,
                    verbose = TRUE) {

  grid <- create_grid(grid)

  if (sim_buffer) {
    grid <- add_buffer(grid)
  }

  if (is.null(rho_s)) {
    rho_s <- mean(grid$cellsize) / 0.125
  }

  cov0 <- cov <- vector("list", nt)
  cov0[[1]] <- .get_cov_one(grid, sd, h, nu, rho_s, delta, matern, simple)
  if (nt > 1) {
    for (i in 2:nt) {
      cov0[[i]] <- rho_t * cov0[[i-1]] + sqrt(1 - rho_t^2) *
        .get_cov_one(grid, sd, h, nu, rho_s, delta, matern, simple)
    }
  }

  ## Rescale
  xcen <- x_centers(grid)
  ycen <- y_centers(grid)
  for (i in 1:nt) {
    cov[[i]] <- .rescale_cov(cov0, zrange, i)
    rownames(cov[[i]]) <- xcen
    colnames(cov[[i]]) <- ycen
  }

  if(is.null(trange)) trange <- c(0, nt-1)
  times <- sprintf("%.2f",
                   seq(trange[1],
                       trange[2],
                       length.out = nt))

  tref <- create_tref(tref$origin, tref$units, tref$period)

  cov <- prep_cov(cov,
                  x_centers = xcen,
                  y_centers = ycen,
                  times = times,
                  sref = sref(grid),
                  tref = tref,
                  verbose = FALSE)

  return(cov)
}


##' Simulate tagging data
##'
##' @description
##' Simulates tagging data under the \code{admove} movement model for one of the
##' supported tag types: data-storage tags, mark-resight tags, or conventional
##' mark-recapture tags.
##'
##' @param tag_type Character string specifying the tag type to simulate.
##'   Supported values are \code{"d"} or \code{"dtags"} for data-storage tags,
##'   \code{"s"} or \code{"stags"} for mark-resight tags, and \code{"c"} or
##'   \code{"ctags"} for conventional mark-recapture tags.
##' @param grid Optional spatial grid used for simulation. If not already of
##'   class \code{"admove_grid"}, it is converted using [create_grid()].
##' @param cov Optional covariate fields used to define spatially and temporally
##'   varying movement rates.
##' @param par Optional named list of simulation parameters. Missing parameters
##'   are filled using [default_sim_par()].
##' @param dat Optional \code{admove_data} object. If \code{NULL}, a default data
##'   object is constructed internally from the supplied inputs.
##' @param conf Optional configuration list controlling the movement model. If
##'   \code{NULL}, [default_conf()] is used.
##' @param n_tags Number of tags to simulate.
##' @param n_resightings Integer vector giving the minimum and maximum number of
##'   resightings for mark-resight tags.
##' @param trange Numeric vector of length 2 giving the simulation time range.
##'   If \code{NULL}, it is derived from the covariate time dimension when
##'   possible, or defaults to \code{c(0, 1)}.
##' @param dt_tags Optional simulation time step for tag trajectories. If
##'   \code{NULL}, a default value of \code{0.1} is used.
##' @param trange_rel Optional time range within which release times are
##'   generated.
##' @param xrange_rel Optional x-range within which release positions are
##'   generated.
##' @param yrange_rel Optional y-range within which release positions are
##'   generated.
##' @param n_release_events Number of release events to simulate if
##'   \code{release_events} is not supplied.
##' @param release_events Optional data frame specifying release events. If
##'   \code{NULL}, release events are generated internally.
##' @param trange_rec Optional time range within which recapture or final
##'   observation times are generated.
##' @param knots_tax Optional knot locations for the taxis component.
##' @param knots_dif Optional knot locations for the diffusion component.
##' @param funcs Optional named list of simulation functions. If \code{NULL},
##'   defaults are created with [default_sim_funcs()].
##' @param use_reject Logical; if \code{TRUE}, invalid movement proposals are
##'   rejected and resampled.
##' @param n_reject Maximum number of rejection attempts per step if
##'   \code{use_reject = TRUE}.
##' @param sim_engine Integer specifying the simulation engine: \code{1} for
##'   continuous-space simulation and \code{2} for CTMC-based grid simulation.
##' @param ctmc_method Integer controlling the matrix-exponential method used for
##'   CTMC simulation.
##' @param sref Optional spatial reference to attach to the simulated data.
##' @param tref Optional temporal reference to attach to the simulated data.
##' @param target_dif_frac Target diffusion strength as a fraction of the
##'   characteristic spatial scale squared per unit time, used by
##'   [default_sim_par()].
##' @param target_tax_frac Target taxis strength as a fraction of the
##'   characteristic spatial scale per unit time, used by [default_sim_par()].
##' @param target_sdO_frac Target observation error as a fraction of the
##'   characteristic spatial scale, used by [default_sim_par()].
##' @param add_obs_unc Add observation uncertainty to tag locations? By default
##'   (NULL), observation uncertainty is added to archival tags (tag_type =
##'   "d"), but not to other tags.
##' @param plot Logical; if \code{TRUE}, the simulated tags are plotted.
##' @param plot_land Logical; if \code{TRUE}, land masses are added to the plot.
##' @param verbose Logical; if \code{TRUE}, informative messages are printed.
##'
##' @details
##' The function first sets up the spatial grid, time range, covariates, and
##' model configuration required for simulation. If no data object is supplied,
##' one is created internally from the provided inputs.
##'
##' Tag trajectories are then simulated from generated or user-supplied release
##' events. The full simulated trajectories are retained for data-storage tags.
##' For mark-resight tags, only the release and a subset of subsequent
##' observations are kept. For conventional mark-recapture tags, only the
##' release and final observation are retained.
##'
##' The returned object contains both the simulated tags and the associated
##' simulation setup, including the grid, covariates, simulation parameters, and
##' internally constructed \code{admove_data} object.
##'
##' @return
##' An object of class \code{"admove_sim"} containing the simulated tags, grid,
##' covariates, simulation parameters, data object, and default fitting
##' components.
##'
##' @examples
##' data(skjepo)
##' sim <- sim_tags("ctags", skjepo$grid)
##'
##' @export
sim_tags <- function(tag_type,
                     grid = NULL,
                     cov = NULL,
                     par = NULL,
                     dat = NULL,
                     conf = NULL,
                     n_tags = 1,
                     n_resightings = c(1,5),
                     trange = NULL,
                     dt_tags = NULL,
                     trange_rel = NULL,
                     xrange_rel = NULL,
                     yrange_rel = NULL,
                     n_release_events = 1,
                     release_events = NULL,
                     trange_rec = NULL,
                     knots_tax = NULL,
                     knots_dif = NULL,
                     funcs = NULL,
                     use_reject = FALSE,
                     n_reject = 20,
                     sim_engine = 1,
                     ctmc_method = 2,
                     sref = NULL,
                     tref = NULL,
                     target_dif_frac = 1/300,
                     target_tax_frac = 1/10,
                     target_sdO_frac = 1/30,
                     add_obs_unc = NULL,
                     plot = FALSE,
                     plot_land = FALSE,
                     verbose = TRUE) {

  if (n_tags == 0) return(NULL)


  tag_type <- .get_tag_type(tag_type)
  if (is.na(tag_type)) {
    stop("Tag type not implemented, use tag_type = 'd', 's', or 'c' for data-storage, mark-resight, and mark-recapture tags, respectively.")
  }

  if (is.null(n_tags) || is.na(n_tags[1])) stop("Please provide a valid number of tags (n_tags)!")
  n_tags <- floor(n_tags[1])

  if (length(n_resightings) == 1) {
    n_resightings <- rep(n_resightings, 2)
  }

  cov <- .make_cov_list(cov)

  ## Dimensions
  if (is.null(trange) && is.null(cov)) {
    trange <- c(0,1)
  }else if(is.null(trange) && !is.null(cov)){
    if (!is.null(attributes(cov[[1]])$dimnames[[3]])) {
      ## stop("Implement to simulate tags from dates in cov, but dates could have any format.")
      trange <- range(as.numeric(attributes(cov[[1]])$dimnames[[3]]))
      if (trange[1] == trange[2]) trange[2] <- trange[1] + 1
    }else{
      trange <- c(0, max(sapply(cov, function(x) dim(x)[3])))
      if(verbose) warning("Time range extracted from number of covariate fields.")
    }
  }

  if (is.null(trange_rel)) {
    trange_rel <- trange
  }
  if (length(trange_rel) == 1) trange_rel <- rep(trange_rel, 2)

  if (trange_rel[2] > trange[2]) stop("trange_rel[2] > trange[2]!")

  if (is.null(trange_rec)) {
    trange_rec <- trange
  }
  if (length(trange_rec) == 1) trange_rec <- rep(trange_rec, 2)

  grid <- create_grid(grid)

  xrange <- grid$xrange
  yrange <- grid$yrange

  if (is.null(xrange_rel)) {
    xrange_rel <- xrange
  }
  if (length(xrange_rel) == 1) xrange_rel <- rep(xrange_rel, 2)

  if (is.null(yrange_rel)) {
    yrange_rel <- yrange
  }
  if (length(yrange_rel) == 1) yrange_rel <- rep(yrange_rel, 2)

  ## Checks
  trange_rel <- sort(trange_rel)
  trange_rec <- sort(trange_rec)
  xrange_rel <- sort(xrange_rel)
  yrange_rel <- sort(yrange_rel)

  if (is.null(sref)) {
    sref_target <- sref(grid)
  } else {
    sref_target <- create_sref(sref$crs,
                                 sref$units,
                                 sref$crs_scale)
  }

  if (is.null(tref)) {
    tref_target <- tref(cov)
  } else {
    tref_target <- create_tref(tref$origin,
                               tref$units,
                               tref$period)
  }

  ## Setup default data and conf
  if (is.null(dat)) {
    dat <- setup_data(cov = cov,
                      grid = grid,
                      knots_tax = knots_tax,
                      knots_dif = knots_dif,
                      trange = trange,
                      sref = sref_target,
                      tref = tref_target,
                      transform_sref = TRUE,
                      shift_tref = TRUE,
                      verbose = FALSE)
  }


  if (is.null(conf)) {
    conf <- default_conf(dat)
  }

  if (is.null(dat$cov)) {
    conf$use_taxis <- FALSE
    conf$use_advection <- FALSE
  }

  ## Parameters
  par <- default_sim_par(par, dat,
                         target_dif_frac = target_dif_frac,
                         target_tax_frac =  target_tax_frac,
                         target_sdO_frac = target_sdO_frac)
  par_units <- attr(par, "units")

  if (is.null(dt_tags)) {
    dt_tags <- 0.1 ## diff(trange) / 10
  }

  ## Functions
  funcs <- default_sim_funcs(dat, conf, par, funcs)

  ## Release events
  if (is.null(release_events)) {
    release_events <- sim_release_events(grid = grid,
                                         trange_rel = trange_rel,
                                         xrange_rel = xrange_rel,
                                         yrange_rel = yrange_rel,
                                         n_release_events = n_release_events,
                                         use_reject = use_reject,
                                         n_reject = n_reject)
  } else {
    n_release_events <- nrow(release_events)
  }

  nextTo <- get_neighbours(dat$grid)
  next_dist <- c(dat$grid$cellsize[1], dat$grid$cellsize[1],
                 dat$grid$cellsize[2], dat$grid$cellsize[2])
  xcen <- x_centers(dat$grid)
  ycen <- y_centers(dat$grid)

  ## to keep release events unique
  id <- .get_random_id()

  res_list <- vector("list", n_tags)
  if (n_tags < n_release_events) n_release_events <- n_tags
  n_by_rel_event <- ceiling(n_tags / n_release_events)
  count <- 1

  for (i in 1:n_release_events) {
    for (j in 1:n_by_rel_event) {
      x0 <- release_events[i,"x0"]
      y0 <- release_events[i,"y0"]
      t0 <- release_events[i,"t0"]
      trec1 <- ifelse(trange_rec[1] < t0, t0+dt_tags, trange_rec[1])
      t1 <- round(runif(1, trec1,
                        ifelse(trange_rec[2]-dt_tags < trec1,
                               trec1, trange_rec[2]-dt_tags))/dt_tags) * dt_tags

      tmp <- .sim_one_tag(conf,
                         funcs,
                         par,
                         x0, y0, t0, t1,
                         dt = dt_tags,
                         id = paste0(id,"-",count),
                         tag_type = tag_type,
                         xygrid = dat$grid$xygrid,
                         nextTo = nextTo,
                         next_dist = next_dist,
                         xgr = dat$grid$xgr,
                         ygr = dat$grid$ygr,
                         celltable = dat$grid$celltable,
                         xcen = xcen,
                         ycen = ycen,
                         sim_engine = sim_engine,
                         use_reject = use_reject,
                         n_reject = n_reject,
                         ctmc_method = ctmc_method)

      if (tag_type == "d") {
        res_list[[count]] <- tmp
      } else if (tag_type == "s") {
        if (nrow(tmp) == 1) next()
        ind <- try(sample(2:nrow(tmp),
                          round(runif(1, n_resightings[1], n_resightings[2]))),
                   silent = TRUE)
        if (inherits(ind, "try-error")) {
          ind <- sample(2:nrow(tmp), 1)
        }
        res_list[[count]] <- tmp[c(1,ind),]
      } else if (tag_type == "c") {
        if (nrow(tmp) == 1) next()
        res_list[[count]] <- tmp[c(1,nrow(tmp)),]
      }
      count <- count + 1
    }
  }

  tags <- as.data.frame(do.call(rbind, res_list))
  if (nrow(tags) == 0) stop("No tags left after simulation (and subsampling). Run again, might just be a random failure or revise the settings of your simulation.")
  tags$tag_type <- tag_type

  tags <- .add_class(tags, "admove_tags")
  tags <- add_sref(tags, sref(dat))
  tags <- add_tref(tags, tref(dat))

  res <- list()
  res$grid <- grid
  res$cov <- cov
  res$par_sim <- par
  res$tags <- tags
  res$dat <- dat

  res$conf <- default_conf(dat)
  res$par <- default_par(dat, res$conf)
  res$map <- default_map(dat, res$conf, res$par)

  res <- .add_class(res, "admove_sim")

  if (plot) plot(res$tags, plot_land = plot_land)

  return(res)
}



##' Plot a summary of simulated admove data
##'
##' @description
##' Produces a set of summary plots for an object of class \code{"admove_sim"}.
##' Depending on the contents of the simulation object, the plots may include a
##' covariate field, habitat preference function, taxis, diffusion, and the
##' simulated tag tracks.
##'
##' @param x An object of class \code{"admove_sim"}, as returned by
##'   [sim_data()].
##' @param auto_layout Logical; if \code{TRUE}, the function automatically sets
##'   and restores graphical parameters.
##' @param plot_land Logical; if \code{TRUE}, land masses are added to spatial
##'   plots.
##' @param cor_taxis Optional scaling factor for taxis arrows.
##' @param cor_diffusion Optional scaling factor for diffusion arrows.
##' @param asp Positive numeric value giving the target aspect ratio
##'   (columns / rows) of the plot arrangement.
##' @param by_tag_type Logical; if \code{TRUE}, simulated tracks are plotted
##'   separately for each tag type. If \code{FALSE}, all tracks are shown in a
##'   single panel.
##' @param ... Additional arguments passed to lower-level plotting functions.
##'
##' @details
##' When \code{auto_layout = TRUE}, the function arranges multiple panels in a
##' suitable plotting layout. The exact set of panels depends on the simulation
##' object and whether simulated tag data are available.
##'
##' If simulated tags are present and \code{by_tag_type = TRUE}, a separate plot
##' is produced for each tag type.
##'
##' @return
##' Invisibly returns \code{NULL}. The function is called for its plotting side
##' effects.
##'
##' @name plot_sim
##' @export
plot_sim <- function(x,
                     auto_layout = TRUE,
                     plot_land = FALSE,
                     cor_taxis = NULL,
                     cor_diffusion = NULL,
                     asp = 2,
                     by_tag_type = TRUE,
                     ...) {


  .check_class(x, "admove_sim")

  if(auto_layout){
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar))
    n <- 4 + length(unique(x$tags$tag_type))
    par(mfrow = n2mfrow(n, asp = asp),
        mar = c(4,4,1,1), oma = c(1,1,1,1))
  }

  i = 1
  if(inherits(x$cov, "list")) {
    tmp <- x$cov[[1]][,,1, drop=FALSE]
  } else {
    tmp <- x$cov[,,1, drop=FALSE]
  }
  plot_cov(tmp, auto_layout = FALSE, main = "", plot_land = plot_land, ...)
  add_lab(LETTERS[i])
  i = i + 1
  plot_pref_func(x, auto_layout = FALSE, main = "", ...)
  add_lab(LETTERS[i])
  i = i + 1
  ## plot_pref_grid(x, auto_layout = TRUE, main = "")
  plot_taxis(x, auto_layout = FALSE, main = "", cor = cor_taxis,
             plot_land = plot_land, ...)
  add_lab(LETTERS[i])
  i = i + 1
  plot_diffusion(x, auto_layout = FALSE, main = "", cor = cor_diffusion,
                 plot_land = plot_land, ...)
  add_lab(LETTERS[i])
  i = i + 1

  if (!is.null(x$tags)) {
    if (by_tag_type) {
      tag_types <- unique(x$tags$tag_type)
      n_tag_types <- length(tag_types)
      tags <- x$tags
      for (j in seq_len(n_tag_types)) {
        x$tags <- tags[tags$tag_type == tag_types[j],]
        plot_grid(x, auto_layout = FALSE, main = "",
                  plot_land = plot_land, labels = FALSE,
                  plot_grid = FALSE, plot_bg = FALSE)
        plot_tags(x, auto_layout = FALSE, add = TRUE, ...)
        add_lab(LETTERS[i])
        i = i + 1
      }
    } else {
        plot_grid(x, auto_layout = FALSE, main = "",
                  plot_land = plot_land, labels = FALSE,
                  plot_grid = FALSE, plot_bg = FALSE)
      plot_tags(x, auto_layout = FALSE, add = TRUE, ...)
      add_lab(LETTERS[i])
      i = i + 1
    }
  }

}


##' Default parameter values for simulation
##'
##' @description
##' Creates a named list of simulation parameter values for \code{admove}. The
##' defaults are chosen to give movement rates and observation error that are
##' scaled to the spatial and temporal extent of the supplied data, where
##' available.
##'
##' @param par An optional named list of parameter values used to overwrite the
##'   corresponding defaults.
##' @param dat An optional data object of class \code{"admove_data"} used to
##'   scale the default simulation parameters to the spatial domain, time range,
##'   and available covariates.
##' @param time_unit Optional time unit label attached to the returned parameter
##'   object when \code{dat} is not provided.
##' @param alpha_template Numeric vector giving the template spline coefficients
##'   used to construct default taxis parameters.
##' @param target_tax_per_time Optional target taxis speed in distance units per
##'   time unit. If \code{NULL}, \code{target_tax_frac} is used instead.
##' @param target_tax_frac Target taxis speed as a fraction of the characteristic
##'   spatial scale per unit time. Used only if \code{target_tax_per_time} is
##'   \code{NULL}.
##' @param target_dif_per_time Optional target diffusion coefficient in distance
##'   squared per time unit. If \code{NULL}, \code{target_dif_frac} is used
##'   instead.
##' @param target_dif_frac Target diffusion coefficient as a fraction of the
##'   squared characteristic spatial scale per unit time. Used only if
##'   \code{target_dif_per_time} is \code{NULL}.
##' @param target_sdO_frac Target observation error as a fraction of the
##'   characteristic spatial scale.
##'
##' @details
##' If \code{dat} is provided, default values are scaled to the spatial extent
##' of the grid and the total time range of the data. If covariates are
##' available, the taxis spline coefficients are additionally scaled so that the
##' resulting taxis strength is approximately consistent with the requested
##' target movement rate.
##'
##' The returned parameter list includes defaults for taxis (\code{alpha}),
##' diffusion (\code{beta}), advection (\code{gamma}), taxis scaling
##' (\code{logKappa}), and observation error (\code{logSdO}). User-supplied
##' values in \code{par} overwrite the corresponding defaults.
##'
##' @return
##' A named list of simulation parameter values. The returned object also
##' carries a \code{"units"} attribute containing the distance and time units,
##' when available.
##'
##' @export
default_sim_par <- function(par = NULL,
                            dat = NULL,
                            time_unit = NULL,
                            alpha_template = c(0, 5, 4),  ## c(0, 4, 3),
                            target_tax_per_time = NULL,
                            target_tax_frac =  1/10,
                            target_dif_per_time = NULL,
                            target_dif_frac =  1/300,
                            target_sdO_frac = 1/30) {

  D_default <- 0.05 * alpha_template
  alpha_default <- c(0,5,1)
  kappa_default <- 1


  if (!is.null(dat)) {
    Lx <- diff(dat$grid$xrange)
    Ly <- diff(dat$grid$yrange)
    dist_diff <- sqrt(Lx * Ly)
    time_diff <- diff(dat$trange)
    dist_unit <- units_space(dat)
    time_unit <- units_time(dat)
  } else {
    dist_diff <- 1
    time_diff <- 1
    dist_unit <- NULL
    time_unit <- NULL
  }

  ## diffusion target (e.g. km^2 / day)
  if (is.null(target_dif_per_time)) {
    D_target <- target_dif_frac * dist_diff^2 / time_diff
  } else {
    D_target <- target_dif_per_time
  }

  if (is.na(D_target) || is.null(D_target) ||
        !is.numeric(D_target)) {
    D_target <- D_default
  }

  ## scaling of dH/dx
  kappa_target <- 1 * dist_diff^2 / time_diff

  if (!is.null(dat$cov)) {

    field <- dat$cov[[1]][,,1]
    xyr <- .get_cov_xyrange(dat$cov)
    xr <- xyr$xr[1,]
    yr <- xyr$yr[1,]

    onedx <- (xr[2] - xr[1]) / (nrow(field) - 1)
    onedy <- (yr[2] - yr[1]) / (ncol(field) - 1)

    dx <- (field[3:nrow(field), ] - field[1:(nrow(field)-2), ]) / (2 * onedx)
    dy <- (field[, 3:ncol(field)] - field[, 1:(ncol(field)-2)]) / (2 * onedy)

    gradC_mag <- sqrt(dx[, 2:(ncol(field)-1)]^2 +
                        dy[2:(nrow(field)-1), ]^2)

    field_sub <- field[2:(nrow(field)-1), 2:(ncol(field)-1)]
    field_vec <- as.numeric(field_sub)

    ## target taxis speed (distance / time)
    if (is.null(target_tax_per_time)) {
      tax_target <- target_tax_frac * dist_diff / time_diff
    } else {
      tax_target <- target_tax_per_time
    }

    ## derivative for template at scale 1
    dS1 <- .poly_fun(as.numeric(dat$knots_tax[,1]),
                    alpha_template,
                    deriv = TRUE)

    ## typical |∇h| magnitude for scale 1
    gh1 <- abs(dS1(field_vec)) * gradC_mag
    gh1 <- gh1[is.finite(gh1)]

    eps <- 1e-12
    m <- if (length(gh1)) median(gh1) else 0
    if (!is.finite(m) || m < eps) {
      x_scale <- 0
    } else {
      ## solve for x so that kappa * typical(|∇h|) matches tax_target
      x_scale <- tax_target / (kappa_target * m)
      ## optional clamp to avoid crazy values
      x_scale <- max(min(x_scale, 1e4), -1e4)
    }

    alpha <- x_scale * alpha_template

  } else {

    alpha <- alpha_default

  }

  sdO <- dist_diff * target_sdO_frac

  par_out <- list(alpha = array(alpha, dim = c(3,1,1)),         ## dimensionless
                  beta = array(log(D_target), dim = c(1,1,1)),  ## distance^2 / time
                  gamma = array(0, dim = c(2,1,1)),             ## dimensionless
                  logKappa = log(kappa_target),                 ## distance^2 / time
                  logSdO = matrix(log(sdO),2,3))                ## distance

  if(!is.null(par)){
    for(i in 1:length(par)){
      par_out[names(par)[i]] <- par[names(par)[i]]
    }
  }

  par <- par_out

  attr(par, "units") <- list(
    distance = dist_unit,
    time = time_unit
  )

  return(par)
}



##' Default simulation functions for admove
##'
##' @description
##' Constructs default function objects used for simulation in \code{admove},
##' including functions for taxis, diffusion, the spatial gradient of
##' diffusion, and advection.
##'
##' @param dat An \code{admove_data} object, as produced by [setup_data()].
##' @param conf A configuration list, typically produced by [default_conf()].
##' @param par A parameter list, typically produced by [default_par()] or
##'   [default_sim_par()].
##' @param funcs An optional named list of user-supplied functions that
##'   overwrite the corresponding defaults. Allowed names are \code{"tax"},
##'   \code{"dif"}, \code{"ddif"}, and \code{"adv"}.
##'
##' @details
##' If covariate data are available in \code{dat}, the default functions are
##' constructed from the parameter values and interpolated covariate fields.
##' This allows taxis, diffusion, and advection to vary in space and time.
##'
##' If no covariate data are available, simple default functions are returned:
##' diffusion is constant, and taxis, advection, and the gradient of diffusion
##' are set to zero.
##'
##' User-supplied functions in \code{funcs} replace the corresponding default
##' functions in the returned list.
##'
##' @return
##' A named list of simulation functions with elements \code{tax},
##' \code{dif}, \code{ddif}, and \code{adv}.
##'
##' @export
default_sim_funcs <- function(dat, conf, par, funcs = NULL) {

  if (!is.null(dat$cov)) {

    ## Make preference functions --------------------------
    pref_funcs <- .make_pref_funcs(par$alpha, par$beta, par$gamma,
                                  dat$knots_tax, dat$knots_dif)

    ## Local interpolation --------------------------------
    liv <- .get_liv(dat$cov)

    ## Setup habi objects ---------------------------------
    habi_dif <- .make_habi(liv, dat$xrange_cov,
                          dat$yrange_cov, dat$time_cov,
                          pref_funcs$dif, pref_funcs$ddif,
                          dat$time_spline, period(dat),
                          conf$seasonal_cov,
                          conf$seasonal_spline)
    habi_tax <- .make_habi(liv, dat$xrange_cov,
                          dat$yrange_cov, dat$time_cov,
                          pref_funcs$tax, pref_funcs$dtax,
                          dat$time_spline, period(dat),
                          conf$seasonal_cov,
                          conf$seasonal_spline)
    habi_adv_x <- .make_habi(liv, dat$xrange_cov,
                            dat$yrange_cov, dat$time_cov,
                            pref_funcs$adv_x, pref_funcs$dadv_x,
                            dat$time_spline, period(dat),
                            conf$seasonal_cov,
                            conf$seasonal_spline)
    habi_adv_y <- .make_habi(liv, dat$xrange_cov,
                            dat$yrange_cov, dat$time_cov,
                            pref_funcs$adv_y, pref_funcs$dadv_y,
                            dat$time_spline, period(dat),
                            conf$seasonal_cov,
                            conf$seasonal_spline)

    dif_fun <- function(xy, t){
      habi_dif$val(xy, t)
    }

    ddif_fun <- function(xy, t){
      habi_dif$grad(xy, t)
    }

    tax_fun <- function(xy, t){
      habi_tax$grad(xy, t)
    }

    adv_fun <- function(xy, t) {
      c(habi_adv_x$val(xy, t), habi_adv_y$val(xy, t))
    }

    res <- list(dif = dif_fun,
                ddif = ddif_fun,
                tax = tax_fun)

    if (conf$use_advection) {
      res$adv <- adv_fun
    }


  } else {

    ## diffusion only if no cov info
    dif_fun <- function(xy, t){par$beta[1,1,1]}
    ddif_fun <- function(xy, t){c(0,0)}
    tax_fun <- function(xy, t) {c(0,0)}
    adv_fun <- function(xy, t) {c(0,0)}

  }

  funcs_out <- list(tax = tax_fun,
                    dif = dif_fun,
                    ddif = ddif_fun,
                    adv = adv_fun)

  if (!is.null(funcs)) {
    for (nm in names(funcs)) {
      if (is.function(funcs[[nm]])) {
        funcs_out[[nm]] <- funcs[[nm]]
      }
    }
  }

  return(funcs_out)
}



## Internal functions --------------------------------------------------------------

##' Rescale covariate fields
##'
##' @description
##' Rescales a covariate field, or one element of a list of covariate fields, to
##' a specified numeric range.
##'
##' @param cov A covariate field or a list of covariate fields.
##' @param zrange Numeric vector of length 2 giving the target range.
##' @param i Optional index specifying which element of \code{cov} to rescale if
##'   \code{cov} is a list.
##'
##' @return
##' A rescaled covariate field with values mapped to \code{zrange}.
##'
##' @keywords internal
.rescale_cov <- function(cov, zrange, i = NULL){
  cov_range <- range(unlist(cov), na.rm = TRUE)
  covi <- if(!is.null(i)) cov[[i]] else cov
  res <- (covi - cov_range[1]) /
    (cov_range[2] - cov_range[1]) * (zrange[2] - zrange[1]) + zrange[1]
  return(res)
}



##' Simulate a single covariate field
##'
##' @description
##' Simulates a single spatial covariate field on an \code{admove_grid}.
##'
##' @param grid A grid object of class \code{"admove_grid"}, as returned by
##'   [create_grid()].
##' @param sd Standard deviation used in the simulation of the random field.
##' @param h Parameter controlling the precision matrix or covariance structure.
##' @param nu Smoothness parameter of the Matérn covariance structure.
##' @param rho Spatial range parameter of the Matérn covariance structure.
##' @param delta Small positive value added to the precision matrix diagonal to
##'   improve numerical stability.
##' @param matern Logical; if \code{TRUE}, a Matérn-based covariance structure is
##'   used. Otherwise, a simpler neighbour-based precision matrix is used.
##' @param simple Logical; if \code{TRUE}, a simple deterministic spatial field
##'   is generated instead of a random field.
##'
##' @details
##' If \code{simple = TRUE}, the simulated field is a smooth radial surface
##' centred on the middle of the grid. Otherwise, a random field is generated
##' from a precision matrix defined by [.get_precision_matrix()].
##'
##' @return
##' A matrix containing one simulated covariate field on the grid.
##'
##' @keywords internal
.get_cov_one <- function(grid, sd, h, nu, rho, delta, matern, simple) {

  dims <- dim(grid)
  rf_smooth <- matrix(NA_real_,
                      nrow = dims["nx"],
                      ncol = dims["ny"])

  if (simple) {

    x0 <- mean(grid$xrange)
    y0 <- mean(grid$yrange)
    ## half of the smaller map width/height (distance from center to edge)
    r_edge <- 0.5 * min(diff(grid$xrange), diff(grid$yrange))

    p_edge <- 0.5
    sigma <- r_edge / sqrt(2 * log(1 / p_edge))
    S <- exp(-((grid$xygrid[,1] - x0)^2 + (grid$xygrid[,2] - y0)^2) / (2 * sigma^2))

  } else {

    ## Generate a random field
    rf <- rnorm(nrow(grid$xygrid), 0, sd = sd)

    ## GMRF
    Q <- .get_precision_matrix(grid, h, nu, rho, delta, matern)
    L <- chol(Q)
    S <- solve(L, rf)

  }

  ind <- which(!is.na(grid$celltable))
  rf_smooth[ind] <- as.numeric(S)

  stopifnot(is.matrix(rf_smooth))

  return(rf_smooth)
}


##' Construct a precision matrix for a grid
##'
##' @description
##' Constructs a spatial precision matrix for an \code{admove_grid}, either from
##' a Matérn covariance structure or from a simple neighbour-based structure.
##'
##' @param grid A grid object of class \code{"admove_grid"}, as returned by
##'   [create_grid()].
##' @param h Distance argument used internally in the Matérn covariance
##'   structure.
##' @param nu Smoothness parameter of the Matérn covariance structure.
##' @param rho Spatial range parameter of the Matérn covariance structure.
##' @param delta Small positive value added to the diagonal to ensure numerical
##'   stability and positive definiteness.
##' @param matern Logical; if \code{TRUE}, a Matérn-based covariance structure is
##'   used. Otherwise, a neighbour-based graph Laplacian is constructed.
##'
##' @details
##' When \code{matern = TRUE}, pairwise distances between grid-cell centres are
##' used to build a Matérn covariance matrix, which is then inverted to obtain a
##' precision matrix. When \code{matern = FALSE}, a sparse neighbour-based
##' precision matrix is constructed from the grid adjacency structure.
##'
##' In both cases, \code{delta * I} is added to the diagonal.
##'
##' @return
##' A precision matrix for the active cells of the grid.
##'
##' @keywords internal
.get_precision_matrix <- function(grid, h, nu, rho, delta,
                                 matern = TRUE) {
  n <- nrow(grid$xygrid)

  if(matern){

    ## Compute the pairwise distance matrix between all points in the grid
    dist_matrix <- as.matrix(dist(grid$xygrid))

    ## Define the Matérn covariance function for distances
    matern_covariance <- function(h, nu, rho) {
      if(h == 0) return(1)  ## Variance at distance 0
      scale_factor <- (2^(1 - nu)) / gamma(nu)
      distance_factor <- (sqrt(2 * nu) * h / rho)^nu
      bessel_part <- besselK(sqrt(2 * nu) * h / rho, nu)
      return(scale_factor * distance_factor * bessel_part)
    }

    ## Construct the covariance matrix using the Matérn function
    C <- matrix(0, n, n)
    for (i in 1:n) {
      for (j in 1:n) {
        C[i, j] <- matern_covariance(dist_matrix[i, j], nu, rho)
      }
    }

    ## Invert the covariance matrix to get the precision matrix Q
    Q0 <- solve(C)

  }else{

    nextTo <- get_neighbours(grid)
    Q0 <- Matrix::sparseMatrix(1:n, 1:n, x=0, dims = c(n, n))
    diag(Q0) <- rowSums(!is.na(nextTo[,-1]), na.rm=TRUE)
    for(i in 1:nrow(nextTo)){
      Q0[i,nextTo[,-1][i,!is.na(nextTo[,-1][i,])]] <- -1
    }

  }

  ## Add delta * I to ensure positive definiteness
  I <- Matrix::sparseMatrix(1:n, 1:n, x=0, dims = c(n, n))
  diag(I) <- 1
  Q <- Q0 + delta * I

  return(Q)
}


##' Simulate a single tag track
##'
##' @description
##' Simulates the movement trajectory of a single tagged individual under the
##' \code{admove} movement model. Depending on \code{sim_engine}, movement is
##' simulated either in continuous space using an Euler-type update or on a
##' discrete spatial grid using a continuous-time Markov chain (CTMC)
##' formulation.
##'
##' @param conf A configuration list controlling which movement components are
##'   active, typically produced by [default_conf()].
##' @param funcs A named list of simulation functions, typically produced by
##'   [default_sim_funcs()]. These functions define taxis, diffusion, diffusion
##'   gradient, and advection.
##' @param par A parameter list containing model parameters used during
##'   simulation.
##' @param x0 Numeric scalar giving the release x-coordinate.
##' @param y0 Numeric scalar giving the release y-coordinate.
##' @param t0 Numeric scalar giving the release time.
##' @param t1 Numeric scalar giving the final time of the simulated track.
##' @param dt Time step used for simulation.
##' @param id Optional tag identifier. If \code{NULL}, a random identifier is
##'   generated.
##' @param tag_type Character string giving the tag type. Used to select the
##'   appropriate observation-error component.
##' @param xygrid Matrix or data frame of active grid-cell centres. Required for
##'   CTMC simulation.
##' @param nextTo Neighbour structure of the grid, as produced for example by
##'   [get_neighbours()]. Required for CTMC simulation.
##' @param next_dist Numeric vector of neighbour distances. Required for CTMC
##'   simulation.
##' @param xgr Numeric vector of x-direction grid breaks.
##' @param ygr Numeric vector of y-direction grid breaks.
##' @param celltable Matrix linking grid-cell positions to cell indices.
##' @param xcen Numeric vector of x-coordinates of grid-cell centres.
##' @param ycen Numeric vector of y-coordinates of grid-cell centres.
##' @param sim_engine Integer selecting the simulation engine: \code{1} for
##'   continuous-space simulation and \code{2} for CTMC-based grid simulation.
##' @param use_reject Logical; if \code{TRUE}, rejected diffusion proposals are
##'   redrawn when simulated moves leave the valid domain.
##' @param n_reject Maximum number of rejection attempts if \code{use_reject =
##'   TRUE}.
##' @param ctmc_method Integer controlling the matrix-exponential method used in
##'   CTMC simulation.
##' @param add_obs_unc Add observation uncertainty to tag locations? By default
##'   (NULL), observation uncertainty is added to archival tags (tag_type =
##'   "d"), but not to other tags.
##'
##' @details
##' For \code{sim_engine = 1}, movement is simulated in continuous space by
##' combining taxis, advection, and diffusion increments at each time step.
##'
##' For \code{sim_engine = 2}, movement is simulated on the spatial grid by
##' constructing a transition-rate matrix from diffusion, taxis, and advection,
##' and then propagating the tag distribution forward over one time step.
##'
##' After simulation, observation error is added to all intermediate positions,
##' while the first and last positions are left unchanged.
##'
##' @return
##' A data frame with simulated tag positions over time. The returned data frame
##' contains at least the columns \code{t}, \code{x}, \code{y}, and \code{id}.
##'
##' @export
.sim_one_tag <- function(conf,
                        funcs,
                        par,
                        x0 = NULL,
                        y0 = NULL,
                        t0 = NULL,
                        t1 = NULL,
                        dt = 0.1,
                        id = NULL,
                        tag_type = "d",
                        xygrid = NULL,
                        nextTo = NULL,
                        next_dist = NULL,
                        xgr = NULL,
                        ygr = NULL,
                        celltable = NULL,
                        xcen = NULL,
                        ycen = NULL,
                        sim_engine = 1,
                        use_reject = FALSE,
                        n_reject = 20,
                        ctmc_method = 2,
                        add_obs_unc = NULL) {

  kappa <- exp(par$logKappa)
  sdO <- exp(par$logSdO)

  tag_type_int <- .get_tag_type_integer(tag_type)

  if (is.null(add_obs_unc)) {
    if (tag_type_int %in% c(1)) {
      add_obs_unc <- TRUE
    } else {
      add_obs_unc <- FALSE
    }
  }

  if (is.null(id)) id <- round(runif(1, 0, 1e4))

  xy <- matrix(c(x0, y0), 1, 2)

  ret <- data.frame(
    t = numeric(0),
    x = numeric(0),
    y = numeric(0),
    id = character(0),
    stringsAsFactors = FALSE
  )

  ret[1,] <- list(t = t0, x = x0, y = y0, id = as.character(id))

  t <- unname(t0)
  nc <- nrow(xygrid)

  flag_expm_uni <- ifelse(ctmc_method == 2, TRUE, FALSE)

  if (ctmc_method > 0) {
    mstar_template <- make_mstar_template(nextTo, ad = FALSE)
  }

  while ((t + dt) < t1) {

    if (sim_engine == 1) {

      moveT <- moveA <- moveD <- c(0,0)

      ## taxis
      if (conf$use_taxis) {
        moveT <- kappa * funcs$tax(xy, t) * dt
      }

      ## advection
      if (conf$use_advection) {
        moveA <- funcs$adv(xy, t) * dt
      }

      ## diffusion
      D <- exp(funcs$dif(xy, t))
      moveD <- rnorm(2, mean = 0,
                     sd = sqrt(2 * D * dt))

      ## Rejection method
      if (use_reject && !is.na(n_reject)) {
        cntr <- 0
        while (cntr < n_reject && (
          any(is.na(funcs$tax(xy + moveT + moveA + moveD, t) * dt)) ||
            any(is.na(funcs$adv(xy + moveT + moveA + moveD, t) * dt)) ||
            is.na(celltable[cbind(cut((xy + moveT + moveA + moveD)[1],
                                      xgr, include.lowest = TRUE),
                                  cut((xy + moveT + moveA + moveD)[2],
                                      ygr, include.lowest = TRUE))]))) {
                                        moveD <- rnorm(2, mean = 0,
                                                       sd = sqrt(2 * D * dt))
                                        cntr <- cntr + 1
                                      }
      }

      ## New position
      xy_new <- xy + moveT + moveA + moveD

    } else if (sim_engine == 2) {

      dist_prob <- numeric(nc)
      dist_prob[celltable[cbind(cut(xy[1], xgr, include.lowest = TRUE),
                                cut(xy[2], ygr, include.lowest = TRUE))]] <- 1

      ## Set to zero
      if (ctmc_method > 0) {
        Dstar <- Zstar <- Astar <- mstar_template
        Dstar@x[] <- Zstar@x[] <- Astar@x[] <- 0
      } else {
        Dstar <- Zstar <- Astar <- matrix(0, nc, nc)
      }

      ## diffusion
      D <- exp(funcs$dif(xygrid, t)) ## distance^2 / time
      ggrad <- funcs$ddif(xygrid, t)  # ∇g = ∇log D
      hD <- D * dt  ## (distance^2)
      for (k in 1:4) {
        j <- k + 1
        ind <- which(!is.na(nextTo[, j]))
        Dstar[cbind(ind, nextTo[ind, j])] <- hD[ind] / next_dist[k]^2
      }

      ## taxis
      if (conf$use_taxis) {
        move <- (D * (funcs$tax(xygrid, t) + ggrad)) * dt
        Zstar <- fill_inst_mat(Zstar, move, nextTo, next_dist)
      }

      ## advection
      if (conf$use_advection) {
        move <- funcs$adv(xygrid, t) * dt
        Astar <- fill_inst_mat(Astar, move, nextTo, next_dist)
      }

      ## movement rates
      Mstar <- Dstar + Zstar + Astar

      ## mass balance
      Mstar[cbind(1:nc, 1:nc)] <- 0
      Mstar[cbind(1:nc, 1:nc)] <- -RTMB::rowSums(Mstar)

      ## Check
      if (any(is.na(Mstar))) stop("NaN in Mstar!")

      if (ctmc_method > 0) {

        p <- as.vector(RTMB::expAv(Mstar,
                                   dist_prob,
                                   transpose = TRUE,
                                   uniformization = flag_expm_uni,
                                   rescale_freq = 1))

      } else {

        M <- Matrix::expm(Mstar)
        p <- as.vector(matrix(dist_prob, 1, nc) %*% M)

      }

      p <- guard_neg(p)
      p <- p / sum(p)

      ## New position
      ind_c <- sample.int(nc, 1, prob = p)
      ind_xy <- which(celltable == ind_c, arr.ind = TRUE)
      xy_new <- c(xcen[ind_xy[1]], ycen[ind_xy[2]])

    } else stop("Only sim_engine = 1 (Kalman filter) and 2 (CTMC) implemented!")

    xy <- xy_new
    t <- t + dt
    ret <- rbind(ret, list(t = t,
                           x = xy[1],
                           y = xy[2],
                           id = as.character(id)))
  }

  if (nrow(ret) >= 3 && isTRUE(add_obs_unc)) {
    ret[-c(1,nrow(ret)),"x"] <- as.numeric(ret[-c(1,nrow(ret)),"x"]) +
      rnorm(nrow(ret)-2, 0, sdO[1, tag_type_int])
    ret[-c(1,nrow(ret)),"y"] <- as.numeric(ret[-c(1,nrow(ret)),"y"]) +
      rnorm(nrow(ret)-2, 0, sdO[2, tag_type_int])
  }

  return(ret)
}




## s3 methods ----------------------------------------------------------------------

##' @rdname print-admove
##' @method print admove_tags
##' @export
print.admove_sim <- function(x, ...) {
  tmp <- x
  attributes(tmp) <- NULL
  NextMethod("print", tmp, ...)
}


##' @rdname plot_sim
##' @export
plot.admove_sim <- function(x, ...) {
  plot_sim(x, ...)
  return(invisible(NULL))
}
