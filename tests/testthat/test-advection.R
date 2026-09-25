## require(admove); require(testthat)


## A uniform field on the grid and time steps of `cov`, in model units.
uniform_adv <- function(cov, u, v) {
  fu <- cov
  fv <- cov
  fu[] <- u
  fv[] <- v
  prep_adv(fu, fv, units = NULL)
}

## Unit square with one simulated covariate and a uniform current.
adv_dat <- function(u = 0.4, v = 0.2, ...) {
  grid <- create_grid(cellsize = 0.25, verbose = FALSE)
  cov <- withr::with_seed(1, sim_cov(grid, nt = 2))
  suppressWarnings(suppressMessages(
    setup_data(grid = grid, cov = list(temp = cov),
               adv = list(cur = uniform_adv(cov, u, v)),
               trange = c(0, 1), verbose = FALSE, ...)
  ))
}

## A covariate field with a real CRS, for the unit conversion.
crs_cov <- function(xs, ys, val, crs, units) {
  a <- array(val, c(length(xs), length(ys), 2),
             dimnames = list(xs, ys, c(0, 1)))
  prep_cov(a, sref = list(crs = crs, units = units),
           tref = list(origin = as.Date("2020-01-01"), units = "month"),
           verbose = FALSE)
}


test_that("prep_adv checks its components", {

  grid <- create_grid(cellsize = 0.25, verbose = FALSE)
  cov <- withr::with_seed(1, sim_cov(grid, nt = 2))

  a <- uniform_adv(cov, 1, 2)
  expect_s3_class(a, "admove_adv")
  expect_null(attr(a, "units"))

  expect_error(prep_adv(cov, cov[, , 1, drop = FALSE]), "share the grid")
  expect_error(prep_adv(cov, cov, units = "furlongs"), "Unknown speed units")
  expect_equal(admove:::.adv_speed_factor("knots"), 1852 / 3600)
  expect_equal(admove:::.adv_speed_factor("cm s-1"), 0.01)
})


test_that("setup_data keeps advection fields apart from the covariates", {

  dat <- adv_dat()

  expect_length(dat$cov, 1L)
  expect_named(dat$adv, "cur")
  expect_equal(ncol(dat$knots_tax), 1L)
  expect_length(dat$time_adv, 2L)
  expect_equal(nrow(dat$xrange_adv), 2L)
  expect_true(all(unclass(dat$adv$cur$u) == 0.4))

  conf <- default_conf(dat, verbose = FALSE)
  expect_true(conf$use_advection)
  expect_equal(conf$adv_gamma, "shared")

  par <- default_par(dat, conf, verbose = FALSE)
  expect_equal(dim(par$gamma), c(2L, 1L, 1L))
  expect_equal(dim(par$adv_const), c(2L, 1L))
  expect_equal(dim(par$alpha)[2], 1L)
})


test_that("default_map sets up gamma and the constant drift from conf", {

  dat <- adv_dat()
  conf <- default_conf(dat, verbose = FALSE)
  par <- default_par(dat, conf, verbose = FALSE)

  map <- default_map(dat, conf, par)
  expect_equal(as.integer(map$gamma), c(1L, 1L))
  expect_true(all(is.na(map$adv_const)))

  conf$adv_gamma <- "xy"
  conf$adv_const <- TRUE
  map <- default_map(dat, conf, par)
  expect_equal(nlevels(map$gamma), 2L)
  expect_false(any(is.na(map$gamma)))
  expect_equal(nlevels(map$adv_const), 2L)

  conf$use_advection <- FALSE
  conf$adv_const <- FALSE
  map <- default_map(dat, conf, par)
  expect_true(all(is.na(map$gamma)))
  expect_true(all(is.na(map$adv_const)))
})


test_that(".make_adv combines fields, gamma, constant drift and seasons", {

  dat <- adv_dat(u = 0.4, v = 0.2)
  xy <- cbind(c(0.3, 0.6), c(0.4, 0.5))

  gamma <- array(c(0.5, 0.25), c(2, 1, 1))
  adv_const <- matrix(c(1, -1), 2, 1)
  a <- admove:::.make_adv(dat$adv, dat$time_adv, gamma, adv_const)
  expect_equal(unname(a$val(xy, 0.5)),
               cbind(rep(0.5 * 0.4 + 1, 2), rep(0.25 * 0.2 - 1, 2)))

  ## a single position given as a vector (as in the CTMC simulator)
  expect_equal(unname(a$val(c(0.3, 0.4), 0.5)), cbind(1.2, -0.95))

  ## constant drift only, two seasons over a period of 1
  a2 <- admove:::.make_adv(NULL, NULL, NULL, cbind(c(1, 0), c(0, 2)),
                           period = 1)
  expect_equal(unname(a2$val(xy, 0.25)), cbind(c(1, 1), c(0, 0)))
  expect_equal(unname(a2$val(xy, 0.75)), cbind(c(0, 0), c(2, 2)))
})


test_that("the conversion from m/s follows the geometry of the grid", {

  skip_if_not_installed("sf")

  ## lon/lat: an east speed at 60N covers 1 / cos(60) as many degrees,
  ## with the WGS84 radius of curvature in the prime vertical
  xs <- 0:10
  ys <- 55:65
  a <- prep_adv(crs_cov(xs, ys, 1, 4326, "degree"),
                crs_cov(xs, ys, 0, 4326, "degree"), units = "m/s")
  a <- suppressMessages(admove:::.adv_convert_field(a, "m/s"))
  s_month <- 365.25 / 12 * 86400
  n60 <- 6378137 / sqrt(1 - 0.00669437999014 * sin(pi / 3)^2)
  expect_equal(unclass(a$u)[1, which(ys == 60), 1],
               s_month * 180 / (pi * n60 * cos(pi / 3)))
  expect_equal(unclass(a$v)[1, which(ys == 60), 1], 0)

  ## azimuthal equidistant in km: a north current keeps its speed along the
  ## central meridian and turns away from grid north off it
  crs <- "+proj=aeqd +lat_0=0 +lon_0=60 +datum=WGS84 +units=m +no_defs"
  xs <- seq(-3000, 3000, 1000)
  ys <- seq(-3000, 3000, 1000)
  b <- prep_adv(crs_cov(xs, ys, 0, crs, "km"), crs_cov(xs, ys, 1, crs, "km"),
                units = "m/s")
  b <- suppressMessages(admove:::.adv_convert_field(b, "m/s"))
  i0 <- which(xs == 0)
  j <- which(ys == 2000)
  expect_equal(unclass(b$u)[i0, j, 1], 0, tolerance = 1e-6)
  expect_equal(unclass(b$v)[i0, j, 1], s_month / 1000, tolerance = 1e-6)
  expect_lt(unclass(b$u)[length(xs), j, 1], -100)
  expect_gt(unclass(b$u)[1, j, 1], 100)

  ## setup_data() converts fields with units and leaves the others alone
  u <- crs_cov(xs, ys, 1, crs, "km")
  v <- crs_cov(xs, ys, 0, crs, "km")
  dat <- suppressWarnings(suppressMessages(setup_data(
    cov = list(temp = u), adv = list(a = prep_adv(u, v, units = "m/s"),
                                     b = prep_adv(u, v, units = NULL)),
    verbose = FALSE)))
  expect_gt(unclass(dat$adv$a$u)[i0, j, 1], 2000)
  expect_equal(unclass(dat$adv$b$u)[i0, j, 1], 1)
})


test_that("advection with zero coefficients leaves the likelihood unchanged", {

  sim <- withr::with_seed(2, suppressWarnings(suppressMessages(
    sim_data(cov = adv_dat()$cov, grid = adv_dat()$grid,
             adv = adv_dat()$adv, trange = c(0, 1),
             par = list(gamma = array(0, c(2, 1, 1))),
             n_dtags = 1, n_ctags = 10, verbose = FALSE))))

  conf1 <- sim$conf
  conf0 <- conf1
  conf0$use_advection <- FALSE
  par <- sim$par

  f0 <- suppressWarnings(admove(sim$dat, conf0, par, run = FALSE,
                                verbose = FALSE))
  f1 <- suppressWarnings(admove(sim$dat, conf1, par, run = FALSE,
                                verbose = FALSE))
  p1 <- f1$obj$par
  expect_equal(unname(p1[names(p1) == "gamma"]), 0)
  expect_equal(f1$obj$fn(p1), f0$obj$fn(f0$obj$par))

  ## and a non-zero gamma does change it
  p1[names(p1) == "gamma"] <- 0.5
  expect_false(isTRUE(all.equal(f1$obj$fn(p1), f0$obj$fn(f0$obj$par))))
})


test_that("the entrainment coefficient is recovered from simulated tags", {

  grid <- create_grid(cellsize = 0.1, verbose = FALSE)
  cov <- withr::with_seed(3, sim_cov(grid, nt = 2))
  sim <- withr::with_seed(3, suppressWarnings(suppressMessages(
    sim_data(cov = cov, grid = grid, adv = uniform_adv(cov, 0.4, 0.2),
             trange = c(0, 1), par = list(gamma = array(0.5, c(2, 1, 1))),
             xrange_rel = c(0.1, 0.4), yrange_rel = c(0.1, 0.4),
             n_dtags = 3, n_ctags = 80, verbose = FALSE))))

  ## Taxis and advection both explain directed movement, and from alpha = 0,
  ## gamma = 0 the optimizer can end in a local optimum where taxis takes the
  ## drift. Starting gamma from a fit with the taxis fixed avoids that.
  map1 <- sim$map
  map1$alpha <- factor(rep(NA, length(sim$par$alpha)))
  fit1 <- suppressWarnings(admove(sim$dat, sim$conf, sim$par, map1,
                                  do_predictions = FALSE, do_report = FALSE,
                                  do_sdreport = FALSE, verbose = FALSE))
  fit <- suppressWarnings(admove(sim$dat, sim$conf, fit1$pl, sim$map,
                                 do_predictions = FALSE, do_report = FALSE,
                                 verbose = FALSE))

  expect_equal(fit$opt$convergence, 0L)
  g <- fit$pl$gamma[1, 1, 1]
  se <- fit$plsd$gamma[1, 1, 1]
  expect_lt(abs(g - 0.5), 3 * se)
  expect_equal(fit$pl$gamma[1, 1, 1], fit$pl$gamma[2, 1, 1])
})


test_that("plot_adv_field draws the input fields", {

  dat <- adv_dat()
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())

  expect_null(plot_adv_field(dat$adv$cur))
  expect_null(plot_adv_field(dat$adv, i = "cur", select = 1:2))
  expect_null(plot_adv_field(dat, select = 1))
  expect_null(plot(dat$adv$cur, select = 2))
  expect_null(plot_data(dat))

  dat0 <- dat
  dat0$adv <- NULL
  expect_error(plot_adv_field(dat0), "no advection field")
  expect_error(plot_adv_field(dat, i = "wind"), "No advection field 'wind'")
  expect_error(plot_adv_field(dat, select = 5), "select")
})
