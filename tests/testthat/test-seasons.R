
## require(admove); require(testthat)


make_dat <- function(ncov = 2, period = 12) {
  dat <- list(
    tags = NULL,
    cov = rep(list(array(1, dim = c(2, 2, 1))), ncov),
    time_spline = rep(list(0), ncov),
    period = period
  )
  names(dat$cov) <- paste0("cov", seq_len(ncov))
  dat
}


test_that("default_conf carries a non-seasonal n_seasons by default", {

  conf <- default_conf(make_dat(), verbose = FALSE)

  expect_equal(conf$n_seasons, c(1L, 1L))
  expect_equal(conf$seasonal_spline, c(FALSE, FALSE))
})


test_that("n_seasons drives seasonal_spline and the derived breakpoints", {

  dat <- make_dat()
  conf <- default_conf(dat, n_seasons = 4, verbose = FALSE)

  expect_equal(conf$n_seasons, c(4L, 4L))
  expect_equal(conf$seasonal_spline, c(TRUE, TRUE))

  res <- admove:::.resolve_seasons(dat, conf)

  expect_equal(as.numeric(res$dat$time_spline[[1]]), c(0, 3, 6, 9))
  expect_true(isTRUE(attr(res$dat$time_spline[[1]], "seasonal")))
})


test_that("seasonality can differ between covariates", {

  dat <- make_dat()
  conf <- default_conf(dat, n_seasons = c(2, 1), verbose = FALSE)
  res <- admove:::.resolve_seasons(dat, conf)

  expect_equal(as.numeric(res$dat$time_spline[[1]]), c(0, 6))
  expect_equal(as.numeric(res$dat$time_spline[[2]]), 0)
  expect_equal(res$conf$seasonal_spline, c(TRUE, FALSE))
})


test_that("resolving seasons leaves the data list alone and is idempotent", {

  dat <- make_dat()
  conf <- default_conf(dat, n_seasons = 3, verbose = FALSE)

  res <- admove:::.resolve_seasons(dat, conf)

  ## the caller's data list is untouched: seasonality lives in conf
  expect_equal(as.numeric(dat$time_spline[[1]]), 0)

  res2 <- admove:::.resolve_seasons(res$dat, res$conf)
  expect_equal(res2$dat$time_spline, res$dat$time_spline)
  expect_equal(res2$conf$n_seasons, res$conf$n_seasons)
})


test_that("parameter dimensions follow conf$n_seasons", {

  dat <- make_dat()
  conf <- default_conf(dat, n_seasons = 4, verbose = FALSE)
  par <- default_par(dat, conf, verbose = FALSE)

  expect_equal(dim(par$alpha)[3L], 4L)
  expect_equal(dim(par$beta)[3L], 4L)
  expect_equal(dim(par$gamma)[3L], 4L)

  map <- default_map(dat, conf, par)
  expect_length(map$alpha, length(par$alpha))
  expect_length(map$beta, length(par$beta))
})


test_that("a seasonal model without a period fails with an actionable error", {

  dat <- make_dat()
  dat$period <- NULL

  expect_error(default_conf(dat, n_seasons = 4, verbose = FALSE),
               "no seasonal cycle length is defined")
  expect_error(check_conf(list(n_seasons = 4), dat, verbose = FALSE),
               "no seasonal cycle length is defined")
})


test_that("n_seasons is validated", {

  dat <- make_dat()

  expect_error(default_conf(dat, n_seasons = 0, verbose = FALSE), "whole numbers")
  expect_error(default_conf(dat, n_seasons = 2.5, verbose = FALSE), "whole numbers")
  expect_error(default_conf(dat, n_seasons = c(2, 2, 2), verbose = FALSE),
               "one value per covariate")
})


test_that("set_seasons targets covariates by name and switches off again", {

  dat <- make_dat()
  conf <- default_conf(dat, verbose = FALSE)

  conf <- set_seasons(conf, dat, n = 4, cov = "cov2", verbose = FALSE)
  expect_equal(conf$n_seasons, c(1L, 4L))
  expect_equal(conf$seasonal_spline, c(FALSE, TRUE))

  conf <- set_seasons(conf, dat, n = 1, verbose = FALSE)
  expect_equal(conf$n_seasons, c(1L, 1L))
  expect_equal(conf$seasonal_spline, c(FALSE, FALSE))

  expect_error(set_seasons(conf, dat, n = 2, cov = "sst", verbose = FALSE),
               "Unknown covariate")
})


test_that("manual breakpoints still work and conf takes precedence", {

  dat <- make_dat()
  conf <- default_conf(dat, verbose = FALSE)

  ## legacy specification: breakpoints in dat, flag in conf
  dat$time_spline[[1]] <- c(0, 6)
  conf$seasonal_spline[1] <- TRUE

  res <- admove:::.resolve_seasons(dat, conf)
  expect_equal(as.numeric(res$dat$time_spline[[1]]), c(0, 6))
  ## marked as phases, so shifting the origin cannot rotate them
  expect_true(isTRUE(attr(res$dat$time_spline[[1]], "seasonal")))

  conf$n_seasons <- c(4L, 1L)
  expect_warning(res2 <- admove:::.resolve_seasons(dat, conf),
                 "configuration takes precedence")
  expect_equal(as.numeric(res2$dat$time_spline[[1]]), c(0, 3, 6, 9))
})


test_that("setting a period keeps the rest of the time reference intact", {

  dat <- skjepo$sim$dat
  tr0 <- tref(dat)

  period(dat) <- 6

  expect_equal(period(dat), 6)
  expect_s3_class(tref(dat), "admove_tref")
  expect_equal(origin(dat), tr0$origin)
  expect_equal(units_time(dat), tr0$units)

  ## and the spatial reference is untouched
  expect_equal(sref(dat), sref(skjepo$sim$dat))

  ## an object without a time reference is refused, not given a broken one
  expect_error(`period<-`(list(a = 1), 12), "tref")
})


test_that("shift_tref does not rotate seasonal breakpoints", {

  dat <- skjepo$sim$dat
  conf <- default_conf(dat, n_seasons = 4, verbose = FALSE)
  dat <- admove:::.resolve_seasons(dat, conf)$dat

  breaks0 <- as.numeric(dat$time_spline[[1]])
  trange0 <- dat$trange

  shifted <- suppressMessages(
    shift_tref(dat, origin = as.POSIXct("2000-06-15", tz = "UTC"),
               verbose = FALSE)
  )

  ## seasonal phases are anchored on the origin and stay put ...
  expect_equal(as.numeric(shifted$time_spline[[1]]), breaks0)
  expect_equal(as.numeric(shifted$time_spline[[1]])[1L], 0)

  ## ... while absolute times do move
  expect_false(isTRUE(all.equal(shifted$trange, trange0)))

  ## and the marker survives, so a second shift is safe too
  expect_true(isTRUE(attr(shifted$time_spline[[1]], "seasonal")))
})


test_that("scale_tref rescales seasonal breakpoints with the period", {

  dat <- skjepo$sim$dat
  conf <- default_conf(dat, n_seasons = 4, verbose = FALSE)
  dat <- admove:::.resolve_seasons(dat, conf)$dat

  breaks0 <- as.numeric(dat$time_spline[[1]])
  per0 <- period(dat)

  scaled <- scale_tref(dat, scale = 2, verbose = FALSE)

  expect_equal(as.numeric(scaled$time_spline[[1]]), breaks0 * 2)
  expect_equal(period(scaled), per0 * 2)
  expect_true(isTRUE(attr(scaled$time_spline[[1]], "seasonal")))
})
