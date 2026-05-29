

test_that("nll builds and evaluates for KF engine", {
  skip_if_not_installed("RTMB")

  dat <- sim_data()$dat
  conf <- default_conf(dat)
  conf$engine <- 1

  par <- default_par(dat, conf)
  map <- default_map(dat, conf, par)

  ## mimic what admove() does before calling nll()
  dat$tags <- check_tags(dat$tags, dat$grid, dat, conf, TRUE, verbose = FALSE)
  dat$period <- period(dat)

  tmb_dat <- c(dat, conf)
  tmb_dat$tags <- split(dat$tags, dat$tags$id)
  tmb_dat$dbg <- FALSE

  obj <- RTMB::MakeADFun(
    func = function(p) admove:::nll(p, tmb_dat),
    parameters = par,
    map = map,
    silent = TRUE
  )

  val <- obj$fn(obj$par)

  expect_type(val, "double")
  expect_length(val, 1)
  expect_true(is.finite(val))
})


test_that("nll builds and evaluates for CTMC engine when ic is available", {
  skip_if_not_installed("RTMB")

  dat <- sim_data()$dat
  conf <- default_conf(dat)
  conf$engine <- 2

  par <- default_par(dat, conf)
  map <- default_map(dat, conf, par)

  ## ensure tags are matched to grid cells
  dat$tags <- check_tags(dat$tags, dat$grid, dat, conf, TRUE, verbose = FALSE)
  expect_true("ic" %in% names(dat$tags))

  dat$period <- period(dat)

  tmb_dat <- c(dat, conf)
  tmb_dat$tags <- split(dat$tags, dat$tags$id)
  tmb_dat$dbg <- FALSE

  obj <- RTMB::MakeADFun(
    func = function(p) admove:::nll(p, tmb_dat),
    parameters = par,
    map = map,
    silent = TRUE
  )

  val <- obj$fn(obj$par)

  expect_type(val, "double")
  expect_length(val, 1)
  expect_true(is.finite(val))
})


test_that("nll errors for unsupported engine", {
  skip_if_not_installed("RTMB")

  dat <- sim_data()$dat
  conf <- default_conf(dat)
  conf$engine <- 99

  par <- default_par(dat, conf)
  map <- default_map(dat, conf, par)

  dat$tags <- check_tags(dat$tags, dat$grid, dat, conf, TRUE, verbose = FALSE)
  dat$period <- period(dat)

  tmb_dat <- c(dat, conf)
  tmb_dat$tags <- split(dat$tags, dat$tags$id)
  tmb_dat$dbg <- FALSE

  expect_error(
  obj <- RTMB::MakeADFun(
    func = function(p) admove:::nll(p, tmb_dat),
    parameters = par,
    map = map,
    silent = TRUE
  ),
  "This engine is not yet implemented"
  )
})


test_that("nll responds to parameter changes", {
  skip_if_not_installed("RTMB")

  dat <- sim_data()$dat
  conf <- default_conf(dat)
  conf$engine <- 1

  par <- default_par(dat, conf)
  map <- default_map(dat, conf, par)

  dat$tags <- check_tags(dat$tags, dat$grid, dat, conf, TRUE, verbose = FALSE)
  dat$period <- period(dat)

  tmb_dat <- c(dat, conf)
  tmb_dat$tags <- split(dat$tags, dat$tags$id)
  tmb_dat$dbg <- FALSE

  obj <- RTMB::MakeADFun(
    func = function(p) admove:::nll(p, tmb_dat),
    parameters = par,
    map = map,
    silent = TRUE
  )

  val1 <- obj$fn(obj$par)

  par2 <- obj$par
  par2[1] <- par2[1] + 0.1

  val2 <- obj$fn(par2)

  expect_true(is.finite(val1))
  expect_true(is.finite(val2))
  expect_false(isTRUE(all.equal(val1, val2)))
})
