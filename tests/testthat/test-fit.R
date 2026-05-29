
## require(admove); require(testthat)


test_that("admove() returns an admove object for a minimal fit", {
  skip_if_not_installed("RTMB")

  data("skjepo", package = "admove")

  fit <- admove(
    skjepo$sim,
    do_predictions = FALSE,
    do_sdreport = FALSE,
    do_report = FALSE,
    verbose = FALSE
  )

  expect_s3_class(fit, "admove")
  expect_type(fit, "list")
  expect_true(all(c("dat", "conf", "par", "map", "opt", "obj") %in% names(fit)))
  expect_s3_class(fit$dat, "admove_data")
})


test_that("admove() engine argument overrides conf$engine", {
  skip_if_not_installed("RTMB")

  data("skjepo", package = "admove")

  conf <- default_conf(skjepo$sim)
  conf$engine <- 2

  fit <- admove(
    skjepo$sim$dat,
    conf = conf,
    engine = 1,
    do_predictions = FALSE,
    do_sdreport = FALSE,
    do_report = FALSE,
    verbose = FALSE
  )

  expect_identical(admove:::.get_engine_integer(fit$conf$engine), 1L)

})


test_that("admove() without sdreport stores pl when do_sdreport = FALSE", {
  skip_if_not_installed("RTMB")

  data("skjepo", package = "admove")

  fit <- admove(
    skjepo$sim,
    do_predictions = FALSE,
    do_sdreport = FALSE,
    do_report = FALSE,
    verbose = FALSE
  )

  expect_true("pl" %in% names(fit))
  expect_false("sdrep" %in% names(fit))
  expect_true(is.list(fit$pl) || is.vector(fit$pl))
})


test_that("admove() errors when conf$obs_var_type and map$logSdO do not agree", {
  skip_if_not_installed("RTMB")

  data("skjepo", package = "admove")

  conf <- default_conf(skjepo$sim$dat)
  par  <- default_par(skjepo$sim$dat, conf)
  map  <- default_map(skjepo$sim$dat, conf, par)

  ## Force a mismatch: fix all observation-SD parameters in the map
  ## but keep conf$obs_var_type as estimated
  map$logSdO[] <- factor(NA)
  conf$obs_var_type[] <- TRUE

  expect_error(
    admove(
      skjepo$sim$dat,
      conf = conf,
      par = par,
      map = map,
      do_predictions = FALSE,
      do_sdreport = FALSE,
      do_report = FALSE,
      verbose = FALSE
    ),
    "conf\\$obs_var_type and mapped parameters \\(map\\) do not agree"
  )

})


test_that("admove() returns objective without optimisation when run = FALSE", {
  skip_if_not_installed("RTMB")

  dat <- sim_data()

  obj <- admove(
    dat,
    run = FALSE,
    verbose = FALSE
  )

  expect_type(obj, "list")
  expect_true("obj" %in% names(obj))
  expect_true(is.function(obj$obj$fn))
})
