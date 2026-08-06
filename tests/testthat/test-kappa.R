## require(admove); require(testthat)


test_that("default_map fixes logKappa", {

  map <- with(skjepo$sim, default_map(dat, conf, par))

  expect_true(all(is.na(map$logKappa)))
})


test_that(".check_kappa_map is silent when logKappa is fixed", {

  sim <- skjepo$sim

  expect_silent(
    admove:::.check_kappa_map(sim$par, sim$map, sim$conf)
  )
})


test_that(".check_kappa_map warns when logKappa is estimated", {

  sim <- skjepo$sim
  map <- sim$map
  map$logKappa <- factor(1)

  expect_warning(
    admove:::.check_kappa_map(sim$par, map, sim$conf),
    "confounded with alpha"
  )

  ## logKappa missing from map is also estimated by RTMB
  map$logKappa <- NULL

  expect_warning(
    admove:::.check_kappa_map(sim$par, map, sim$conf),
    "confounded with alpha"
  )
})


test_that(".check_kappa_map is silent when alpha is anchored at a non-zero value", {

  sim <- skjepo$sim
  map <- sim$map
  map$logKappa <- factor(1)

  ## alpha[1,,] is fixed by default_map(); holding it at a non-zero value pins
  ## the scale of alpha and makes kappa identifiable
  par <- sim$par
  par$alpha[1, 1, 1] <- 0.5

  expect_silent(
    admove:::.check_kappa_map(par, map, sim$conf)
  )
})


test_that(".check_kappa_map is silent when taxis is switched off", {

  sim <- skjepo$sim
  map <- sim$map
  map$logKappa <- factor(1)
  conf <- sim$conf
  conf$use_taxis <- FALSE

  expect_silent(
    admove:::.check_kappa_map(sim$par, map, conf)
  )
})


test_that("summary reports the fixed kappa scale", {

  ## nlminb may warn about NA/NaN function evaluations on this example
  fit <- suppressWarnings(
    admove(skjepo$sim, do_sdreport = FALSE, do_predictions = FALSE,
           do_report = FALSE, verbose = FALSE)
  )

  out <- capture.output(summary(fit))

  expect_true(any(grepl("kappa (fixed scale)", out, fixed = TRUE)))
})
