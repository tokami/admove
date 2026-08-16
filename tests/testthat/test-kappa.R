## require(admove); require(testthat)


test_that("default_map fixes logKappa", {

  map <- with(skjepo$sim, default_map(dat, conf, par))

  expect_true(all(is.na(map$logKappa)))
})


test_that("the default kappa is the geometric mean of the per-tag-type values", {

  dat <- skjepo$sim$dat
  conf <- skjepo$sim$conf

  kt <- admove:::.kappa_by_tag_type(dat,
                                    admove:::.get_tags_in_use(dat, conf),
                                    seq_along(dat$cov))

  ## the example carries both data-storage and mark-recapture tags
  expect_gt(length(kt), 1)

  par <- default_par(dat, conf, verbose = FALSE)

  expect_equal(exp(par$logKappa), exp(mean(log(kt))))

  ## a pooled median could fall outside the per-type values; a geometric mean
  ## cannot
  expect_gt(exp(par$logKappa), min(kt))
  expect_lt(exp(par$logKappa), max(kt))
})


test_that("cov_taxis restricts which covariates scale kappa", {

  dat <- skjepo$sim$dat
  conf <- skjepo$sim$conf

  ## add a second covariate with a different range-to-gradient ratio, as an
  ## advection input would have, and check that it only affects kappa when it is
  ## included. (Rescaling a covariate would not do: kappa = L * R / (G * T) is
  ## invariant to that, since R and G scale together.)
  dat2 <- dat
  co <- dat$cov[[1]]
  co[] <- array(rep(seq_len(dim(co)[1]), prod(dim(co)[-1])), dim(co))
  dat2$cov <- c(dat$cov, list(co))
  dat2$time_cov <- c(dat$time_cov, dat$time_cov[1])

  k_all <- admove:::.kappa_by_tag_type(dat2,
                                       admove:::.get_tags_in_use(dat2, conf),
                                       seq_along(dat2$cov))
  k_one <- admove:::.kappa_by_tag_type(dat2,
                                       admove:::.get_tags_in_use(dat2, conf),
                                       1L)
  k_ref <- admove:::.kappa_by_tag_type(dat,
                                       admove:::.get_tags_in_use(dat, conf),
                                       1L)

  expect_false(isTRUE(all.equal(unname(k_all), unname(k_one))))
  expect_equal(k_one, k_ref)

  ## invariance to a pure rescaling of a covariate
  dat3 <- dat
  dat3$cov[[1]][] <- unclass(dat$cov[[1]]) * 100
  expect_equal(admove:::.kappa_by_tag_type(dat3, admove:::.get_tags_in_use(dat3, conf), 1L),
               k_ref)
})


test_that("cov_taxis rejects covariates that do not exist", {

  dat <- skjepo$sim$dat

  expect_error(default_par(dat, cov_taxis = 99, verbose = FALSE),
               "must index covariates")
  expect_error(default_par(dat, cov_taxis = "nope", verbose = FALSE),
               "Unknown covariate")
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
