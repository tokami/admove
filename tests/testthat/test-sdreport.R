## add_sdreport(): the AD Hessian, and the '...' passthrough to RTMB::sdreport()


test_that("add_sdreport passes ... through to RTMB::sdreport", {

  fit <- small_fit()

  ## skip.delta.method is the clearest evidence that ... arrives: it leaves the
  ## ADREPORTed values in place but drops their standard errors
  a <- suppressWarnings(add_sdreport(fit))
  b <- suppressWarnings(add_sdreport(fit, skip.delta.method = TRUE))

  expect_false(all(is.na(a$sdrep$sd)))
  expect_true(all(is.na(b$sdrep$sd)))
  expect_equal(length(b$sdrep$value), length(a$sdrep$value))

  ## and the fixed effects are untouched by it
  expect_equal(b$sdrep$cov.fixed, a$sdrep$cov.fixed)
  expect_equal(b$plsd$beta, a$plsd$beta)
})


test_that("add_sdreport keeps its own arguments out of ...", {

  fit <- small_fit()

  ## save_covariance and ad_hessian are matched before ..., so they are still
  ## the local arguments and never reach RTMB::sdreport()
  expect_false("cov" %in% names(suppressWarnings(add_sdreport(fit))$sdrep))
  expect_true("cov" %in%
                names(suppressWarnings(add_sdreport(fit, save_covariance = TRUE))$sdrep))

  ## ad_hessian = FALSE drops the obj$he() call (the ADHess tape, which is where
  ## the memory goes on a large model -- see ?add_sdreport). The estimates are
  ## the same; only the route to the Hessian differs.
  s_ad <- suppressWarnings(add_sdreport(fit, ad_hessian = TRUE))
  s_fd <- suppressWarnings(add_sdreport(fit, ad_hessian = FALSE))
  expect_equal(s_fd$pl, s_ad$pl)
  expect_equal(dim(s_fd$sdrep$cov.fixed), dim(s_ad$sdrep$cov.fixed))
})


test_that("the AD Hessian differs from TMB's default differencing, and is the exact one", {

  ## This is the reason ad_hessian = TRUE exists. RTMB::sdreport() builds the
  ## Hessian with optimHess(par, obj$fn, obj$gr) at a *fixed absolute* step
  ## (ndeps, default 1e-3), which is far too coarse when the parameters sit on
  ## a small scale. On small_fit() that route is ~59% off the exact Hessian in
  ## its largest entry; the finite differences only converge to the AD value as
  ## the step shrinks, which is what identifies the AD one as correct.
  fit <- small_fit()
  p <- fit$opt$par

  h_ad <- admove:::.get_ad_hessian(fit$obj)
  expect_true(is.matrix(h_ad))
  expect_equal(dim(h_ad), c(length(p), length(p)))

  rel <- function(h) max(abs(h - h_ad) / pmax(abs(h_ad), 1e-8))

  ## sdreport's own route, i.e. what ad_hessian = FALSE falls back to
  h_sdrep <- stats::optimHess(p, fit$obj$fn, fit$obj$gr,
                              control = list(ndeps = rep(1e-3, length(p))))
  expect_gt(rel(h_sdrep), 0.1)

  ## differencing the objective converges to the AD Hessian as the step shrinks
  step <- function(nd) rel(stats::optimHess(p, fit$obj$fn,
                                            control = list(ndeps = rep(nd, length(p)))))
  expect_gt(step(1e-2), step(1e-4))
  expect_lt(step(1e-4), 0.1)
})
