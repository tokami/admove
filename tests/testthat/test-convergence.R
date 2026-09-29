## require(admove); require(testthat)

## A fit passes only if every convergence check passes: optimizer status,
## finite objective, gradient, bounds, positive-definite Hessian, finite SEs.


test_that(".max_abs_gradient returns the largest absolute gradient component", {

  obj <- list(gr = function(p) c(-3, 1, 2))

  expect_equal(admove:::.max_abs_gradient(obj, c(0, 0, 0)), 3)
})


test_that(".max_abs_gradient is NA when the gradient cannot be evaluated", {

  expect_true(is.na(admove:::.max_abs_gradient(NULL, 1)))
  expect_true(is.na(admove:::.max_abs_gradient(list(gr = function(p) stop("no")), 1)))

  ## nothing to differentiate: everything fixed in the map
  expect_equal(admove:::.max_abs_gradient(list(gr = function(p) numeric(0)),
                                          numeric(0)), 0)
})


## minimal fit that passes every check; tests break one thing at a time
.good_fit <- function() {
  cov <- matrix(c(1, 0.2, 0.2, 1), 2, dimnames = list(c("alpha", "alpha"),
                                                      c("alpha", "alpha")))
  list(opt = list(convergence = 0L, message = "relative convergence (4)",
                  objective = 100, par = c(alpha = 0.1, alpha = 0.2)),
       max_gradient = 1e-6,
       grad_tol = 1e-3,
       at_bound = character(0),
       sdrep = structure(list(pdHess = TRUE, cov.fixed = cov), class = "sdreport"))
}


test_that("a fit passes when every check passes", {

  checks <- admove:::.convergence_checks(.good_fit())

  expect_true(all(checks))
  expect_equal(names(checks),
               c("optimizer", "objective", "gradient", "bounds", "hessian", "se"))
  expect_equal(admove:::.convergence_verdict(checks), "pass")
})


test_that("each check fails on its own", {

  break_fit <- list(
    optimizer = function(f) { f$opt$convergence <- 1L; f },
    objective = function(f) { f$opt$objective <- NaN; f },
    gradient = function(f) { f$max_gradient <- 60; f },
    bounds = function(f) { f$at_bound <- "alpha1"; f },
    hessian = function(f) { f$sdrep$pdHess <- FALSE; f },
    se = function(f) { f$sdrep$cov.fixed[1, 1] <- -1; f }
  )

  for (nm in names(break_fit)) {
    checks <- admove:::.convergence_checks(break_fit[[nm]](.good_fit()))
    expect_equal(names(checks)[!checks], nm, info = nm)
    expect_equal(admove:::.convergence_verdict(checks), "fail", info = nm)
  }
})


test_that("false convergence fails even with a zero gradient", {

  fit <- .good_fit()
  fit$opt$convergence <- 1L
  fit$opt$message <- "false convergence (8)"
  fit$max_gradient <- 1e-8

  checks <- admove:::.convergence_checks(fit)
  expect_equal(admove:::.convergence_verdict(checks), "fail")
  expect_match(admove:::.convergence_message(fit, checks),
               "did not pass convergence checks: optimizer: false convergence (8).",
               fixed = TRUE)
})


test_that("the gradient tolerance is absolute", {

  fit <- .good_fit()
  fit$opt$objective <- 1e6
  fit$max_gradient <- 0.01

  expect_false(admove:::.convergence_checks(fit)[["gradient"]])
  fit$grad_tol <- 0.1
  expect_true(admove:::.convergence_checks(fit)[["gradient"]])
})


test_that("without an sdreport the Hessian and SEs are not checked", {

  fit <- .good_fit()
  fit$sdrep <- NULL

  checks <- admove:::.convergence_checks(fit)
  expect_true(is.na(checks[["hessian"]]))
  expect_true(is.na(checks[["se"]]))
  expect_equal(admove:::.convergence_verdict(checks), "partial")
  expect_match(admove:::.convergence_message(fit, checks), "not checked: hessian, se",
               fixed = TRUE)
})


test_that(".at_bound finds estimates on a finite bound", {

  par <- c(a = 1, a = 2, b = 0)

  expect_equal(admove:::.at_bound(par, rep(-Inf, 3), rep(Inf, 3)), character(0))
  expect_equal(admove:::.at_bound(par, c(-Inf, -Inf, 0), rep(Inf, 3)), "b1")
  expect_equal(admove:::.at_bound(par, rep(-Inf, 3), c(5, 2, Inf)), "a2")
  ## bounds that do not line up with the estimated parameters: not checked
  expect_null(admove:::.at_bound(par, c(0, 0), c(Inf, Inf)))
})


test_that("a fit stores its convergence checks and summary lists them", {

  fit <- small_fit()

  expect_true(is.numeric(fit$max_gradient))
  expect_equal(fit$grad_tol, 1e-3)
  expect_true(is.logical(fit$convergence))

  out <- capture.output(summary(fit))
  expect_true(any(grepl("Convergence checks", out, fixed = TRUE)))
  expect_true(any(grepl("max|gradient|", out, fixed = TRUE)))
})


test_that("the small fit passes all checks with an sdreport", {

  fit <- small_fit(sdreport = TRUE)

  expect_equal(admove:::.convergence_verdict(fit$convergence), "pass")
})


test_that("summary warns when a check fails", {

  fit <- small_fit()
  fit$opt$convergence <- 1L
  fit$opt$message <- "false convergence (8)"
  fit$max_gradient <- 1e-8

  out <- capture.output(summary(fit))
  expect_true(any(grepl("[FAIL] optimizer: false convergence (8)", out, fixed = TRUE)))
  expect_true(any(grepl("did not pass convergence checks", out)))
})


test_that("Newton steps close the gradient gap left by nlminb", {

  ## quadratic with an off-diagonal Hessian: one Newton step is exact
  target <- c(1, -2, 3)
  A <- matrix(c(4, 1, 0, 1, 3, 1, 0, 1, 2), 3)
  obj <- RTMB::MakeADFun(function(p) {
    d <- p$a - target
    0.5 * sum(d * as.vector(A %*% d))
  }, list(a = c(0, 0, 0)), silent = TRUE)
  start <- target + c(1e-3, -2e-3, 1e-3)
  opt <- list(par = setNames(start, rep("a", 3)), objective = obj$fn(start))

  nt <- admove:::.newton_steps(obj, opt, 2, rep(-Inf, 3), rep(Inf, 3))
  expect_equal(unname(nt$opt$par), target, tolerance = 1e-10)
  expect_lt(nt$steps$max_gradient[nrow(nt$steps)], 1e-10)
  expect_equal(unname(obj$env$last.par.best), unname(nt$opt$par))

  ## the same with the finite-difference Hessian
  nt2 <- admove:::.newton_steps(obj, opt, 1, rep(-Inf, 3), rep(Inf, 3),
                                ad_hessian = FALSE)
  expect_equal(unname(nt2$opt$par), target, tolerance = 1e-6)

  ## a step that would leave the bounds is not taken
  nt3 <- admove:::.newton_steps(obj, opt, 1,
                                c(-Inf, -Inf, target[3] + 5e-4), rep(Inf, 3))
  expect_identical(nt3$opt$par, opt$par)
  expect_equal(nrow(nt3$steps), 1L)
  expect_equal(unname(obj$env$last.par.best), start)

  ## switched off
  nt4 <- admove:::.newton_steps(obj, opt, 0, rep(-Inf, 3), rep(Inf, 3))
  expect_identical(nt4$opt, opt)
})


test_that("admove() records the Newton steps and keeps last.par.best on the estimates", {

  fit <- small_fit()
  expect_s3_class(fit$newton, "data.frame")
  expect_equal(fit$max_gradient,
               fit$newton$max_gradient[nrow(fit$newton)], tolerance = 1e-8)
  expect_equal(unname(fit$obj$env$last.par.best), unname(fit$opt$par))
  expect_equal(fit$opt$objective, fit$newton$objective[nrow(fit$newton)])
})


test_that("admove() passes 'control' on to nlminb", {

  sim <- tiny_sim()
  fit <- suppressWarnings(
    admove(sim, control = list(iter.max = 1), newton_steps = 0,
           do_sdreport = FALSE, do_predictions = FALSE, do_report = FALSE,
           verbose = FALSE)
  )
  expect_lte(fit$opt$iterations, 1)
  expect_match(fit$opt$message, "iteration limit")
  expect_equal(nrow(fit$newton), 1L)

  expect_error(suppressWarnings(admove(sim, control = list(1), verbose = FALSE)),
               "named list")
})
