## require(admove); require(testthat)

## nlminb can report convergence = 0 ("relative convergence") at a point that is
## not stationary, where sdreport() still returns a positive-definite Hessian and
## small standard errors. The fit therefore judges on the gradient as well.


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


test_that(".not_converged flags a large gradient even when nlminb is happy", {

  ok <- list(convergence = 0L, message = "relative convergence (4)")

  expect_false(admove:::.not_converged(ok, 1e-4, 1e-2))
  expect_true(admove:::.not_converged(ok, 60, 1e-2))
  expect_true(admove:::.not_converged(ok, NA_real_, 1e-2))
  expect_true(admove:::.not_converged(list(convergence = 1L), 1e-8, 1e-2))
})


test_that("a fit stores the maximum gradient and summary reports it", {

  fit <- suppressWarnings(
    admove(skjepo$sim, do_sdreport = FALSE, do_predictions = FALSE,
           do_report = FALSE, verbose = FALSE)
  )

  expect_true(is.numeric(fit$max_gradient))
  expect_equal(fit$grad_tol, 1e-4)

  out <- capture.output(summary(fit))
  expect_true(any(grepl("Max. gradient component", out, fixed = TRUE)))
})


test_that("summary warns when the optimizer stopped away from a stationary point", {

  fit <- suppressWarnings(
    admove(skjepo$sim, do_sdreport = FALSE, do_predictions = FALSE,
           do_report = FALSE, verbose = FALSE)
  )

  ## pretend the optimizer stopped with a large gradient but reported success
  fit$opt$convergence <- 0L
  fit$max_gradient <- 60

  out <- capture.output(summary(fit))
  expect_true(any(grepl("did not obtain proper convergence", out)))
  expect_true(any(grepl("gradient is not zero", out)))
})


test_that("n_restarts = 0 skips the restart loop", {

  ## only checks that the argument is accepted and the fit is unchanged in shape
  fit <- suppressWarnings(
    admove(skjepo$sim, n_restarts = 0, do_sdreport = FALSE,
           do_predictions = FALSE, do_report = FALSE, verbose = FALSE)
  )

  expect_true(is.numeric(fit$max_gradient))
})


test_that("the gradient tolerance scales with the objective", {

  ## the same absolute gradient is fine for a large objective and not for a
  ## small one: the gradient scales with the number of observations
  big   <- list(convergence = 0L, objective = 11244)
  small <- list(convergence = 0L, objective = 12)

  expect_equal(admove:::.grad_threshold(1e-4, 11244), 1.1244)
  expect_equal(admove:::.grad_threshold(1e-4, 12), 12e-4)

  ## the values the EPO skipjack fits actually stop at
  for (g in c(0.013, 0.034, 0.087)) expect_false(admove:::.not_converged(big, g, 1e-4))
  expect_true(admove:::.not_converged(big, 60, 1e-4))
  expect_true(admove:::.not_converged(small, 0.087, 1e-4))

  ## objectives below 1 do not shrink the threshold further
  expect_equal(admove:::.grad_threshold(1e-4, 0.001), 1e-4)
  expect_equal(admove:::.grad_threshold(1e-4, NULL), 1e-4)
})


test_that("a zero gradient outranks the optimizer's status code", {

  big <- list(convergence = 1L, message = "false convergence (8)", objective = 15050)

  ## nlminb unhappy but the gradient is at zero: a stalled line search AT the
  ## optimum, not a failed fit
  expect_equal(admove:::.convergence_status(big, 0.017, 1e-4), "gradient_ok")

  ## nlminb happy but the gradient is not: the dangerous case
  expect_equal(admove:::.convergence_status(
    list(convergence = 0L, objective = 11468), 60, 1e-4), "bad")

  expect_equal(admove:::.convergence_status(
    list(convergence = 0L, objective = 11244), 0.087, 1e-4), "ok")

  ## a restart is still attempted for anything short of "ok"
  expect_true(admove:::.not_converged(big, 0.017, 1e-4))
  expect_false(admove:::.not_converged(
    list(convergence = 0L, objective = 11244), 0.087, 1e-4))
})


test_that("summary notes, rather than condemns, a stalled line search", {

  fit <- suppressWarnings(
    admove(skjepo$sim, do_sdreport = FALSE, do_predictions = FALSE,
           do_report = FALSE, verbose = FALSE)
  )
  fit$opt$convergence <- 1L
  fit$opt$message <- "false convergence (8)"
  fit$max_gradient <- 1e-8

  out <- capture.output(summary(fit))
  expect_false(any(grepl("did not obtain proper convergence", out)))
  expect_true(any(grepl("stalled line search", out)))
})
