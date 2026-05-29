
## require(admove); require(testthat)

make_test_dat_conf <- function(seasonal = FALSE) {
  dat <- list(
    knots_tax = matrix(1:6, nrow = 3, ncol = 2),
    knots_dif = matrix(1:4, nrow = 2, ncol = 2),
    cov = array(1:5, dim = c(5)),
    time_spline = if (seasonal) list(1:4, 1:12) else list(1)
  )

  conf <- list(
    seasonal_spline = if (seasonal) c(FALSE, TRUE) else FALSE
  )

  list(dat = dat, conf = conf)
}





test_that("default_par returns expected parameter names", {

  x <- make_test_dat_conf()
  par <- default_par(x$dat, x$conf)

  expect_true(is.list(par))
  expect_named(par, c("alpha", "beta", "gamma", "logKappa", "logSdO"))
})


test_that("default_par creates correct dimensions in non-seasonal case", {

  x <- make_test_dat_conf()
  par <- default_par(x$dat, x$conf)

  expect_equal(dim(par$alpha), c(3, 2, 1))
  expect_equal(dim(par$beta),  c(2, 2, 1))
  expect_equal(dim(par$gamma), c(2, 5, 1))
  expect_equal(dim(par$logSdO), c(2, 3))
})


test_that("default_par expands seasonal dimension correctly", {

  dat <- list(
    knots_tax = matrix(1:6, nrow = 3, ncol = 2),
    knots_dif = matrix(1:4, nrow = 2, ncol = 2),
    cov = array(1:5, dim = c(5)),
    time_spline = list(
      1:4,
      1:12
    )
  )

  conf <- list(
    seasonal_spline = c(FALSE, TRUE)
  )

  par <- default_par(dat, conf)

  expect_equal(dim(par$alpha), c(3, 2, 12))
  expect_equal(dim(par$beta),  c(2, 2, 12))
  expect_equal(dim(par$gamma), c(2, 5, 12))
})


test_that("default_par handles NULL knots and covariates", {

  dat <- list(
    knots_tax = NULL,
    knots_dif = NULL,
    cov = NULL,
    time_spline = list(1)
  )

  conf <- list(
    seasonal_spline = FALSE
  )

  par <- default_par(dat, conf)

  expect_equal(dim(par$alpha), c(1, 1, 1))
  expect_equal(dim(par$beta),  c(1, 1, 1))
  expect_equal(dim(par$gamma), c(2, 1, 1))
  expect_equal(dim(par$logSdO), c(2, 3))
})


test_that("check_par accepts valid parameter list", {

  dat <- list(
    knots_tax = matrix(1:6, nrow = 3, ncol = 2),
    knots_dif = matrix(1:4, nrow = 2, ncol = 2),
    cov = array(1:5, dim = c(5)),
    time_spline = list(1)
  )

  conf <- list(
    seasonal_spline = FALSE
  )

  par <- default_par(dat, conf)

  expect_invisible(check_par(par, dat, conf))
})


test_that("check_par errors for missing parameters", {

  x <- make_test_dat_conf()
  par <- default_par(x$dat, x$conf)
  par$beta <- NULL

  expect_error(
    check_par(par, x$dat, x$conf),
    "Missing parameter"
  )
})


test_that("check_par errors for wrong parameter dimensions", {

  x <- make_test_dat_conf()
  par <- default_par(x$dat, x$conf)
  par$alpha <- array(0, dim = c(99, 2, 1))

  expect_error(
    check_par(par, x$dat, x$conf),
    "wrong dimensions"
  )
})


test_that("check_par errors for unexpected parameters", {

  x <- make_test_dat_conf()
  par <- default_par(x$dat, x$conf)
  par$delta <- 1

  expect_error(
    check_par(par, x$dat, x$conf),
    "Unknown parameter"
  )
})
