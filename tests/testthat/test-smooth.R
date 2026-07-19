
## require(admove); require(testthat)

## Internal helper under test
poly_fun <- admove:::.poly_fun


test_that("natural spline interpolates the knot values exactly", {

  xp <- c(0, 1, 2, 3, 4)
  yp <- c(0, 1, 0, -1, 2)

  f <- poly_fun(xp, yp, method = "natural")

  expect_equal(as.numeric(f(xp)), yp, tolerance = 1e-8)
})


test_that("natural spline analytic derivative matches finite differences", {

  xp <- c(0, 2, 4, 6)
  yp <- c(1, -0.5, 0.5, 2)

  f <- poly_fun(xp, yp, method = "natural")
  df <- poly_fun(xp, yp, method = "natural", deriv = TRUE)

  xx <- c(0.7, 1.9, 3.3, 5.1)
  h <- 1e-5
  fd <- (f(xx + h) - f(xx - h)) / (2 * h)

  expect_equal(as.numeric(df(xx)), as.numeric(fd), tolerance = 1e-5)
})


test_that("natural spline extrapolates linearly (constant tail slope)", {

  xp <- c(0, 2, 4)
  yp <- c(0, 1, -0.5)

  df <- poly_fun(xp, yp, method = "natural", deriv = TRUE)

  ## below the first knot
  left <- as.numeric(df(c(-3, -2, -1)))
  ## above the last knot
  right <- as.numeric(df(c(5, 6, 7)))

  expect_equal(left, rep(left[1], 3), tolerance = 1e-8)
  expect_equal(right, rep(right[1], 3), tolerance = 1e-8)
})


test_that("legacy polynomial does not extrapolate linearly", {

  xp <- c(0, 2, 4)
  yp <- c(0, 1, -0.5)

  df <- poly_fun(xp, yp, method = "poly", deriv = TRUE)

  right <- as.numeric(df(c(5, 6, 7)))

  ## quadratic through 3 knots -> derivative is not constant in the tail
  expect_false(isTRUE(all.equal(right, rep(right[1], 3))))
})


test_that("special cases are method independent", {

  ## single knot -> constant (returns the knot value; scalar recycles when added
  ## to a location vector inside .make_habi)
  fc <- poly_fun(c(1.5), c(2.3), method = "natural")
  expect_equal(as.numeric(fc(-1)), 2.3)
  expect_equal(as.numeric(fc(5)), 2.3)

  ## advection -> linear through the origin
  fa <- poly_fun(NULL, 0.7, adv = TRUE, method = "natural")
  expect_equal(as.numeric(fa(c(0, 1, 2))), c(0, 0.7, 1.4))
})


test_that("natural spline is AD-differentiable through knot values and eval point", {

  skip_if_not_installed("RTMB")

  xp <- c(0, 1, 2, 3, 4)

  ## gradient wrt knot values (the estimated parameters)
  f1 <- function(p) {
    fn <- poly_fun(xp, p$y, method = "natural")
    sum(fn(c(0.5, 1.5, 3.2)))
  }
  obj1 <- RTMB::MakeADFun(f1, list(y = c(0, 1, 0, -1, 2)), silent = TRUE)
  g1 <- obj1$gr(obj1$par)
  expect_true(all(is.finite(g1)))
  expect_true(any(g1 != 0))

  ## gradient wrt an AD-typed evaluation point (the KF case), checked vs FD
  f2 <- function(p) {
    fn <- poly_fun(xp, c(0, 1, 0, -1, 2), method = "natural")
    dfn <- poly_fun(xp, c(0, 1, 0, -1, 2), method = "natural", deriv = TRUE)
    fn(p$loc) + dfn(p$loc)
  }
  obj2 <- RTMB::MakeADFun(f2, list(loc = 1.7), silent = TRUE)
  g2 <- as.numeric(obj2$gr(obj2$par))
  h <- 1e-5
  fd <- (obj2$fn(1.7 + h) - obj2$fn(1.7 - h)) / (2 * h)
  expect_equal(g2, fd, tolerance = 1e-4)
})


test_that("natural spline composes with interpol2Dfun inside a tape", {

  ## Regression guard: inside nll() the spline is evaluated at the output of the
  ## interpol2Dfun habitat interpolator (itself an RTMB atomic). The hand-rolled
  ## truncated-power spline must tape cleanly in that composition.
  skip_if_not_installed("RTMB")

  kn <- c(21.9, 24.6, 26.5)
  set.seed(1)
  field <- matrix(runif(100, 21, 27), 10, 10)
  liv <- RTMB::interpol2Dfun(field, xlim = c(0, 1), ylim = c(0, 1))

  f <- function(p) {
    val <- poly_fun(kn, p$a, method = "natural")
    dval <- poly_fun(kn, p$a, method = "natural", deriv = TRUE)
    cov <- liv(p$x, p$y)
    val(cov)^2 + dval(cov)^2
  }

  obj <- RTMB::MakeADFun(
    f, list(a = c(0, 5, 1), x = 0.5, y = 0.5),
    map = list(a = factor(c(NA, 1, 2))), silent = TRUE
  )
  g <- obj$gr(obj$par)
  expect_true(all(is.finite(g)))
})
