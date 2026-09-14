
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


test_that("\"rtmb\" is the default smooth and matches \"natural\"", {

  skip_if_not_installed("RTMB")

  expect_equal(default_conf(skjepo$sim$dat)$smooth_method, "rtmb")

  xp <- c(0, 1, 2, 3, 4)
  yp <- c(0, 1, 0, -1, 2)
  xx <- seq(-2, 6, length.out = 50)

  for (d in c(FALSE, TRUE)) {
    a <- as.numeric(poly_fun(xp, yp, method = "rtmb", deriv = d)(xx))
    b <- as.numeric(poly_fun(xp, yp, method = "natural", deriv = d)(xx))
    expect_equal(a, b, tolerance = 1e-8)
    expect_equal(a, as.numeric(stats::splinefun(xp, yp, method = "natural")(
                                 xx, deriv = as.integer(d))),
                 tolerance = 1e-8)
  }

  ## degenerate knots are rejected the same way whatever the method
  expect_null(poly_fun(c(2, 2, 2), c(1, 2, 3), method = "rtmb"))
})


test_that("\"rtmb\" smooth composes with interpol2Dfun inside a tape", {

  ## The composition that used to abort with "'*this' is not a valid 'advector'":
  ## the spline derivative evaluated at the output of another RTMB atomic. The
  ## cached derivative tape keeps the evaluation point a bare tape variable.
  skip_if_not_installed("RTMB")

  kn <- c(21.9, 24.6, 26.5)
  set.seed(1)
  field <- matrix(runif(100, 21, 27), 10, 10)
  liv <- RTMB::interpol2Dfun(field, xlim = c(0, 1), ylim = c(0, 1))

  f <- function(p) {
    val <- poly_fun(kn, p$a, method = "rtmb")
    dval <- poly_fun(kn, p$a, method = "rtmb", deriv = TRUE)
    cov <- liv(p$x, p$y)
    val(cov)^2 + dval(cov)^2
  }

  obj <- RTMB::MakeADFun(
    f, list(a = c(0, 5, 1), x = 0.5, y = 0.5),
    map = list(a = factor(c(NA, 1, 2))), silent = TRUE
  )
  g <- obj$gr(obj$par)
  expect_true(all(is.finite(g)))
  expect_true(any(g != 0))
})


test_that("\"rtmb\" derivative is correct and its tape cache tracks length(x)", {

  skip_if_not_installed("RTMB")

  xp <- c(0, 1, 2, 3, 4)
  yv <- c(0, 1, 0, -1, 2)

  ## gradient wrt an AD evaluation point, vs finite differences
  f <- function(p) {
    dfn <- poly_fun(xp, yv, method = "rtmb", deriv = TRUE)
    ## scalar calls then a vector call: the cached tape must be rebuilt for the
    ## new length rather than reused at the wrong size
    sum(dfn(p$loc)^2) + sum(dfn(c(0.5, 1.5, 3.2) * p$loc)^2) + sum(dfn(p$loc)^2)
  }
  obj <- RTMB::MakeADFun(f, list(loc = 1.7), silent = TRUE)
  g <- as.numeric(obj$gr(obj$par))
  h <- 1e-5
  fd <- (obj$fn(1.7 + h) - obj$fn(1.7 - h)) / (2 * h)
  expect_equal(g, fd, tolerance = 1e-4)

  ## gradient wrt the knot values (the estimated parameters)
  f2 <- function(p) {
    dfn <- poly_fun(xp, p$y, method = "rtmb", deriv = TRUE)
    sum(dfn(c(0.5, 1.5, 3.2))^2)
  }
  obj2 <- RTMB::MakeADFun(f2, list(y = yv), silent = TRUE)
  g2 <- as.numeric(obj2$gr(obj2$par))
  fd2 <- sapply(seq_along(yv), function(i) {
    yu <- yl <- yv; yu[i] <- yu[i] + h; yl[i] <- yl[i] - h
    (obj2$fn(yu) - obj2$fn(yl)) / (2 * h)
  })
  expect_equal(g2, fd2, tolerance = 1e-4)
})


test_that("\"rtmb\" smooth propagates NA off-tape like \"natural\"", {

  ## Regression guard: off the tape the "rtmb" derivative falls back to
  ## stats::splinefun, whose method = "natural" branch errors on NA. The
  ## covariate interpolant returns NaN outside the field and sim_tags() relies
  ## on those gaps coming back as missing rather than as an error.
  skip_if_not_installed("RTMB")

  xp <- c(0, 1, 2, 3, 4)
  yp <- c(0, 1, 0, -1, 2)
  xx <- c(0.5, NA, 2.5, NaN)

  for (d in c(FALSE, TRUE)) {
    a <- poly_fun(xp, yp, method = "rtmb", deriv = d)(xx)
    b <- poly_fun(xp, yp, method = "natural", deriv = d)(xx)
    expect_equal(is.na(a), is.na(b))
    expect_equal(as.numeric(a)[!is.na(a)], as.numeric(b)[!is.na(b)],
                 tolerance = 1e-8)
  }

  ## all-missing input must not error either
  expect_true(all(is.na(poly_fun(xp, yp, method = "rtmb",
                                 deriv = TRUE)(c(NA_real_, NA_real_)))))
})
