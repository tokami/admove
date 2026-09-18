
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


test_that("\"rtmb\" smooth is not frozen at the taping point", {

  ## Regression guard for the cached derivative tape in .poly_fun(). D is built
  ## once per smooth, inside the active nll tape; `yp` enters it as a reference
  ## to the outer tape's nodes rather than as a numeric copy, so a new alpha must
  ## shift both the spline and its derivative when the tape is replayed. A tape
  ## that froze the knot values would keep returning the value at the taping
  ## point, which finite differences of the same tape could not detect - so the
  ## reference is the "natural" branch, which is plain R with no nested tape.
  skip_if_not_installed("RTMB")

  xp <- c(0, 1, 2, 3, 4)
  yv <- c(0, 1, 0, -1, 2)
  xx <- c(0.5, 1.5, 3.2)

  ## the taping point is a flat spline (all knot values zero, as default_par()
  ## starts alpha), so a frozen tape returns exactly zero for every yv below
  f <- function(p) sum(poly_fun(xp, p$y, method = "rtmb", deriv = TRUE)(xx))
  obj <- RTMB::MakeADFun(f, list(y = rep(0, 5)), silent = TRUE)

  ref <- sum(stats::splinefun(xp, yv, method = "natural")(xx, deriv = 1))
  expect_equal(obj$fn(yv), ref, tolerance = 1e-8)
  expect_false(isTRUE(all.equal(ref, 0)))

  ## and again far from both the taping point and the first replay
  yv2 <- c(-1, 0.5, 2, 0, 1)
  expect_equal(obj$fn(yv2),
               sum(stats::splinefun(xp, yv2, method = "natural")(xx, deriv = 1)),
               tolerance = 1e-8)
})


test_that("perturbed alpha gives the same nll and gradient under both smooths", {

  ## Same guard at model level: the taxis drift is kappa * d/dx of the taxis
  ## smooth, so a derivative tape stuck at the taping point (alpha = 0, i.e. no
  ## drift at all) would change the nll. "natural" builds the same spline without
  ## a nested tape, so it is an independent reference for both value and gradient.
  skip_if_not_installed("RTMB")

  obj_rtmb <- tiny_obj("rtmb")
  obj_nat <- tiny_obj("natural")

  p0 <- obj_rtmb$par
  expect_equal(names(p0), names(obj_nat$par))
  ia <- which(names(p0) == "alpha")
  expect_gt(length(ia), 0)

  ## far enough from alpha = 0 to make the taxis drift bite, but still inside the
  ## region where the KF on a single archival tag returns a finite nll
  p1 <- p0
  p1[ia] <- c(0.2, -0.3)[seq_along(ia)]
  p1[names(p1) == "beta"] <- 0.4
  p1[names(p1) == "logSdO"] <- p0[names(p0) == "logSdO"] + 0.2

  ## the perturbation must actually move the smooth, or the test is vacuous
  expect_false(isTRUE(all.equal(obj_rtmb$fn(p0), obj_rtmb$fn(p1))))

  expect_equal(obj_rtmb$fn(p1), obj_nat$fn(p1), tolerance = 1e-8)

  g_rtmb <- as.numeric(obj_rtmb$gr(p1))
  g_nat <- as.numeric(obj_nat$gr(p1))
  expect_true(all(is.finite(g_rtmb)))
  expect_false(isTRUE(all.equal(g_rtmb[ia], rep(0, length(ia)))))
  expect_equal(g_rtmb, g_nat, tolerance = 1e-6)

  ## "natural" has no nested tape, so its own gradient can be checked against
  ## finite differences of its value - which pins the "rtmb" gradient too
  h <- 1e-5
  fd <- sapply(seq_along(p1), function(i) {
    pu <- pl <- p1
    pu[i] <- pu[i] + h
    pl[i] <- pl[i] - h
    (obj_nat$fn(pu) - obj_nat$fn(pl)) / (2 * h)
  })
  expect_equal(g_nat, fd, tolerance = 1e-4)
})
