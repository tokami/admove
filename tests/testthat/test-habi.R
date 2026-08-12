## require(admove); require(testthat)


## Minimal stand-in for a fitted object carrying everything .get_habi() needs
## to rebuild the habitat objects from scratch.
make_habi_fit <- function() {

  nx <- 5
  ny <- 4
  xgr <- seq(0, 1, length.out = nx)
  ygr <- seq(0, 1, length.out = ny)
  cov1 <- array(outer(xgr, ygr, "+"), dim = c(nx, ny, 1),
                dimnames = list(xgr, ygr, 0))

  dat <- list(
    cov = list(cov1 = cov1),
    xrange_cov = matrix(range(xgr), 1, 2),
    yrange_cov = matrix(range(ygr), 1, 2),
    time_cov = list(0),
    time_spline = list(0),
    knots_tax = matrix(c(0, 1, 2), 3, 1),
    knots_dif = matrix(c(0, 1, 2), 3, 1),
    period = 1,
    pred = list(grid = list(xygrid = cbind(x = c(0.25, 0.5),
                                           y = c(0.25, 0.5))),
                time = 0)
  )

  conf <- list(seasonal_cov = FALSE, seasonal_spline = FALSE)

  par <- list(alpha = array(c(0, 0.5, 1), c(3, 1, 1)),
              beta = array(c(0, 0.2, 0.4), c(3, 1, 1)),
              gamma = array(0, c(2, 1, 1)),
              logKappa = 0)

  x <- list(dat = dat,
            conf = conf,
            par = par,
            map = list(),
            opt = list(par = numeric(0)))

  attr(x, "tref") <- structure(list(origin = 0, units = "day", period = 1),
                               class = "admove_tref")

  x$pred <- list(habi = admove:::.build_habi(dat, conf, par, 1)$habi)

  x
}


test_that(".get_habi returns the stored objects when they are alive", {

  x <- make_habi_fit()

  expect_true(admove:::.habi_usable(x$pred$habi, x))
  expect_identical(admove:::.get_habi(x), x$pred$habi)
})


test_that("habi objects are detected as stale after a save/load round trip", {

  x <- make_habi_fit()

  ## RTMB's interpolants hold an external pointer, which serialisation drops
  x2 <- unserialize(serialize(x, NULL))

  expect_false(admove:::.habi_usable(x2$pred$habi, x2))
  expect_error(x2$pred$habi$tax$grad(x2$dat$pred$grid$xygrid, 0),
               "external pointer")
})


test_that(".get_habi rebuilds stale habi objects and reproduces the values", {

  x <- make_habi_fit()
  xy <- x$dat$pred$grid$xygrid

  val0 <- x$pred$habi$tax$val(xy, 0)
  grad0 <- x$pred$habi$tax$grad(xy, 0)

  x2 <- unserialize(serialize(x, NULL))
  habi <- admove:::.get_habi(x2)

  expect_equal(habi$tax$val(xy, 0), val0)
  expect_equal(habi$tax$grad(xy, 0), grad0)
})


test_that(".get_habi errors informatively without covariates", {

  x <- make_habi_fit()
  x$dat$cov <- NULL
  x$pred$habi <- NULL

  expect_error(admove:::.get_habi(x), "no covariates")
})
