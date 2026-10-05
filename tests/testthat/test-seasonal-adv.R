## setup_data(seasonal_adv): advection fields that repeat every period (current
## climatologies). Layer m of u holds 0.002 * m, so the slice read is visible in
## the velocity; v is zero.

sadv_inputs <- function(n_years = NULL) {
  tr <- create_tref(origin = "2002-01-01", units = "month")
  xc <- seq(0.05, 0.95, by = 0.1)
  arr <- array(rep(1:12, each = 100), c(10, 10, 12),
               dimnames = list(xc, xc, 0:11))
  cov <- prep_cov(arr, tref = tr, verbose = FALSE)
  ## the climatology, or the same twelve months repeated over n_years
  nl <- if (is.null(n_years)) 12L else 12L * n_years
  ua <- array(0.002 * rep(rep(1:12, length.out = nl), each = 100),
              c(10, 10, nl), dimnames = list(xc, xc, seq_len(nl) - 1))
  va <- ua
  va[] <- 0
  adv <- prep_adv(prep_cov(ua, tref = tr, verbose = FALSE),
                  prep_cov(va, tref = tr, verbose = FALSE), units = NULL)
  grid <- create_grid(cov, cellsize = 0.2)
  n <- 10
  t0 <- withr::with_seed(1, runif(n, 20, 40))
  df <- withr::with_seed(2, data.frame(
    id = seq_len(n), t0 = t0, t1 = t0 + runif(n, 1, 10),
    x0 = runif(n, 0.2, 0.8), y0 = runif(n, 0.2, 0.8),
    x1 = runif(n, 0.2, 0.8), y1 = runif(n, 0.2, 0.8)))
  tags <- suppressWarnings(prep_tags(
    df, tag_type = "c",
    names = c(id = "id", t0 = "t0", t1 = "t1", x0 = "x0", y0 = "y0",
              x1 = "x1", y1 = "y1"),
    tref = tr, sref = sref(cov)))
  list(grid = grid, cov = cov, tags = tags, adv = list(cur = adv))
}

sadv_dat <- function(n_years = NULL, ...) {
  inp <- sadv_inputs(n_years)
  setup_data(inp$grid, inp$cov, inp$tags, adv = inp$adv,
             seasonal_cov = TRUE, verbose = FALSE, ...)
}


test_that("seasonal advection fields are read at the wrapped time", {

  dat <- sadv_dat(seasonal_adv = TRUE)
  expect_equal(dat$seasonal_adv, c(cur = TRUE))

  gamma <- array(1, c(2, 1, 1))
  a <- admove:::.make_adv(dat$adv, dat$time_adv, gamma, NULL,
                          admove:::.get_period(dat), dat$seasonal_adv)
  xy <- cbind(0.5, 0.5)
  ## t = 25.5 is month 1.5 of the third year: slice 2
  expect_equal(unname(a$val(xy, 25.5)[1, 1]), 0.004)
  expect_equal(a$slice(25.5)[2], 2L)

  ## without the flag the last slice is used, with a warning naming the field
  expect_warning(dat0 <- sadv_dat(), "advection field 'cur'.*seasonal_adv")
  expect_equal(unname(dat0$seasonal_adv), FALSE)
  a0 <- admove:::.make_adv(dat0$adv, dat0$time_adv, gamma, NULL,
                           admove:::.get_period(dat0), dat0$seasonal_adv)
  expect_equal(unname(a0$val(xy, 25.5)[1, 1]), 0.024)
})


test_that("seasonal_adv is checked against the fields and the period", {

  inp <- sadv_inputs()
  expect_error(setup_data(inp$grid, inp$cov, inp$tags, adv = inp$adv,
                          seasonal_adv = c(TRUE, FALSE), verbose = FALSE),
               "length 2.*advection field")
  expect_error(setup_data(inp$grid, inp$cov, inp$tags,
                          seasonal_adv = TRUE, verbose = FALSE),
               "no advection fields")

  ## days have no default period
  tr <- create_tref(origin = "2002-01-01", units = "day")
  xc <- seq(0.05, 0.95, by = 0.1)
  ua <- array(0, c(10, 10, 2), dimnames = list(xc, xc, 0:1))
  fd <- prep_cov(ua, tref = tr, verbose = FALSE)
  expect_error(suppressWarnings(setup_data(
    adv = prep_adv(fd, fd, units = NULL), seasonal_adv = TRUE,
    verbose = FALSE)),
    "seasonal_adv = TRUE needs a seasonal period")
})


test_that("a seasonal field gives the likelihood of the tiled field", {

  dat_s <- sadv_dat(seasonal_adv = TRUE)
  dat_t <- sadv_dat(n_years = 8)
  expect_false(any(dat_t$seasonal_adv))

  for (engine in c("kf", "ctmc")) {
    conf <- default_conf(dat_s, verbose = FALSE)
    conf$engine <- engine
    par <- default_par(dat_s, conf, verbose = FALSE)
    par$gamma[] <- 1
    nll <- vapply(list(dat_s, dat_t), function(d) {
      f <- suppressWarnings(admove(d, conf, par, run = FALSE,
                                   verbose = FALSE))
      f$obj$fn(f$obj$par)
    }, numeric(1))
    expect_equal(nll[1], nll[2], tolerance = 1e-10, info = engine)
  }

  ## and the CTMC lattice carries the replicated slice starts
  dat_s$use_advection <- TRUE
  dat_s$period <- admove:::.get_period(dat_s)
  br <- admove:::.ctmc_breaks(dat_s)
  expect_true(all(c(24, 25, 26, 36) %in% br))
})
