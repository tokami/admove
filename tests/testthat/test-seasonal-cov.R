## setup_data(seasonal_cov): covariate fields that repeat every period
## (climatologies). Each layer holds its month number, so the value read at a
## tag is the slice that was used.

clim_inputs <- function(units = "month") {
  tr <- create_tref(origin = "2002-01-01", units = units)
  xc <- seq(0.05, 0.95, by = 0.1)
  arr <- array(rep(1:12, each = 100), c(10, 10, 12),
               dimnames = list(xc, xc, 0:11))
  cov <- prep_cov(arr, tref = tr)
  grid <- create_grid(cov, cellsize = 0.2)
  n <- 20
  t0 <- withr::with_seed(1, runif(n, 20, 60))
  df <- withr::with_seed(2, data.frame(
    id = seq_len(n), t0 = t0, t1 = t0 + runif(n, 1, 20),
    x0 = runif(n, 0.2, 0.8), y0 = runif(n, 0.2, 0.8),
    x1 = runif(n, 0.2, 0.8), y1 = runif(n, 0.2, 0.8)))
  tags <- suppressWarnings(prep_tags(
    df, tag_type = "c",
    names = c(id = "id", t0 = "t0", t1 = "t1", x0 = "x0", y0 = "y0",
              x1 = "x1", y1 = "y1"),
    tref = tr, sref = sref(cov)))
  list(grid = grid, cov = cov, tags = tags)
}

clim_dat <- function(...) {
  inp <- clim_inputs()
  setup_data(inp$grid, inp$cov, inp$tags, verbose = FALSE, ...)
}


test_that("seasonal covariates are read at the wrapped time", {

  dat <- clim_dat(seasonal_cov = TRUE)

  expect_equal(dat$seasonal_cov, TRUE)
  expect_equal(cov_at_tags(dat)[[1]], floor(dat$tags$t %% 12) + 1)
  ## the default knots follow the months the tags were at liberty in
  expect_gt(length(unique(dat$knots_tax[, 1])), 1L)
})


test_that("without seasonal_cov, later tags read the last slice and warn", {

  expect_warning(dat <- clim_dat(), "seasonal_cov = TRUE")
  expect_equal(dat$seasonal_cov, FALSE)
  expect_equal(cov_at_tags(dat)[[1]], rep(12, nrow(dat$tags)))
})


test_that("seasonal_cov is checked against the covariates and the period", {

  inp <- clim_inputs()
  expect_error(setup_data(inp$grid, inp$cov, inp$tags,
                          seasonal_cov = c(TRUE, FALSE), verbose = FALSE),
               "length 2")
  expect_error(setup_data(inp$grid, inp$cov, inp$tags,
                          seasonal_cov = NA, verbose = FALSE),
               "TRUE or FALSE")

  ## monthly units imply period 12; days have no default period
  np <- clim_inputs(units = "day")
  expect_error(setup_data(np$grid, np$cov, np$tags, seasonal_cov = TRUE,
                          verbose = FALSE),
               "seasonal period")
})


test_that("the configuration no longer carries seasonal_cov", {

  dat <- clim_dat(seasonal_cov = TRUE)
  conf <- default_conf(dat, verbose = FALSE)

  expect_false("seasonal_cov" %in% names(conf))
})
