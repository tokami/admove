## require(admove); require(testthat)

## Covariate time slices are labelled by their START time, so slice j covers
## [t_j, t_{j+1}). A time exactly on a start must use the slice starting there:
## dates convert to exact slice starts (first of the month in monthly units).


test_that("t2index assigns a time on a slice start to that slice", {

  tv <- c(0, 1, 2, 3)

  expect_equal(admove:::t2index(c(0, 0.5, 1, 1.5, 2, 2.999, 3, 3.5), tv),
               c(1, 1, 2, 2, 3, 3, 4, 4))
  ## before the first slice: undefined; after the last: clamped
  expect_equal(admove:::t2index(c(-0.1, 10), tv), c(0, 4))
})


test_that("seasonal indices start a new season exactly at its breakpoint", {

  ts <- c(0, 3, 6, 9)   ## four seasons over a period of 12

  expect_equal(admove:::t2index(c(0, 2.9, 3, 6, 9, 11.9, 12, 15), ts,
                                period = 12, seasonal = TRUE),
               c(1, 1, 2, 3, 4, 4, 1, 2))
})


test_that("the first of a month uses the covariate layer of that month", {

  dates <- as.Date(c("2020-01-15", "2020-03-01", "2020-03-31"))
  t <- date_2_time(dates, tref = list(origin = as.Date("2020-01-01"),
                                      units = "month"))
  layers <- 0:11   ## monthly layers labelled by month start

  expect_equal(admove:::t2index(as.numeric(t), layers), c(1, 3, 3))
})
