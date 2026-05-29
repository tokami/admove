
## require(admove); require(testthat)

cov1 <- sim_cov()
cov2 <- sim_cov()


test_that("default tref", {

  tr <- create_tref()

  expect_equal(tr$origin, as.POSIXct(NA))
  expect_equal(tr$units, NA_character_)
  expect_equal(tr$period, NA_integer_)

})


test_that("default periods", {

  tr <- create_tref(units = "year")
  expect_equal(tr$period, 1)

  tr <- create_tref(units = "semester")
  expect_equal(tr$period, 2)

  tr <- create_tref(units = "quarter")
  expect_equal(tr$period, 4)

  tr <- create_tref(units = "month")
  expect_equal(tr$period, 12)

  tr <- create_tref(units = "week")
  expect_equal(tr$period, 52)

  tr <- create_tref(units = "day")
  expect_equal(tr$period, NA_real_) ## could be 365 and 365.25 or something else

})


test_that("custom periods", {

  tr <- create_tref(period = 10)
  expect_equal(tr$period, 10)

})



test_that("auto recognize units via period", {

  tr <- create_tref(period = 12)
  expect_equal(tr$units, "month")

  tr <- create_tref(period = 10)
  expect_equal(tr$units, "custom")

})




tr <- create_tref(origin = as.Date("2025-01-01"), units = "month")
tref(cov1) <- tr

test_that("extract and set trefs", {

  expect_error(tref(tr))
  expect_error(tref(NA))

  expect_equal(tr, tref(cov1))

  tr <- tref(cov1)

  expect_equal(tr, tref(cov1))

})



cov3 <- add_tref(cov2, cov1)

test_that("copy trefs", {

  expect_equal(tref(cov3), tref(cov1))

})



test_that("origins, units, periods", {

  expect_equal(origin(tr), as.POSIXct("2025-01-01", tz = "UTC"))
  expect_equal(units_time(tr), "month")
  expect_equal(period(tr), 12)

  expect_equal(origin(cov3), as.POSIXct("2025-01-01", tz = "UTC"))
  expect_equal(units_time(cov3), "month")
  expect_equal(period(cov3), 12)

})


test_that("add trefs", {

  out <- add_tref(cov2, list(origin = as.Date("2025-01-01")))
  expect_equal(origin(out), as.POSIXct("2025-01-01", tz = "UTC"))

  out <- add_tref(cov2, list(units = "month"))
  expect_equal(units_time(out), "month")
  expect_equal(period(out), 12)

  out <- add_tref(cov2, tr)
  expect_equal(tref(out), tr)

})



test_that("compare trefs", {

  out <- tref_equal(tref(cov3), tref(cov1))
  expect_true(out)

  out <- tref_equal(tref(cov2), tref(cov1))
  expect_false(out)

})



test_that("scale trefs", {

  out <- scale_tref(cov3, 1)
  expect_equal(tref(out), tref(cov3))

  out <- scale_tref(cov3, 0.1)
  expect_equal(period(out), 1.2)

  ## TODO: test actual time vectors (also for data, tags, etc.)

})
