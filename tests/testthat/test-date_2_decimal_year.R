

test_that("works with different dates", {

  out1 <- admove:::.date_2_decimal_year(as.Date("2026-01-01"))
  out2 <- admove:::.date_2_decimal_year(as.POSIXct("2026-01-01 00:00:00", tz = "UTC"))

  expect_equal(out1, out2)
  expect_equal(out1%%1, 0)

})


test_that("accounts for time", {

  out1 <- admove:::.date_2_decimal_year(as.POSIXct("2026-01-01 00:00:00", tz = "UTC"))
  out2 <- admove:::.date_2_decimal_year(as.POSIXct("2026-01-01 01:00:00", tz = "UTC"))

  expect_gt(out2, out1)
  expect_equal(out2%%1, 1 / (24 * 365))

})
