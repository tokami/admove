## require(admove); require(testthat)


test_that(".time_labels gives dates at the resolution of the time units", {

  tr <- create_tref("2003-01-01", units = "months")
  expect_equal(admove:::.time_labels(c(48.45, 49.5, 59.45), tr),
               c("Jan 2007", "Feb 2007", "Dec 2007"))

  tr <- create_tref("2003-01-01", units = "years")
  expect_equal(admove:::.time_labels(c(0, 1, 2), tr),
               c("2003", "2004", "2005"))

  tr <- create_tref("2003-01-01", units = "days")
  expect_equal(admove:::.time_labels(c(0, 40), tr),
               c("01 Jan 2003", "10 Feb 2003"))
})


test_that(".time_labels moves to a finer format rather than repeat a label", {

  ## two times in the same year: month resolution
  tr <- create_tref("2003-01-01", units = "years")
  expect_equal(admove:::.time_labels(c(1, 1.5), tr), c("Jan 2004", "Jul 2004"))

  ## two times on the same day: time of day
  tr <- create_tref("2003-01-01", units = "days")
  expect_equal(admove:::.time_labels(c(10, 10.5), tr),
               c("11 Jan 2003 00:00", "11 Jan 2003 12:00"))

  ## the same time twice is not a clash
  tr <- create_tref("2003-01-01", units = "months")
  expect_equal(admove:::.time_labels(c(1, 1), tr), c("Feb 2003", "Feb 2003"))
})


test_that(".time_labels falls back to model time without an origin", {

  tr <- create_tref(NA, units = "months")
  expect_equal(admove:::.time_labels(c(1, 2.345), tr), c("t = 1", "t = 2.35"))
})
