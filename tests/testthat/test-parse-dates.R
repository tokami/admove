## require(admove); require(testthat)

## prep_tags()/prep_cov() date parsing keeps the time of day and reports
## unparseable dates


test_that(".parse_dates keeps the time of day of character date-times", {

  out <- admove:::.parse_dates(c("2/13/2025 9:00", "2/13/2025 12:30"),
                               date_format = "%m/%d/%Y %H:%M")
  expect_s3_class(out, "POSIXct")
  expect_equal(format(out, "%Y-%m-%d %H:%M", tz = "UTC"),
               c("2025-02-13 09:00", "2025-02-13 12:30"))
})


test_that(".parse_dates keeps fractional days of numeric dates", {

  out <- admove:::.parse_dates(c(45000, 45000.5), date_origin = "1899-12-30")
  expect_s3_class(out, "POSIXct")
  expect_equal(as.numeric(diff(out), units = "hours"), 12)
})


test_that(".parse_dates returns Dates when no time of day is present", {

  expect_s3_class(admove:::.parse_dates(c("2/13/2025", "2/14/2025"),
                                        date_format = "%m/%d/%Y"), "Date")
  expect_s3_class(admove:::.parse_dates(c(45000, 45001),
                                        date_origin = as.Date("1899-12-30")), "Date")
  expect_s3_class(admove:::.parse_dates(c("2025-02-13 00:00", "2025-02-14 00:00"),
                                        date_format = "%Y-%m-%d %H:%M"), "Date")
})


test_that(".parse_dates warns about dates that do not match the format", {

  expect_warning(
    admove:::.parse_dates(c("2/13/2025 9:00", "12/2/2025 9:00"),
                          date_format = "%d/%m/%Y %H:%M"),
    "1 of 2 dates could not be parsed"
  )
})


test_that("prep_tags keeps sub-daily positions of archival tags", {

  d <- data.frame(time = c("2/5/2025 12:00", "2/6/2025 0:00",
                           "2/6/2025 10:25", "2/6/2025 11:00"),
                  lon = c(-66.975, -66.9, -66.825, -66.825),
                  lat = c(17.65, 17.6, 17.575, 17.575),
                  id = 1)
  tags <- prep_tags(d, tag_type = "d",
                    names = c(t = "time", x = "lon", y = "lat", id = "id"),
                    date_format = "%m/%d/%Y %H:%M",
                    tref = create_tref("2025-01-01 UTC", units = "day"),
                    verbose = FALSE)

  expect_false(anyNA(tags$t))
  expect_false(any(duplicated(tags$t)))
  expect_equal(time_2_date(tags$t, tref = tref(tags)),
               as.POSIXct(c("2025-02-05 12:00", "2025-02-06 00:00",
                            "2025-02-06 10:25", "2025-02-06 11:00"), tz = "UTC"))
})
