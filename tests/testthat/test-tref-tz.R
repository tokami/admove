## helper: run an expression with a temporary session time zone
with_tz_session <- function(tz, expr) {
  old <- Sys.getenv("TZ", unset = NA)
  Sys.setenv(TZ = tz)
  on.exit({
    if (is.na(old)) Sys.unsetenv("TZ") else Sys.setenv(TZ = old)
  }, add = TRUE)
  force(expr)
}

test_that("create_tref reads a zone-less origin as UTC, not local time", {
  ref <- as.numeric(as.POSIXct("2003-01-01", tz = "UTC"))

  for (tz in c("UTC", "Europe/Berlin", "America/New_York")) {
    with_tz_session(tz, {
      expect_equal(as.numeric(create_tref("2003-01-01")$origin), ref)
      expect_equal(as.numeric(create_tref(as.Date("2003-01-01"))$origin), ref)
    })
  }
})

test_that("create_tref honours a time zone carried by the origin itself", {
  ref <- as.numeric(as.POSIXct("2003-01-01", tz = "UTC"))

  with_tz_session("Europe/Berlin", {
    ## trailing zone name
    expect_equal(as.numeric(create_tref("2003-01-01 UTC")$origin), ref)
    ## ISO-8601 "Z" and UTC offsets
    expect_equal(as.numeric(create_tref("2003-01-01T00:00:00Z")$origin), ref)
    expect_equal(as.numeric(create_tref("2003-01-01 01:00:00+01:00")$origin), ref)
    expect_equal(as.numeric(create_tref("2002-12-31 19:00:00-0500")$origin), ref)
    ## an explicit zone beats the tz argument
    expect_equal(as.numeric(create_tref("2003-01-01 UTC", tz = "Asia/Tokyo")$origin),
                 ref)
    ## POSIXct input keeps its own zone
    o <- as.POSIXct("2003-01-01", tz = "Asia/Tokyo")
    expect_equal(as.numeric(create_tref(o)$origin), as.numeric(o))
  })
})

test_that("tz argument sets the zone of a zone-less origin", {
  with_tz_session("UTC", {
    expect_equal(as.numeric(create_tref("2003-01-01", tz = "Asia/Tokyo")$origin),
                 as.numeric(as.POSIXct("2003-01-01", tz = "Asia/Tokyo")))
    expect_equal(as.numeric(create_tref(as.Date("2003-01-01"),
                                        tz = "Asia/Tokyo")$origin),
                 as.numeric(as.POSIXct("2003-01-01", tz = "Asia/Tokyo")))
  })
})

test_that("date_2_time accepts unit aliases", {
  d <- as.Date(c("2003-10-11", "2014-07-07"))
  tr <- create_tref("2003-01-01 UTC", units = "months")

  t_alias <- date_2_time(d, tref = tr)
  t_canon <- date_2_time(d, tref = list(origin = tr$origin, units = "month"))

  expect_equal(as.numeric(t_alias), as.numeric(t_canon))
  ## the canonical unit name is what gets stored
  expect_identical(attr(t_alias, "tref")$units, "month")

  expect_identical(
    attr(date_2_time(d, tref = list(origin = tr$origin, units = "wks")),
         "tref")$units,
    "week"
  )
})

test_that("date_2_time / time_2_date do not depend on the session time zone", {
  tr <- create_tref("2003-01-01 UTC", units = "months")
  d <- as.Date(c("2003-10-11", "2014-07-07"))

  res <- lapply(c("UTC", "Europe/Berlin", "America/New_York"), function(tz) {
    with_tz_session(tz, {
      t <- date_2_time(d, tref = tr)
      list(t = as.numeric(t),
           d = as.numeric(as.POSIXct(time_2_date(t, tref = tr))))
    })
  })

  expect_equal(res[[2]], res[[1]])
  expect_equal(res[[3]], res[[1]])

  ## and the round trip recovers the original dates
  expect_equal(as.Date(time_2_date(date_2_time(d, tref = tr), tref = tr)), d)
})

test_that("origins denoting the same instant compare equal across classes", {
  expect_true(admove:::.same_origin(as.Date("2003-01-01"),
                                    as.POSIXct("2003-01-01", tz = "UTC")))
  expect_true(admove:::.same_origin("2003-01-01 UTC",
                                    as.POSIXct("2003-01-01 09:00", tz = "Asia/Tokyo")))
  expect_false(admove:::.same_origin(as.Date("2003-01-01"),
                                     as.Date("2003-01-02")))
})
