test_that("time_2_date inverts date_2_time for fixed-length units", {
  origin <- as.POSIXct("2010-03-01 00:00:00", tz = "UTC")
  set.seed(1)
  d <- origin + runif(500, 0, 400 * 86400)

  for (u in c("second", "minute", "hour", "day", "week")) {
    tr <- list(origin = origin, units = u)
    t <- date_2_time(d, tref = tr)
    drec <- time_2_date(as.numeric(t), tref = tr)
    expect_lt(max(abs(as.numeric(difftime(drec, d, units = "secs")))), 1e-3)
  }
})

test_that("time_2_date round-trips date_2_time for month/year (calendar-aware)", {
  tr <- create_tref("2003-01-01 UTC", units = "months")
  set.seed(2)
  d <- tr$origin + runif(1000, 0, 15 * 365.25) * 86400

  for (u in c("month", "year")) {
    ref <- list(origin = tr$origin, units = u)
    t <- date_2_time(d, tref = ref)
    drec <- time_2_date(as.numeric(t), tref = ref)
    ## forward-consistent to machine precision (true inverse of date_2_time)
    expect_lt(max(abs(as.numeric(date_2_time(drec, tref = ref)) - as.numeric(t))), 1e-6)
    ## day-resolution dates are recovered essentially exactly
    d_day <- as.Date(d)
    t_day <- date_2_time(d_day, tref = ref)
    drec_day <- as.Date(time_2_date(as.numeric(t_day), tref = ref))
    expect_true(all(abs(as.numeric(drec_day - d_day)) <= 1))
  }
})

test_that("time_2_date accepts an admove_tref and a bare origin", {
  tr <- create_tref("2003-01-01 UTC", units = "months")
  out <- time_2_date(c(0, 9, 12), tref = tr)
  expect_s3_class(out, "POSIXct")
  expect_length(out, 3)
  expect_equal(as.numeric(out[1]), as.numeric(tr$origin))

  ## bare POSIXt origin -> units default to "day"
  o <- as.POSIXct("2020-01-01", tz = "UTC")
  expect_equal(time_2_date(2, tref = o), o + as.difftime(2, units = "days"))
})

test_that("time_2_date preserves NA and errors without a reference", {
  tr <- create_tref("2003-01-01 UTC", units = "months")
  out <- time_2_date(c(1, NA, 3), tref = tr)
  expect_true(is.na(out[2]))
  expect_false(anyNA(out[c(1, 3)]))
  expect_error(time_2_date(1), "tref")
})

test_that("time_2_date dispatches on an admove_tags object", {
  ## minimal tags built through the package pipeline
  set.seed(3)
  n <- 20
  df <- data.frame(
    id = seq_len(n),
    t0 = as.Date("2005-01-01") + sample(0:200, n),
    x0 = runif(n, 40, 60), y0 = runif(n, -10, 5)
  )
  df$t1 <- df$t0 + sample(10:200, n, replace = TRUE)
  df$x1 <- df$x0 + rnorm(n)
  df$y1 <- df$y0 + rnorm(n)

  tags <- prep_tags(df, "c",
                    names = c(t0 = "t0", t1 = "t1", x0 = "x0", x1 = "x1",
                              y0 = "y0", y1 = "y1", id = "id"),
                    date_format = "%Y-%m-%d",
                    sref = create_sref(crs = 4326))
  tags <- add_tref(tags, create_tref("2005-01-01 UTC", units = "months"),
                   shift_origin = TRUE)

  dts <- time_2_date(tags)
  expect_s3_class(dts, "POSIXct")
  expect_length(dts, length(tags$t))
  ## converting the stored numeric time via the object == via t + tref
  expect_equal(dts, time_2_date(as.numeric(tags$t), tref = tref(tags)))
})
