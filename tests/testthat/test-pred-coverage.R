## require(admove); require(testthat)


## Minimal data/conf pair carrying only what .check_pred_time_coverage() reads.
make_cov_dat <- function(pred_time, cov_time = 10:20, spline_time = 0) {
  list(
    cov = list(depth = NULL, bt = NULL),
    time_cov = list(cov_time, cov_time),
    time_spline = list(spline_time, spline_time),
    period = 12,
    pred = list(time = pred_time)
  )
}

cov_conf <- list(seasonal_cov = c(FALSE, FALSE),
                 seasonal_spline = c(FALSE, FALSE))


test_that("prediction times inside the covariate coverage pass silently", {

  expect_silent(admove:::.check_pred_time_coverage(make_cov_dat(12:18), cov_conf))
})


test_that("prediction times before the covariate coverage warn", {

  expect_warning(admove:::.check_pred_time_coverage(make_cov_dat(1:5), cov_conf),
                 "precede the covariate time coverage")

  ## the warning has to say what the consequence is, since the fields come back
  ## as a plausible zero rather than NA
  expect_warning(admove:::.check_pred_time_coverage(make_cov_dat(1:5), cov_conf),
                 "diffusion collapses to exp\\(0\\)")
})


test_that("the warning names the covariate and counts the affected times", {

  w <- tryCatch(admove:::.check_pred_time_coverage(make_cov_dat(8:14), cov_conf),
                warning = conditionMessage)

  expect_match(w, "'depth'")
  expect_match(w, "'bt'")
  ## 8 and 9 are before the first covariate time, the other five are covered
  expect_match(w, "undefined at 2 of 7 prediction time\\(s\\), from 8 to 9")
})


test_that("a time before the spline breaks is flagged too", {

  ## covariate coverage is fine, but the spline starts at 15
  dat <- make_cov_dat(12:18, cov_time = 10:20, spline_time = 15)

  expect_warning(admove:::.check_pred_time_coverage(dat, cov_conf),
                 "undefined at 3 of 7")
})


test_that("the check is a no-op without covariates or prediction times", {

  dat <- make_cov_dat(1:5)

  dat_no_cov <- dat; dat_no_cov$cov <- NULL
  expect_silent(admove:::.check_pred_time_coverage(dat_no_cov, cov_conf))

  dat_no_time <- dat; dat_no_time$pred$time <- NULL
  expect_silent(admove:::.check_pred_time_coverage(dat_no_time, cov_conf))

  dat_empty <- dat; dat_empty$pred$time <- numeric(0)
  expect_silent(admove:::.check_pred_time_coverage(dat_empty, cov_conf))
})


test_that("seasonal covariates are judged on the wrapped time", {

  ## with seasonality the lookup is on t %% period, so absolute times far beyond
  ## the covariate range still resolve
  dat <- make_cov_dat(30:35, cov_time = 0:11)
  conf_sea <- list(seasonal_cov = c(TRUE, TRUE), seasonal_spline = c(FALSE, FALSE))

  expect_silent(admove:::.check_pred_time_coverage(dat, conf_sea))
})


test_that("times after the coverage are not flagged: they clamp, not vanish", {

  ## t2index() is findInterval(), which returns the last index above the range
  ## rather than 0, so a later time reuses the last layer. That is stale but not
  ## empty, and it is the intended behaviour for a static covariate carrying a
  ## single time stamp -- warning about it would fire on every such model.
  expect_equal(t2index(30, 0:11, period = 12, seasonal = FALSE), 12L)
  expect_equal(t2index(30, 18, period = 12, seasonal = FALSE), 1L)

  expect_silent(admove:::.check_pred_time_coverage(make_cov_dat(30:35), cov_conf))
  expect_silent(
    admove:::.check_pred_time_coverage(make_cov_dat(20:25, cov_time = 18), cov_conf))
})
