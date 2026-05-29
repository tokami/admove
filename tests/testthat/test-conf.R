
## require(admove); require(testthat)



test_that("default_conf returns expected names", {

  dat <- list(
    tags = data.frame(tag_type = c("d", "s", "c")),
    cov = array(1:4, dim = c(4))
  )

  conf <- default_conf(dat, verbose = FALSE)

  expect_true(is.list(conf))
  expect_named(
    conf,
    c(
      "use_ctags",
      "use_dtags",
      "use_stags",
      "use_taxis",
      "use_advection",
      "obs_var_type",
      "do_update",
      "engine",
      "ctmc_method",
      "seasonal_cov",
      "seasonal_spline"
    )
  )
})


test_that("default_conf detects available tag types correctly", {

  dat <- list(
    tags = data.frame(tag_type = c("d", "d", "s")),
    cov = array(1:3, dim = c(3))
  )

  conf <- default_conf(dat, verbose = FALSE)

  expect_true(conf$use_dtags)
  expect_true(conf$use_stags)
  expect_false(conf$use_ctags)
})


test_that("default_conf sets sensible default variance options", {

  dat <- list(
    tags = data.frame(tag_type = c("d", "s", "c")),
    cov = array(1:2, dim = c(2))
  )

  conf <- default_conf(dat, verbose = FALSE)

  expect_equal(conf$obs_var_type, c(1, 0, 0))
  expect_equal(conf$do_update, c(TRUE,TRUE,FALSE))
  expect_equal(conf$engine, 1)
  expect_equal(conf$ctmc_method, 0)
})


test_that("default_conf creates non-seasonal defaults with correct length", {

  dat <- list(
    tags = NULL,
    cov = array(1:5, dim = c(5))
  )

  conf <- default_conf(dat, verbose = FALSE)

  expect_equal(conf$seasonal_cov, rep(FALSE, 5))
  expect_equal(conf$seasonal_spline, rep(FALSE, 5))
})


test_that("default_conf uses length 1 seasonal defaults when cov is NULL", {

  dat <- list(
    tags = NULL,
    cov = NULL
  )

  conf <- default_conf(dat, verbose = FALSE)

  expect_equal(conf$seasonal_cov, FALSE)
  expect_equal(conf$seasonal_spline, FALSE)
})


test_that("check_conf returns defaults when conf is NULL", {

  dat <- list(
    tags = data.frame(tag_type = c("d")),
    cov = array(1:3, dim = c(3))
  )

  conf1 <- check_conf(NULL, dat, verbose = FALSE)
  conf2 <- default_conf(dat, verbose = FALSE)

  expect_equal(conf1, conf2)
})


test_that("check_conf fills in missing settings", {

  dat <- list(
    tags = data.frame(tag_type = c("d")),
    cov = array(1:3, dim = c(3))
  )

  conf <- list(
    use_taxis = FALSE
  )

  conf_checked <- check_conf(conf, dat, verbose = FALSE)

  expect_false(conf_checked$use_taxis)
  expect_true("use_dtags" %in% names(conf_checked))
  expect_true("engine" %in% names(conf_checked))
  expect_true("seasonal_cov" %in% names(conf_checked))
  expect_true("seasonal_spline" %in% names(conf_checked))
})


test_that("check_conf preserves user supplied values", {

  dat <- list(
    tags = data.frame(tag_type = c("d", "c")),
    cov = array(1:2, dim = c(2))
  )

  conf <- list(
    use_taxis = FALSE,
    use_advection = TRUE,
    engine = 2
  )

  conf_checked <- check_conf(conf, dat, verbose = FALSE)

  expect_false(conf_checked$use_taxis)
  expect_true(conf_checked$use_advection)
  expect_equal(conf_checked$engine, 2)
})


test_that("check_conf keeps extra user settings unchanged", {

  dat <- list(
    tags = NULL,
    cov = NULL
  )

  conf <- list(
    my_custom_option = 123
  )

  conf_checked <- check_conf(conf, dat, verbose = FALSE)

  expect_true("my_custom_option" %in% names(conf_checked))
  expect_equal(conf_checked$my_custom_option, 123)
})


test_that("check_conf errors when conf is not a list", {

  dat <- list(
    tags = NULL,
    cov = NULL
  )

  expect_error(
    check_conf(conf = 1, dat = dat, verbose = FALSE),
    "'conf' must be a list or NULL"
  )
})
