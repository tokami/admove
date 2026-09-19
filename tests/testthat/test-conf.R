
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
      "drift_scheme",
      "n_seasons",
      "seasonal_cov",
      "seasonal_spline",
      "seasonal_dif",
      "smooth_method",
      "kf_boundary"
    )
  )
})


test_that("default_conf switches seasonal diffusion off by default", {

  dat <- list(
    tags = NULL,
    cov = array(1:3, dim = c(3))
  )

  conf <- default_conf(dat, verbose = FALSE)

  expect_false(conf$seasonal_dif)
})


test_that("check_conf rejects an invalid seasonal_dif", {

  dat <- list(
    tags = NULL,
    cov = array(1:3, dim = c(3))
  )

  expect_error(
    check_conf(list(seasonal_dif = c(TRUE, FALSE)), dat, verbose = FALSE),
    "single TRUE or FALSE"
  )
  expect_error(
    check_conf(list(seasonal_dif = NA), dat, verbose = FALSE),
    "single TRUE or FALSE"
  )
})


test_that("default_conf sets natural spline as default smooth", {

  dat <- list(
    tags = NULL,
    cov = array(1:5, dim = c(5))
  )

  conf <- default_conf(dat, verbose = FALSE)

  expect_equal(conf$smooth_method, "rtmb")
})


test_that("check_conf accepts valid smooth_method and rejects invalid", {

  dat <- list(
    tags = data.frame(tag_type = c("d")),
    cov = array(1:3, dim = c(3))
  )

  for (m in c("rtmb", "natural", "poly")) {
    expect_equal(
      check_conf(list(smooth_method = m), dat, verbose = FALSE)$smooth_method,
      m
    )
  }

  expect_error(
    check_conf(list(smooth_method = "bspline"), dat, verbose = FALSE),
    "smooth_method"
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

  expect_equal(conf$obs_var_type, rep("none", 3))
  expect_equal(conf$do_update, c(TRUE,TRUE,FALSE))
  expect_equal(conf$engine, "kf")
  expect_equal(conf$ctmc_method, "expav")
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
  ## numbers are accepted and normalised to the names
  expect_equal(conf_checked$engine, "ctmc")
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


test_that("check_conf warns when observation error is estimated from ctags alone", {

  dat <- skjepo$sim$dat
  conf <- default_conf(dat, verbose = FALSE)
  conf$obs_var_type[3] <- 1L

  ## archival and mark-recapture tags: no warning
  expect_no_warning(check_conf(conf, dat, verbose = FALSE))

  ## mark-recapture tags only
  dat_c <- dat
  dat_c$tags <- dat$tags[dat$tags$tag_type == "c", ]
  expect_warning(check_conf(conf, dat_c, verbose = FALSE),
                 "only mark-recapture tags")

  conf$obs_var_type[3] <- 0L
  expect_no_warning(check_conf(conf, dat_c, verbose = FALSE))
})


test_that("ctmc_method takes the method names and rejects the old numbers", {

  expect_silent(admove:::.check_ctmc_method("expm"))
  expect_silent(admove:::.check_ctmc_method("expav"))
  expect_error(admove:::.check_ctmc_method("uniformization"), "must be either")

  ## the old numbering is not remapped: a 1 written for it would otherwise
  ## silently select the other method
  expect_error(admove:::.check_ctmc_method(0), 'replace 0 with "expm"')
  expect_error(admove:::.check_ctmc_method(1), 'replace 1 with "expav"')
  expect_error(admove:::.check_ctmc_method(2), "now given as a string")

  conf <- default_conf(skjepo$sim$dat, verbose = FALSE)
  conf$ctmc_method <- 1
  expect_error(check_conf(conf, skjepo$sim$dat, verbose = FALSE),
               'replace 1 with "expav"')
})


test_that("both CTMC matrix-exponential methods give the same likelihood", {

  grid <- create_grid(xrange = c(0, 1), yrange = c(0, 1), cellsize = 0.25,
                      verbose = FALSE)
  sim <- withr::with_seed(1, suppressMessages(suppressWarnings(
    sim_data(grid = grid, n_dtags = 2, n_ctags = 5, trange = c(0, 1),
             verbose = FALSE)
  )))

  nll <- function(method) {
    conf <- sim$conf
    conf$engine <- 2
    conf$ctmc_method <- method
    obj <- suppressMessages(suppressWarnings(
      admove(sim$dat, conf, sim$par, sim$map, run = FALSE, verbose = FALSE)
    ))$obj
    obj$fn(obj$par)
  }

  expect_equal(nll("expav"), nll("expm"), tolerance = 1e-6)
})


test_that("engine and obs_var_type are configured by name, numbers still work", {

  dat <- skjepo$sim$dat

  conf <- default_conf(dat, verbose = FALSE)
  expect_equal(conf$engine, "kf")
  expect_equal(conf$obs_var_type, rep("none", 3))

  ## names, in any case, and the numbers they replaced
  expect_equal(admove:::.get_engine_integer("ctmc"), 2L)
  expect_equal(admove:::.get_engine_integer("CTMC"), 2L)
  expect_equal(admove:::.get_engine_name(2), "ctmc")
  expect_equal(admove:::.get_engine_name("kf"), "kf")

  expect_equal(admove:::.get_obs_var_type_integer(rep("none", 3)), rep(0L, 3))
  expect_equal(admove:::.get_obs_var_type_name(c(1, 2, 3)),
               c("all_but_last", "all", "data"))

  ## a number assigned into a character setting arrives as "1"
  conf$obs_var_type[1] <- 1L
  expect_equal(admove:::.get_obs_var_type_integer(conf$obs_var_type),
               c(1L, 0L, 0L))
  expect_equal(check_conf(conf, dat, verbose = FALSE)$obs_var_type,
               c("all_but_last", "none", "none"))

  ## check_conf normalises whatever spelling was used
  conf2 <- check_conf(list(engine = 2, obs_var_type = c(0, 0, 2)), dat,
                      verbose = FALSE)
  expect_equal(conf2$engine, "ctmc")
  expect_equal(conf2$obs_var_type, c("none", "none", "all"))

  expect_error(check_conf(list(engine = "kalman"), dat, verbose = FALSE),
               "must be \"kf\"")
  expect_error(check_conf(list(obs_var_type = c("none", "none", "sometimes")),
                          dat, verbose = FALSE),
               "must be one of")
})


test_that("the map follows obs_var_type given by name", {

  dat <- skjepo$sim$dat
  conf <- default_conf(dat, verbose = FALSE)
  par <- suppressMessages(default_par(dat, conf, verbose = FALSE))

  conf_name <- conf
  conf_name$obs_var_type <- c("all", "none", "none")
  conf_num <- conf
  conf_num$obs_var_type <- c(2L, 0L, 0L)

  expect_equal(default_map(dat, conf_name, par)$logSdO,
               default_map(dat, conf_num, par)$logSdO)
})
