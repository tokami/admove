
## require(admove); require(testthat)


make_inputs <- function(time_spline, nsea_dim, nknot = 4, ncov = 2,
                        seasonal_dif = FALSE) {

  dat <- list(time_spline = time_spline)

  conf <- list(
    use_dtags = TRUE,
    use_stags = FALSE,
    use_ctags = FALSE,
    use_advection = TRUE,
    obs_var_type = c(0L, 0L, 0L),
    seasonal_dif = seasonal_dif
  )

  par <- list(
    alpha = array(0, dim = c(nknot, ncov, nsea_dim)),
    beta = array(0, dim = c(nknot, ncov, nsea_dim)),
    gamma = array(0, dim = c(2, ncov, nsea_dim)),
    logSdO = matrix(0, 2, 3)
  )

  list(dat = dat, conf = conf, par = par)
}

as_map_array <- function(m, par_arr) {
  array(as.integer(m), dim = dim(par_arr))
}

n_estimated <- function(m) {
  v <- as.integer(m)
  length(unique(v[!is.na(v)]))
}


test_that("map entries match the parameter arrays in length", {

  x <- make_inputs(list(c(0, 3, 6, 9), 0), nsea_dim = 4)
  map <- default_map(x$dat, x$conf, x$par)

  expect_length(map$alpha, length(x$par$alpha))
  expect_length(map$beta, length(x$par$beta))
  expect_length(map$gamma, length(x$par$gamma))
})


test_that("diffusion is coupled across seasons by default", {

  x <- make_inputs(list(c(0, 3, 6, 9), 0), nsea_dim = 4)
  map <- default_map(x$dat, x$conf, x$par)

  b <- as_map_array(map$beta, x$par$beta)

  ## every non-intercept coefficient shares one id across the four seasons
  for (k in 2:4) {
    expect_equal(length(unique(b[k, 1, ])), 1L)
  }

  ## intercept stays coupled across seasons as well
  expect_equal(length(unique(b[1, 1, ])), 1L)

  ## taxis, in contrast, is estimated season by season
  a <- as_map_array(map$alpha, x$par$alpha)
  expect_equal(length(unique(a[2, 1, ])), 4L)
})


test_that("conf$seasonal_dif frees season-specific diffusion", {

  x <- make_inputs(list(c(0, 3, 6, 9), 0), nsea_dim = 4, seasonal_dif = TRUE)
  map <- default_map(x$dat, x$conf, x$par)

  b <- as_map_array(map$beta, x$par$beta)

  for (k in 2:4) {
    expect_equal(length(unique(b[k, 1, ])), 4L)
  }

  ## the intercept -- the overall diffusion level -- is freed as well
  expect_equal(length(unique(b[1, 1, ])), 4L)

  x0 <- make_inputs(list(c(0, 3, 6, 9), 0), nsea_dim = 4, seasonal_dif = FALSE)
  map0 <- default_map(x0$dat, x0$conf, x0$par)
  expect_gt(n_estimated(map$beta), n_estimated(map0$beta))
})


test_that("slices a covariate never evaluates are fixed", {

  ## covariate 1 has four spline periods, covariate 2 has one, so slices 2:4
  ## of covariate 2 never enter the likelihood
  x <- make_inputs(list(c(0, 3, 6, 9), 0), nsea_dim = 4)
  map <- default_map(x$dat, x$conf, x$par)

  a <- as_map_array(map$alpha, x$par$alpha)
  b <- as_map_array(map$beta, x$par$beta)
  g <- as_map_array(map$gamma, x$par$gamma)

  expect_true(all(is.na(a[, 2, 2:4])))
  expect_true(all(is.na(b[, 2, 2:4])))
  expect_true(all(is.na(g[, 2, 2:4])))

  ## covariate 1 keeps all of its slices
  expect_true(all(!is.na(a[2:4, 1, ])))
  expect_true(all(!is.na(g[, 1, ])))
})


test_that("non-seasonal defaults are unchanged", {

  x <- make_inputs(list(0, 0), nsea_dim = 1)
  map <- default_map(x$dat, x$conf, x$par)

  a <- as_map_array(map$alpha, x$par$alpha)
  b <- as_map_array(map$beta, x$par$beta)
  g <- as_map_array(map$gamma, x$par$gamma)

  ## first taxis coefficient fixed, the rest estimated independently
  expect_true(all(is.na(a[1, , 1])))
  expect_equal(n_estimated(map$alpha), 6L)

  ## diffusion intercept plus one coefficient per remaining knot and covariate
  expect_equal(n_estimated(map$beta), 7L)

  ## advection couples x and y within each covariate
  expect_equal(g[1, , 1], g[2, , 1])
  expect_equal(n_estimated(map$gamma), 2L)
})


test_that("seasonal_dif applies to a single-knot (constant) diffusion", {

  x <- make_inputs(list(c(0, 6), 0), nsea_dim = 2, nknot = 1)
  map0 <- default_map(x$dat, x$conf, x$par)
  expect_equal(n_estimated(map0$beta), 1L)

  x$conf$seasonal_dif <- TRUE
  map1 <- default_map(x$dat, x$conf, x$par)
  expect_equal(n_estimated(map1$beta), 2L)
})


test_that("advection is fixed when use_advection is FALSE", {

  x <- make_inputs(list(c(0, 6), 0), nsea_dim = 2)
  x$conf$use_advection <- FALSE
  map <- default_map(x$dat, x$conf, x$par)

  expect_length(map$gamma, length(x$par$gamma))
  expect_true(all(is.na(as.integer(map$gamma))))
})


test_that(".get_nsea reports one slice per covariate", {

  expect_equal(admove:::.get_nsea(list(time_spline = list(c(0, 3, 6), 0))),
               c(3L, 1L))
  expect_equal(admove:::.get_nsea(list(time_spline = NULL)), 1L)
  expect_equal(admove:::.get_nsea(list(time_spline = list())), 1L)
})
