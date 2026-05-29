
## require(admove); require(testthat)



test_that("default_map returns expected named factor list", {

  par <- list(
    alpha = array(0, dim = c(3, 1, 2)),
    beta = array(0, dim = c(3, 1, 2)),
    gamma = matrix(0, nrow = 2, ncol = 2)
  )

  conf <- list(
    use_advection = TRUE,
    use_dtags = TRUE,
    use_stags = TRUE,
    obs_var_type = c(TRUE, TRUE)
  )

  dat <- list()

  map <- default_map(dat, conf, par)

  expect_type(map, "list")
  expect_named(map, c("alpha", "logKappa", "beta", "gamma", "logSdO"))

  expect_s3_class(map$alpha, "factor")
  expect_s3_class(map$logKappa, "factor")
  expect_s3_class(map$beta, "factor")
  expect_s3_class(map$gamma, "factor")
  expect_s3_class(map$logSdO, "factor")
})


test_that("default_map couples gamma by covariate when advection is enabled", {

  par <- list(
    alpha = array(0, dim = c(3, 1, 2)),
    beta = array(0, dim = c(3, 1, 2)),
    gamma = matrix(0, nrow = 2, ncol = 3)
  )

  conf <- list(
    use_advection = TRUE,
    use_dtags = FALSE,
    use_stags = FALSE,
    obs_obs_type = c(FALSE, FALSE)
  )

  dat <- list()

  map <- default_map(dat, conf, par)

  ## Expected layout from:
  ## factor(sapply(1:ncol(par$gamma), function(x) rep(x, nrow(par$gamma))))
  expect_equal(
    as.integer(map$gamma),
    c(1, 1, 2, 2, 3, 3)
  )
})


test_that("default_map fixes all gamma parameters when advection is disabled", {

  par <- list(
    alpha = array(0, dim = c(3, 1, 2)),
    beta = array(0, dim = c(3, 1, 2)),
    gamma = matrix(0, nrow = 2, ncol = 3)
  )

  conf <- list(
    use_advection = FALSE,
    use_dtags = FALSE,
    use_stags = FALSE,
    obs_var_type = c(FALSE, FALSE)
  )

  dat <- list()

  map <- default_map(dat, conf, par)

  expect_true(all(is.na(map$gamma)))
  expect_length(map$gamma, length(par$gamma))
})


test_that("default_map sets default logSdO mapping for d-tags and s-tags", {

  par <- list(
    alpha = array(0, dim = c(3, 1, 2)),
    beta = array(0, dim = c(3, 1, 2)),
    gamma = matrix(0, nrow = 2, ncol = 1)
  )

  conf <- list(
    use_advection = FALSE,
    use_dtags = TRUE,
    use_stags = TRUE,
    obs_var_type = c(TRUE, TRUE)
  )

  dat <- list()

  map <- default_map(dat, conf, par)

  ## d-tags: positions 1:2 coupled as 1
  ## s-tags: positions 3:4 coupled as 2
  ## remaining positions fixed
  expect_equal(
    as.integer(map$logSdO),
    c(1, 1, 2, 2, NA, NA)
  )
})


test_that("default_map sets default logSdO mapping for d-tags and s-tags", {

  par <- list(
    alpha = array(0, dim = c(3, 1, 2)),
    beta = array(0, dim = c(3, 1, 2)),
    gamma = matrix(0, nrow = 2, ncol = 1)
  )

  conf <- list(
    use_advection = FALSE,
    use_dtags = TRUE,
    use_stags = TRUE,
    obs_var_type = c(TRUE, TRUE)
  )

  dat <- list()

  map <- default_map(dat, conf, par)

  ## d-tags: positions 1:2 coupled as 1
  ## s-tags: positions 3:4 coupled as 2
  ## remaining positions fixed
  expect_equal(
    as.integer(map$logSdO),
    c(1, 1, 2, 2, NA, NA)
  )
})
