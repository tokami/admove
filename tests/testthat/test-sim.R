
## require(admove); require(testthat)


make_sim_dat <- function(ncov = 1, nt = 2, cellsize = 0.25, units = NULL) {
  grid <- create_grid(cellsize = cellsize, verbose = FALSE)
  tref <- if (is.null(units)) NULL else list(origin = as.Date("2020-01-01"),
                                             units = units)
  cov <- lapply(seq_len(ncov), function(i) sim_cov(grid, nt = nt, tref = tref))
  names(cov) <- paste0("cov", seq_len(ncov))
  suppressWarnings(suppressMessages(
    setup_data(grid = grid, cov = cov, trange = c(0, 1), verbose = FALSE)
  ))
}


test_that("default_sim_par has the same dimensions as default_par", {

  for (ncov in 1:2) {

    dat <- make_sim_dat(ncov)
    conf <- default_conf(dat, verbose = FALSE)

    par_fit <- default_par(dat, conf, verbose = FALSE)
    par_sim <- default_sim_par(NULL, dat, conf = conf)

    for (nm in c("alpha", "beta", "gamma")) {
      expect_equal(dim(par_sim[[nm]]), dim(par_fit[[nm]]))
    }
    expect_equal(dim(par_sim$logSdO), c(2L, 3L))
  }
})


test_that("default_sim_par follows the seasonal structure in conf", {

  dat <- make_sim_dat(2, nt = 4, units = "month")
  conf <- default_conf(dat, n_seasons = 2, verbose = FALSE)

  par_fit <- default_par(dat, conf, verbose = FALSE)
  par_sim <- default_sim_par(NULL, dat, conf = conf)

  expect_equal(dim(par_sim$alpha)[3], 2L)
  expect_equal(dim(par_sim$alpha), dim(par_fit$alpha))
  expect_equal(dim(par_sim$beta), dim(par_fit$beta))
})


test_that("the simulation functions accept several covariates", {

  dat <- make_sim_dat(2)
  conf <- default_conf(dat, verbose = FALSE)

  funcs <- default_sim_funcs(dat, conf, default_sim_par(NULL, dat, conf = conf))

  expect_true(all(c("tax", "dif", "ddif") %in% names(funcs)))

  xy <- matrix(c(mean(dat$grid$xrange), mean(dat$grid$yrange)), 1, 2)
  expect_length(funcs$tax(xy, 0), 2)
  expect_true(all(is.finite(funcs$tax(xy, 0))))
})


test_that("default_sim_par overwrites the defaults and checks the dimensions", {

  dat <- make_sim_dat(1)
  conf <- default_conf(dat, verbose = FALSE)

  alpha <- array(c(0, 2, 1), dim = c(3, 1, 1))
  par <- default_sim_par(list(alpha = alpha), dat, conf = conf)

  expect_equal(par$alpha, alpha)

  expect_error(
    default_sim_par(list(alpha = array(0, dim = c(3, 2, 1))), dat, conf = conf),
    "knots x covariates x seasons"
  )
})


test_that("sim_data keeps every covariate of the data it is given", {

  dat <- make_sim_dat(2)

  set.seed(42)
  res <- sim_data(dat, simulate_cov = FALSE,
                  n_ctags = 10, n_dtags = 1, verbose = FALSE)

  expect_length(res$cov, 2)
  expect_equal(ncol(res$dat$pred$cov), 2)
  expect_equal(dim(res$par_true$alpha)[2], 2L)
  expect_equal(unclass(res$cov[[1]]), unclass(dat$cov[[1]]), ignore_attr = TRUE)
})


test_that("sim_data accepts a data or simulation object as first argument", {

  dat <- make_sim_dat(2)

  set.seed(1)
  res <- sim_data(dat, simulate_cov = FALSE,
                  n_ctags = 5, n_dtags = 1, verbose = FALSE)
  set.seed(1)
  res_named <- sim_data(dat = dat, simulate_cov = FALSE,
                        n_ctags = 5, n_dtags = 1, verbose = FALSE)

  expect_length(res$cov, 2)
  expect_equal(res$par_true, res_named$par_true)
  expect_equal(res$tags, res_named$tags)

  ## an admove_sim brings its own simulation parameters along
  res_sim <- sim_data(res, simulate_cov = FALSE,
                      n_ctags = 5, n_dtags = 1, verbose = FALSE)

  expect_equal(res_sim$par_true$alpha, res$par_true$alpha)
  expect_equal(res_sim$par_true$beta, res$par_true$beta)
})


test_that("re-simulated covariates are the ones the tags are simulated from", {

  dat <- make_sim_dat(2)

  res <- sim_data(dat, simulate_cov = TRUE,
                  n_ctags = 5, n_dtags = 1, verbose = FALSE)

  expect_length(res$cov, 2)
  expect_equal(unclass(res$dat$cov[[1]]), unclass(res$cov[[1]]),
               ignore_attr = TRUE)
  expect_false(isTRUE(all.equal(unclass(res$cov[[1]]),
                                unclass(dat$cov[[1]]),
                                check.attributes = FALSE)))
})


test_that("sim_data routes 'x' to the argument it belongs to", {

  dat <- make_sim_dat(2)

  ## a data object
  res_dat <- sim_data(dat, simulate_cov = FALSE,
                      n_ctags = 5, n_dtags = 1, verbose = FALSE)
  expect_length(res_dat$cov, 2)

  ## covariate fields: the grid follows from them and the fields are kept
  res_cov <- sim_data(dat$cov, n_ctags = 5, n_dtags = 1, verbose = FALSE)
  expect_length(res_cov$cov, 2)
  expect_equal(res_cov$grid$cellsize, dat$grid$cellsize)
  expect_equal(unclass(res_cov$cov[[1]]), unclass(dat$cov[[1]]),
               ignore_attr = TRUE)

  ## a grid
  res_grid <- sim_data(dat$grid, n_ctags = 5, n_dtags = 1, verbose = FALSE)
  expect_equal(res_grid$grid$cellsize, dat$grid$cellsize)

  ## a simulation, which brings its own parameters
  res_sim <- sim_data(res_dat, simulate_cov = FALSE,
                      n_ctags = 5, n_dtags = 1, verbose = FALSE)
  expect_equal(res_sim$par_true$alpha, res_dat$par_true$alpha)

  ## 'x' does not silently overrule an argument that was given as well
  expect_warning(
    sim_data(dat, dat = dat, simulate_cov = FALSE,
             n_ctags = 5, n_dtags = 1, verbose = FALSE),
    "was also supplied"
  )
})


test_that("arguments given explicitly win over the object in 'x'", {

  dat <- make_sim_dat(2)
  src <- sim_data(dat, simulate_cov = FALSE,
                  n_ctags = 5, n_dtags = 1, verbose = FALSE)

  ## a new grid replaces the one of the simulation, everywhere
  grid2 <- create_grid(cellsize = 0.5, verbose = FALSE)
  res <- sim_data(src, grid = grid2, simulate_cov = FALSE,
                  n_ctags = 5, n_dtags = 1, verbose = FALSE)

  expect_equal(res$grid$cellsize, grid2$cellsize)
  expect_equal(res$dat$grid$cellsize, grid2$cellsize)
  expect_length(res$cov, 2)
  ## the parameters still describe the same covariates, so they are kept
  expect_equal(res$par_true$alpha, src$par_true$alpha)

  ## a new configuration and new parameters are used as given
  conf2 <- src$conf
  conf2$use_advection <- TRUE
  expect_true(sim_data(src, conf = conf2, simulate_cov = FALSE,
                       n_ctags = 5, n_dtags = 1, verbose = FALSE)$conf$use_advection)

  par2 <- src$par_true
  par2$alpha[] <- 0
  res_par <- sim_data(src, par = par2, simulate_cov = FALSE,
                      n_ctags = 5, n_dtags = 1, verbose = FALSE)
  expect_true(all(res_par$par_true$alpha == 0))
})


test_that("parameters and configuration are dropped when the covariates change", {

  dat <- make_sim_dat(2)
  src <- sim_data(dat, simulate_cov = FALSE,
                  n_ctags = 5, n_dtags = 1, verbose = FALSE)

  expect_message(
    res <- sim_data(src, cov = dat$cov[[1]], n_ctags = 5, n_dtags = 1),
    "do not match the covariate structure"
  )

  expect_length(res$cov, 1)
  expect_equal(dim(res$par_true$alpha)[2], 1L)
  expect_equal(ncol(res$dat$knots_tax), 1L)
  expect_length(res$conf$n_seasons, 1L)
})


test_that("simulated tags keep the admove_tags class", {

  res <- sim_data(n_ctags = 5, n_dtags = 1, verbose = FALSE)

  expect_s3_class(res$tags, "admove_tags")
  expect_false(is.null(attr(res$tags, "sref")))
  expect_true(all(c("d", "c") %in% res$tags$tag_type))

  ## and are therefore accepted as the object to simulate from
  expect_s3_class(sim_data(res$tags, n_ctags = 5, n_dtags = 1, verbose = FALSE),
                  "admove_sim")
})


test_that("sim_tags resolves 'x' and 'fit' like sim_data", {

  ## release early and recapture late, so that every tag has a track to keep
  ## (sim_tags() drops mark-recapture tags that end up with a single row)
  set.seed(3)
  rel <- c(0, 0.2)
  rec <- c(0.8, 1)

  dat <- make_sim_dat(2)
  src <- sim_data(dat, simulate_cov = FALSE,
                  n_ctags = 5, n_dtags = 1, verbose = FALSE)

  ## a data object: its grid and covariates are used, not a default grid
  res_dat <- sim_tags("c", dat, n_tags = 10, trange_rel = rel, trange_rec = rec,
                      verbose = FALSE)
  expect_length(res_dat$cov, 2)
  expect_equal(res_dat$grid$cellsize, dat$grid$cellsize)

  ## a simulation brings its parameters along
  res_sim <- sim_tags("c", src, n_tags = 10, trange_rel = rel, trange_rec = rec,
                      verbose = FALSE)
  expect_equal(res_sim$par_true$alpha, src$par_true$alpha)

  ## and an argument given explicitly still wins
  grid2 <- create_grid(cellsize = 0.5, verbose = FALSE)
  res_grid <- sim_tags("c", src, grid = grid2, n_tags = 10, trange_rel = rel,
                       trange_rec = rec, verbose = FALSE)
  expect_equal(res_grid$grid$cellsize, grid2$cellsize)
  expect_equal(res_grid$dat$grid$cellsize, grid2$cellsize)
  expect_equal(res_grid$par_true$alpha, src$par_true$alpha)

  ## covariate fields as 'x': the grid follows from them
  res_cov <- sim_tags("c", dat$cov, n_tags = 10, trange_rel = rel, trange_rec = rec,
                      verbose = FALSE)
  expect_length(res_cov$cov, 2)
  expect_equal(res_cov$grid$cellsize, dat$grid$cellsize)

  ## the vignette's positional form is unchanged
  par <- list(alpha = array(c(0, 20, 10), dim = c(3, 1, 1)),
              beta = array(log(0.05), dim = c(1, 1, 1)))
  knots_tax <- matrix(round(quantile(as.numeric(dat$cov[[1]]),
                                     c(0.05, 0.5, 0.95), na.rm = TRUE), 2), 3, 1)
  res_pos <- sim_tags("c", dat$grid, dat$cov[[1]], par, n_tags = 10,
                      knots_tax = knots_tax, knots_dif = matrix(0, 1, 1),
                      trange_rel = rel, trange_rec = rec, verbose = FALSE)
  expect_equal(res_pos$par_true$alpha, par$alpha)
})


test_that("simulated observation error is estimated by the returned conf", {

  conf <- list(obs_var_type = c(0L, 0L, 0L))

  ## archival tags got noise -> estimate it
  expect_equal(admove:::.sim_obs_var_conf(conf, NULL, "d")$obs_var_type,
               c(1L, 0L, 0L))
  ## no noisy tag types -> unchanged
  expect_equal(admove:::.sim_obs_var_conf(conf, NULL, NULL)$obs_var_type,
               c(0L, 0L, 0L))
  ## an explicit user setting wins
  expect_equal(admove:::.sim_obs_var_conf(conf, conf, "d")$obs_var_type,
               c(0L, 0L, 0L))
  ## mark-recapture tags are never switched on: their kept positions carry no
  ## simulated noise, and ctags alone cannot estimate it
  expect_equal(admove:::.sim_obs_var_conf(conf, NULL, c("c"))$obs_var_type,
               c(0L, 0L, 0L))
  expect_equal(admove:::.sim_obs_var_conf(conf, NULL, c("d", "c"))$obs_var_type,
               c(1L, 0L, 0L))
  ## a stronger existing setting is kept
  conf2 <- list(obs_var_type = c(2L, 0L, 0L))
  expect_equal(admove:::.sim_obs_var_conf(conf2, NULL, "d")$obs_var_type,
               c(2L, 0L, 0L))

  sim <- small_sim()
  expect_equal(sim$conf$obs_var_type[1], 1L)
  expect_false(all(is.na(sim$map$logSdO)))
})
