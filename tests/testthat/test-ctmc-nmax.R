## conf$ctmc_nmax: cap on the expAv uniformization terms, and the post-fit
## check that the cap does not bind at the estimates.


## Same shape as ctmc_pin_sim_obj() in test-ctmc-events.R: 2 archival and 5
## conventional tags on a 4 x 4 grid.
nmax_sim <- function() {
  grid <- create_grid(xrange = c(0, 1), yrange = c(0, 1), cellsize = 0.25,
                      verbose = FALSE)
  sim <- withr::with_seed(1, suppressMessages(suppressWarnings(
    sim_data(grid = grid, n_dtags = 2, n_ctags = 5, trange = c(0, 1),
             verbose = FALSE)
  )))
  sim$conf$engine <- "ctmc"
  sim
}

nmax_obj <- function(nmax) {
  sim <- nmax_sim()
  sim$conf["ctmc_nmax"] <- list(nmax)
  suppressMessages(suppressWarnings(
    admove(sim$dat, sim$conf, sim$par, sim$map, run = FALSE,
           verbose = FALSE)))$obj
}


test_that("default_conf() has no cap and check_conf() validates it", {

  sim <- nmax_sim()
  conf <- default_conf(sim$dat, verbose = FALSE)
  expect_true("ctmc_nmax" %in% names(conf))
  expect_null(conf$ctmc_nmax)

  for (bad in list(0, -5, 2.5, c(10, 20), "100", NA_real_)) {
    conf$ctmc_nmax <- bad
    expect_error(check_conf(conf, sim$dat, verbose = FALSE), "ctmc_nmax")
  }
  conf$ctmc_nmax <- 500
  expect_no_error(suppressMessages(suppressWarnings(
    check_conf(conf, sim$dat, verbose = FALSE))))
})


test_that("a cap that does not bind leaves the likelihood unchanged", {

  obj0 <- nmax_obj(NULL)
  obj1 <- nmax_obj(1e4)
  expect_equal(obj1$fn(), obj0$fn(), tolerance = 1e-12)
  expect_equal(obj1$gr(), obj0$gr(), tolerance = 1e-10)

  ## exit rate x step per cell and cached slice; nothing negative
  exit <- obj0$report()$ctmc_exit
  expect_true(length(exit) %% 16 == 0)
  expect_true(all(exit >= 0))
})


test_that("a binding cap truncates the exponential and is flagged after a fit", {

  obj0 <- nmax_obj(NULL)
  rho <- max(obj0$report()$ctmc_exit)
  need <- stats::qpois(admove:::.ctmc_expav_tol, rho, lower.tail = FALSE)
  skip_if(need < 3, "generator too soft to truncate")

  ## truncation loses mass, so the nll gets worse, never better. (A cap of 1
  ## gives Inf here: one term cannot reach every recapture cell.)
  obj1 <- nmax_obj(2)
  expect_gt(obj1$fn(), obj0$fn())

  sim <- nmax_sim()
  sim$conf$ctmc_nmax <- 2
  expect_warning(
    suppressMessages(admove(sim$dat, sim$conf, sim$par, sim$map,
                            do_sdreport = FALSE, do_predictions = FALSE,
                            do_report = FALSE, verbose = FALSE)),
    "binds at the estimates")
})
