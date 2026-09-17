## Number of default spline knots in setup_data(), sim_data() and sim_tags()

make_knot_cov <- function(ncov = 1) {
  grid <- create_grid(cellsize = 0.25, verbose = FALSE)
  cov <- lapply(seq_len(ncov), function(i) suppressMessages(sim_cov(grid, nt = 2)))
  names(cov) <- paste0("cov", seq_len(ncov))
  list(grid = grid, cov = cov)
}

knot_dat <- function(...) {
  inp <- make_knot_cov(2)
  suppressWarnings(suppressMessages(
    setup_data(grid = inp$grid, cov = inp$cov, trange = c(0, 1),
               verbose = FALSE, ...)
  ))
}


test_that("setup_data keeps 3 taxis knots and 1 diffusion knot by default", {

  dat <- knot_dat()

  expect_equal(dim(dat$knots_tax), c(3L, 2L))
  expect_equal(dim(dat$knots_dif), c(1L, 2L))
  expect_equal(dat$knots_tax[, 1],
               as.numeric(quantile(as.numeric(dat$cov[[1]]),
                                   c(0.05, 0.5, 0.95), na.rm = TRUE)))
})


test_that("n_knots_tax and n_knots_dif set the number of default knots", {

  dat <- knot_dat(n_knots_tax = 6, n_knots_dif = 2)

  expect_equal(dim(dat$knots_tax), c(6L, 2L))
  expect_equal(dim(dat$knots_dif), c(2L, 2L))
  expect_equal(dat$knots_tax[, 2],
               as.numeric(quantile(as.numeric(dat$cov[[2]]),
                                   seq(0.05, 0.95, length.out = 6),
                                   na.rm = TRUE)))
  expect_equal(dat$knots_dif[, 1],
               as.numeric(quantile(as.numeric(dat$cov[[1]]),
                                   c(0.25, 0.75), na.rm = TRUE)))

  ## the parameters and the map follow the knots
  conf <- default_conf(dat, verbose = FALSE)
  par <- suppressMessages(default_par(dat, conf, verbose = FALSE))
  expect_equal(dim(par$alpha)[1:2], c(6L, 2L))
  expect_equal(dim(par$beta)[1:2], c(2L, 2L))
  expect_silent(admove:::.check_knots_dims(dat, par))
})


test_that("a knot matrix wins over the number of knots", {

  inp <- make_knot_cov(1)
  knots <- matrix(c(21, 24, 27), 3, 1)

  dat <- suppressWarnings(suppressMessages(
    setup_data(grid = inp$grid, cov = inp$cov, trange = c(0, 1),
               knots_tax = knots, n_knots_tax = 5, verbose = FALSE)
  ))

  expect_equal(dat$knots_tax, knots)
})


test_that("invalid numbers of knots are rejected", {

  inp <- make_knot_cov(1)
  for (n in list(0, 2.5, c(2, 3), NA_real_, "3")) {
    expect_error(
      setup_data(grid = inp$grid, cov = inp$cov, trange = c(0, 1),
                 n_knots_tax = n, verbose = FALSE),
      "n_knots_tax"
    )
  }
})


test_that("coinciding default knots warn and name the covariate", {

  inp <- make_knot_cov(1)
  cov <- inp$cov
  ## mostly one value, so several quantiles coincide
  vals <- as.numeric(cov[[1]])
  vals[] <- 20
  vals[seq(1, length(vals), by = 10)] <- 25
  cov[[1]][] <- vals

  expect_warning(
    suppressMessages(setup_data(grid = inp$grid, cov = cov, trange = c(0, 1),
                                n_knots_tax = 4, verbose = FALSE)),
    "covariate 'cov1'"
  )
})


test_that("sim_data and sim_tags take the number of knots", {

  inp <- make_knot_cov(1)

  set.seed(1)
  sim <- suppressWarnings(suppressMessages(
    sim_data(grid = inp$grid, cov = inp$cov, n_knots_tax = 4, n_knots_dif = 2,
             n_ctags = 5, n_dtags = 1, verbose = FALSE)
  ))
  expect_equal(nrow(sim$dat$knots_tax), 4L)
  expect_equal(nrow(sim$dat$knots_dif), 2L)
  expect_equal(dim(sim$par_true$alpha)[1], 4L)
  expect_equal(dim(sim$par_true$beta)[1], 2L)
  expect_equal(dim(sim$par$alpha)[1], 4L)

  ## the likelihood can be built and evaluated with these knots
  obj <- suppressWarnings(suppressMessages(
    admove(sim, run = FALSE, verbose = FALSE)
  ))$obj
  expect_true(is.finite(obj$fn(obj$par)))

  ## a different number replaces the knots inherited from a simulation, and
  ## the parameters sized for the old knots with them
  set.seed(2)
  sim5 <- suppressWarnings(suppressMessages(
    sim_data(sim, n_knots_tax = 5, n_ctags = 5, n_dtags = 1, verbose = FALSE)
  ))
  expect_equal(nrow(sim5$dat$knots_tax), 5L)
  expect_equal(nrow(sim5$dat$knots_dif), 2L)
  expect_equal(dim(sim5$par_true$alpha)[1], 5L)

  ## the same number keeps the inherited knots
  set.seed(3)
  sim4 <- suppressWarnings(suppressMessages(
    sim_data(sim, n_knots_tax = 4, n_ctags = 5, n_dtags = 1, verbose = FALSE)
  ))
  expect_equal(sim4$dat$knots_tax, sim$dat$knots_tax)
  expect_equal(sim4$par_true$alpha, sim$par_true$alpha)

  set.seed(4)
  tags <- suppressWarnings(suppressMessages(
    sim_tags("c", grid = inp$grid, cov = inp$cov, n_tags = 5,
             n_knots_tax = 5, verbose = FALSE)
  ))
  expect_equal(nrow(tags$dat$knots_tax), 5L)
  expect_equal(dim(tags$par_true$alpha)[1], 5L)
})
