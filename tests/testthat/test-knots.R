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


## Default knots from the covariate values at the tag observations

tag_knot_dat <- function(tags = NULL, ...) {
  d <- small_sim()$dat
  if (is.null(tags)) tags <- d$tags
  suppressWarnings(suppressMessages(
    setup_data(grid = d$grid, cov = d$cov, tags = tags, trange = d$trange,
               shift_tref = TRUE, verbose = FALSE, ...)
  ))
}

wq <- function(x, w, p) admove:::.wquantile(x, w, p)


test_that(".wquantile equals quantile() for equal weights", {

  set.seed(2)
  x <- rnorm(37)
  p <- c(0, 0.05, 0.3, 0.5, 0.95, 1)
  expect_equal(wq(x, NULL, p), unname(quantile(x, p)))
  expect_equal(wq(x, rep(3, 37), p), unname(quantile(x, p)))
  ## a heavier value pulls the quantiles towards it
  expect_gt(wq(c(1, 2, 3), c(1, 1, 5), 0.5), 2)
  expect_lt(wq(c(1, 2, 3), c(5, 1, 1), 0.5), 2)
})


test_that("default knots are quantiles of the covariate at the tags, per tag", {

  dat <- tag_knot_dat()
  vals <- cov_at_tags(dat)
  w <- admove:::.tag_weights(dat$tags)

  ## every tag has total weight one
  expect_equal(as.numeric(tapply(w, dat$tags$id, sum)),
               rep(1, length(unique(dat$tags$id))))
  expect_equal(dat$knots_tax[, 1], wq(vals[[1]], w, c(0.05, 0.5, 0.95)))
  expect_equal(dat$knots_dif[, 1], wq(vals[[1]], w, 0.5))
  expect_equal(dat$knots_from, c(tax = "tags", dif = "tags"))

  ## within the range the tags experienced, unlike the field quantiles
  expect_true(all(dat$knots_tax >= min(vals[[1]], na.rm = TRUE) &
                    dat$knots_tax <= max(vals[[1]], na.rm = TRUE)))

  ## the likelihood can be built with these knots
  conf <- default_conf(dat, verbose = FALSE)
  par <- suppressMessages(default_par(dat, conf, verbose = FALSE))
  expect_silent(admove:::.check_knots_dims(dat, par))
})


test_that("with only conventional tags, per-tag weights give plain quantiles", {

  d <- small_sim()$dat
  ctags <- d$tags[d$tags$tag_type == "c", ]
  dat <- tag_knot_dat(tags = ctags)
  vals <- cov_at_tags(dat)[[1]]

  expect_equal(dat$knots_tax[, 1],
               unname(quantile(vals, c(0.05, 0.5, 0.95), na.rm = TRUE)))
})


test_that("knots_from = 'cov' and missing tags give field quantiles", {

  d <- small_sim()$dat
  field <- unname(quantile(as.numeric(d$cov[[1]]), c(0.05, 0.5, 0.95),
                           na.rm = TRUE))

  dat <- tag_knot_dat(knots_from = "cov")
  expect_equal(dat$knots_tax[, 1], field)
  expect_equal(dat$knots_from, c(tax = "cov", dif = "cov"))

  dat0 <- suppressWarnings(suppressMessages(
    setup_data(grid = d$grid, cov = d$cov, trange = d$trange, verbose = FALSE)
  ))
  expect_equal(dat0$knots_tax[, 1], field)

  ## a supplied knot matrix wins in either mode
  knots <- matrix(c(22, 24, 26), 3, 1)
  expect_equal(tag_knot_dat(knots_tax = knots)$knots_tax, knots)
  expect_equal(tag_knot_dat(knots_tax = knots)$knots_from[["tax"]], "user")
})


test_that("cov_at_tags matches the tag rows and accepts sim and data objects", {

  sim <- small_sim()
  v <- cov_at_tags(sim)
  expect_named(v, names(sim$dat$cov))
  expect_length(v[[1]], nrow(sim$dat$tags))
  expect_null(attr(v[[1]], "slice"))
  expect_equal(v, cov_at_tags(sim$dat))
  expect_error(cov_at_tags(list()), "admove_data")
})


test_that("candidate positions of an ambiguous recapture share one weight", {

  tags <- data.frame(id = c("a", "a", "a", "b", "b"),
                     event = c(1, 2, 2, 1, 2),
                     prob = c(1, 0.25, 0.75, 1, 1))
  expect_equal(admove:::.tag_weights(tags),
               c(0.5, 0.125, 0.375, 0.5, 0.5))
})
