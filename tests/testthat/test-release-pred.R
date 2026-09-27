test_that("release_predictions steps like the CTMC likelihood", {

  grid <- create_grid(xrange = c(0, 1), yrange = c(0, 1), cellsize = 0.25,
                      verbose = FALSE)
  sim <- withr::with_seed(1, suppressMessages(suppressWarnings(
    sim_data(grid = grid, n_dtags = 2, n_ctags = 5, trange = c(0, 1),
             verbose = FALSE)
  )))
  conf <- sim$conf
  conf$engine <- "ctmc"
  fit <- suppressWarnings(suppressMessages(
    admove(sim$dat, conf, sim$par, sim$map, verbose = FALSE)
  ))

  ## a release at an archival tag's release, predicted at the tag's times,
  ## reproduces the tag's forecast from tag_predictions(): same lattice, same
  ## generator, same exponential
  tags <- split(fit$dat$tags, fit$dat$tags$id)
  k <- which(vapply(tags, function(z) z$tag_type[1] == "d" && nrow(z) > 2,
                    logical(1)))[1]
  tg <- tags[[k]]
  fc <- tag_predictions(fit, type = "forecast", i = names(tags)[k])
  prob <- attr(fc, "ctmc")$dist[[1]]$prob
  tt <- unique(tg$t)
  rp <- release_predictions(fit, x0 = tg$x[1], y0 = tg$y[1], t0 = tg$t[1],
                            times = tt[-1] - tg$t[1])
  rows <- match(tt[-1], tg$t)
  expect_equal(rp$releases[[1]]$prob, prob[, rows, drop = FALSE])

  ## several releases, recycled arguments, probabilities over the grid
  rp2 <- release_predictions(fit, x0 = c(0.3, 0.7), y0 = 0.5, t0 = 0.1,
                             times = c(0.1, 0.3), n = c(100, 10))
  expect_s3_class(rp2, "admove_release_pred")
  expect_length(rp2$releases, 2L)
  for (z in rp2$releases) {
    expect_equal(dim(z$prob), c(nrow(fit$dat$grid$xygrid), 2L))
    expect_equal(unname(colSums(z$prob)), c(1, 1))
    expect_equal(z$t, 0.1 + c(0.1, 0.3))
  }
  expect_equal(rp2$releases[[2]]$n, 10)
  expect_output(print(rp2), "releases: 2")

  pdf(NULL)
  on.exit(dev.off())
  expect_silent(plot_release_pred(rp2))
  expect_silent(plot_release_pred(rp2, type = "count", select = 2,
                                  select_release = 1))
  expect_error(plot_release_pred(rp2, select = 3), "indices of the prediction")

  expect_error(release_predictions(fit, x0 = 5, y0 = 5, t0 = 0.1),
               "not in a grid cell")
  expect_error(release_predictions(fit, x0 = 0.5, y0 = 0.5, t0 = 0.1,
                                   times = c(0, 1)), "positive times")
})


test_that("release_predictions works for a Kalman filter fit", {

  fit <- small_fit()
  xy <- fit$dat$grid$xygrid[1, ]
  rp <- release_predictions(fit, x0 = xy[1], y0 = xy[2],
                            t0 = min(fit$dat$tags$t), times = 1:2)
  expect_equal(unname(colSums(rp$releases[[1]]$prob)), c(1, 1))
})
