## Predicted vs observed tag positions: tag_predictions(), summarise_tag_pred(),
## plot_tag_pred(), plot_tag_resid()

.pred_osa <- function() .cached("pred_osa", tag_predictions(small_fit(report = TRUE)))


test_that("tag_predictions returns one row per predicted observation", {

  fit <- small_fit(report = TRUE)
  pred <- .pred_osa()

  tags <- admove:::.split_tags(fit$dat$tags)

  ## the release of each tag is conditioned on, not predicted
  expect_equal(nrow(pred), sum(vapply(tags, nrow, integer(1))) - length(tags))
  expect_true(all(pred$obs >= 2))
  expect_setequal(unique(pred$id), names(tags))

  expect_true(all(c("pred_x", "sd_x", "res_x", "z_x", "d2", "inside_95",
                    "logdens", "date", "horizon") %in% names(pred)))
  expect_true(all(is.finite(c(pred$pred_x, pred$pred_y, pred$sd_x, pred$sd_y))))
  expect_true(all(pred$sd_x > 0))
  expect_equal(pred$type, rep("osa", nrow(pred)))
})


test_that("the derived residual columns are consistent", {

  pred <- .pred_osa()

  expect_equal(pred$res_x, pred$x - pred$pred_x)
  expect_equal(pred$z_y, (pred$y - pred$pred_y) / pred$sd_y)
  expect_equal(pred$d2, pred$z_x^2 + pred$z_y^2)
  expect_equal(pred$dist, sqrt(pred$res_x^2 + pred$res_y^2))
  expect_equal(pred$inside_95, pred$d2 <= qchisq(0.95, df = 2))
  expect_equal(pred$logdens,
               dnorm(pred$res_x, 0, pred$sd_x, log = TRUE) +
                 dnorm(pred$res_y, 0, pred$sd_y, log = TRUE))

  ## one step ahead: the prediction starts from the previous observation
  d <- pred[pred$id == pred$id[1], ]
  tag <- admove:::.split_tags(small_fit()$dat$tags)[[d$id[1]]]
  expect_equal(d$x_from, as.numeric(tag$x)[d$obs - 1])
  expect_equal(d$horizon, as.numeric(tag$t)[d$obs] - as.numeric(tag$t)[d$obs - 1])
})


test_that("the predictions are the ones the likelihood evaluates", {

  fit <- small_fit(report = TRUE)
  pred <- .pred_osa()

  ## summed log predictive density equals the negative log-likelihood, since
  ## every tag contributes exactly these terms (no ambiguous positions here)
  expect_equal(-sum(pred$logdens), as.numeric(fit$opt$objective), tolerance = 1e-6)
})


test_that("mark-recapture tags are forecasts from release either way", {

  fit <- small_fit(report = TRUE)
  tags <- admove:::.split_tags(fit$dat$tags)
  ctags <- names(tags)[vapply(tags, function(z) z$tag_type[1] == "c", logical(1))]
  ids <- head(ctags, 3)

  osa <- .pred_osa()
  osa <- osa[osa$id %in% ids, ]
  fc <- tag_predictions(fit, type = "forecast", i = ids, verbose = FALSE)

  expect_equal(fc$type, rep("forecast", nrow(fc)))
  expect_equal(fc$pred_x, osa$pred_x)
  expect_equal(fc$sd_y, osa$sd_y)
})


test_that("a forecast of an archival tag ignores the observations", {

  fit <- small_fit(report = TRUE)
  tags <- admove:::.split_tags(fit$dat$tags)
  dtag <- names(tags)[vapply(tags, function(z) z$tag_type[1] == "d", logical(1))][1]

  fc <- tag_predictions(fit, type = "forecast", i = dtag, verbose = FALSE)
  osa <- .pred_osa()
  osa <- osa[osa$id == dtag, ]

  ## uncertainty accumulates from release instead of being reset at every
  ## observation, so it grows and is larger than one step ahead
  expect_true(all(diff(fc$sd_x) > 0))
  expect_gt(tail(fc$sd_x, 1), max(osa$sd_x))
  expect_equal(fc$horizon, fc$t - min(as.numeric(tags[[dtag]]$t)))
})


test_that("tags can be selected by index and by id, and bad input is caught", {

  fit <- small_fit(report = TRUE)
  ids <- names(admove:::.split_tags(fit$dat$tags))

  expect_setequal(unique(tag_predictions(fit, i = 2:3)$id), ids[2:3])
  expect_setequal(unique(tag_predictions(fit, i = ids[2:3])$id), ids[2:3])

  expect_error(tag_predictions(fit, i = 0), "outside")
  expect_error(tag_predictions(fit, i = "nope"), "Unknown tag id")
})


test_that("summarise_tag_pred computes the skill measures", {

  pred <- data.frame(tag_type = c("c", "c", "d", "d"),
                     res_x = c(1, -1, 2, -2),
                     res_y = c(0, 0, 0, 0),
                     sd_x = 1, sd_y = 1,
                     z_x = c(1, -1, 2, -2), z_y = 0,
                     d2 = c(1, 1, 4, 4),
                     inside_50 = c(TRUE, TRUE, FALSE, FALSE),
                     inside_95 = c(TRUE, TRUE, TRUE, TRUE),
                     logdens = c(-1, -1, -3, -3),
                     stringsAsFactors = FALSE)

  s <- summarise_tag_pred(m = pred)

  expect_equal(s$model, c("m", "m"))
  expect_equal(s$tag_type, c("c", "d"))
  expect_equal(s$n, c(2L, 2L))
  expect_equal(s$rmse, c(1, 2))
  expect_equal(s$cover_50, c(1, 0))
  expect_equal(s$logdens, c(-1, -3))

  ## grouping can be switched off, and several tables compared
  s2 <- summarise_tag_pred(a = pred, b = pred, by = NULL)
  expect_equal(s2$model, c("a", "b"))
  expect_equal(s2$n, c(4L, 4L))

  expect_error(summarise_tag_pred(), "No prediction tables")
  expect_error(summarise_tag_pred(data.frame(x = 1)), "neither an 'admove' fit")
  expect_error(summarise_tag_pred(m = pred, by = "nope"), "is missing")
})


test_that("the prediction plots draw without error", {

  fit <- small_fit(report = TRUE)
  pred <- .pred_osa()
  tags <- admove:::.split_tags(fit$dat$tags)
  dtag <- names(tags)[vapply(tags, function(z) z$tag_type[1] == "d", logical(1))][1]

  pdf(NULL)
  on.exit(dev.off())

  ## single tag: map plus coordinates over time
  expect_s3_class(plot_tag_pred(fit, i = dtag, type = "osa", pred = pred,
                                plot_land = FALSE),
                  "data.frame")
  ## several tags: one map
  expect_s3_class(plot_tag_pred(fit, i = 2:5, type = "osa", pred = pred,
                                plot_land = FALSE),
                  "data.frame")
  expect_s3_class(plot_tag_resid(fit, pred = pred), "data.frame")
  expect_s3_class(plot_tag_resid(fit, pred = pred, tag_type = "c"), "data.frame")
  expect_s3_class(plot_tag_resid(fit, pred = pred, map_at = "from"), "data.frame")

  ## the pairing cues: dated pairs only, every pair, or none
  for (lk in list(TRUE, "all", FALSE)) {
    expect_s3_class(plot_tag_pred(fit, i = dtag, type = "osa", pred = pred,
                                  link = lk, plot_land = FALSE),
                    "data.frame")
  }

  expect_error(plot_tag_resid(fit, pred = pred, tag_type = "s"),
               "No predicted positions left")
})


test_that("tag_predictions runs the likelihood's own CTMC pass", {

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
  ll_fit <- fit$obj$report(fit$opt$par)$loglik_tags
  ids <- names(admove:::.split_tags(fit$dat$tags))
  area <- prod(fit$dat$grid$cellsize)

  ## "osa" is what the likelihood scores: the rows' likelihood terms (logdens
  ## is per unit area) add up to every tag's term of the fit, grouped
  ## mark-recapture tags too, up to the rounding of their shared pass
  osa <- tag_predictions(fit, type = "osa")

  ## quantile residuals: reproducible by seed, randomised for these exact
  ## positions, and the caller's random stream is left alone
  set.seed(42)
  r_before <- runif(1)
  set.seed(42)
  expect_identical(tag_predictions(fit, type = "osa")$z_x, osa$z_x)
  expect_identical(runif(1), r_before)
  expect_false(identical(tag_predictions(fit, type = "osa", seed = 2)$z_x,
                         osa$z_x))
  expect_true(all(is.finite(osa$z_x) & is.finite(osa$z_y)))

  ll_osa <- tapply(osa$logdens + log(area), factor(osa$id, levels = ids), sum)
  expect_equal(unname(as.numeric(ll_osa)), ll_fit, tolerance = 1e-8)

  ## the distributions: one column per tag row, each summing to 1, and the
  ## predicted mean is their mean over the cell centres
  info <- attr(osa, "ctmc")
  expect_named(info$dist, ids, ignore.order = TRUE)
  d1 <- info$dist[[ids[1]]]
  expect_equal(ncol(d1$prob), nrow(d1$tag))
  expect_equal(unname(colSums(d1$prob)), rep(1, nrow(d1$tag)))
  r1 <- osa[osa$id == ids[1], ]
  expect_equal(r1$pred_x, colSums(d1$prob[, r1$obs, drop = FALSE] *
                                    fit$dat$grid$xygrid[, 1]))
  expect_true(all(r1$sd_x > 0))

  ## "forecast" switches the updates off: identical for tags that are never
  ## updated, different for an archival tag after its first update
  fc <- tag_predictions(fit, type = "forecast")
  types <- vapply(fc_d <- attr(fc, "ctmc")$dist,
                  function(d) as.character(d$tag$tag_type[1]), character(1))
  k_c <- names(types)[types == "c"][1]
  k_d <- names(types)[types == "d" &
                        vapply(fc_d, function(d) nrow(d$tag), 1L) > 2][1]
  expect_equal(fc_d[[k_c]]$prob, info$dist[[k_c]]$prob)
  p_fc <- fc_d[[k_d]]$prob
  p_osa <- info$dist[[k_d]]$prob
  expect_gt(max(abs(p_fc[, ncol(p_fc)] - p_osa[, ncol(p_osa)])), 1e-3)

  ## selection by id and by index
  expect_setequal(unique(tag_predictions(fit, i = ids[2:3])$id), ids[2:3])
  expect_named(attr(tag_predictions(fit, i = 1), "ctmc")$dist, ids[1])

  ## independent of the prediction times: a window that covers none of the
  ## tags used to fail with "NA/NaN argument"
  fit_w <- fit
  fit_w$dat$pred$time <- max(fit$dat$tags$t) + 1:3
  fit_w$pred$mstar <- NULL
  expect_equal(attr(tag_predictions(fit_w, i = 1), "ctmc")$dist[[1]]$prob,
               info$dist[[ids[1]]]$prob)

  pdf(NULL)
  on.exit(dev.off())

  ## maps for one tag and for several, the mean map on request, a supplied
  ## table (whose maps survive the subsetting), and the residual panels
  expect_silent(plot_tag_pred(fit, i = ids[1], plot_land = FALSE))
  expect_silent(plot_tag_pred(fit, i = 1:3, type = "osa", dist = TRUE,
                              plot_land = FALSE))
  expect_silent(plot_tag_pred(fit, i = 1:3, type = "osa", dist = FALSE,
                              plot_land = FALSE, plot_grid = FALSE))
  expect_silent(plot_tag_pred(fit, i = ids[k_d == names(fc_d)], pred = fc,
                              plot_land = FALSE))
  expect_silent(plot_tag_resid(fit, pred = osa))
})


test_that("CTMC quantile residuals are N(0, 1) under the predicted distribution", {

  grid <- create_grid(xrange = c(0, 1), yrange = c(0, 1), cellsize = 0.2,
                      verbose = FALSE)
  nc <- nrow(grid$xygrid)
  ix <- grid$igrid$idx
  iy <- grid$igrid$idy

  ## skewed and bimodal: most mass in one corner, a second mode opposite, and
  ## one cell without mass (e.g. unreachable)
  p <- exp(-3 * (ix + iy)) + 0.3 * exp(-2 * ((ix - 5)^2 + (iy - 4)^2))
  p[13] <- 0
  p <- p / sum(p)

  n <- 3000
  z <- withr::with_seed(1, {
    ic <- sample.int(nc, n, replace = TRUE, prob = p)
    x <- runif(n, grid$xgr[ix[ic]], grid$xgr[ix[ic] + 1L])
    y <- runif(n, grid$ygr[iy[ic]], grid$ygr[iy[ic] + 1L])

    ## exact positions: only the cell counts, randomised within it
    z0 <- t(vapply(seq_len(n), function(k)
      admove:::.ctmc_quantile_resid(p, grid, ic[k], x[k], y[k]), numeric(2)))

    ## with observation error: continuous, uniform in the cell plus N(0, sd^2)
    sd <- c(0.15, 0.1)
    xo <- x + rnorm(n, 0, sd[1])
    yo <- y + rnorm(n, 0, sd[2])
    z1 <- t(vapply(seq_len(n), function(k)
      admove:::.ctmc_quantile_resid(p, grid, ic[k], xo[k], yo[k], sd),
      numeric(2)))
    list(z0 = z0, z1 = z1)
  })

  for (zz in z) {
    expect_true(all(is.finite(zz)))
    expect_gt(suppressWarnings(ks.test(zz[, 1], "pnorm"))$p.value, 0.01)
    expect_gt(suppressWarnings(ks.test(zz[, 2], "pnorm"))$p.value, 0.01)
    expect_lt(abs(cor(zz[, 1], zz[, 2])), 0.06)
  }
})


test_that("plot_tag_resid handles many residuals, mixed tag types and NA", {

  n <- 6000
  pred <- withr::with_seed(1, data.frame(
    tag_type = sample(c("d", "c"), n, replace = TRUE),
    t = runif(n, 0, 10), date = as.POSIXct(NA), horizon = runif(n, 0, 5),
    x = runif(n), y = runif(n), pred_x = runif(n), pred_y = runif(n),
    z_x = rnorm(n), z_y = rnorm(n), stringsAsFactors = FALSE))
  pred$z_x[1] <- NA

  pdf(NULL)
  on.exit(dev.off())

  ## more than 5000 residuals: Shapiro-Wilk on a subsample, without touching
  ## the caller's random numbers
  set.seed(7)
  r <- runif(1)
  set.seed(7)
  expect_message(out <- plot_tag_resid(NULL, pred = pred), "not shown")
  expect_identical(runif(1), r)
  expect_equal(nrow(out), n - 1L)
  expect_s3_class(plot_tag_resid(NULL, pred = pred[-1, ], tag_type = "c"),
                  "data.frame")
  expect_s3_class(plot_tag_resid(NULL, pred = pred[-1, ], map_at = "obs"),
                  "data.frame")
  expect_error(plot_tag_resid(NULL, pred = pred[-1, ], map_at = "from"),
               "x_from")
})
