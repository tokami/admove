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


test_that("tag_predictions refuses the CTMC engine", {

  fit <- small_fit(report = TRUE)
  fit$conf$engine <- 2L

  expect_error(tag_predictions(fit), "Kalman filter")
  expect_error(add_tag_dist(fit_kf <- within(fit, conf$engine <- 1L)),
               "CTMC engine")
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

  ## the pairing cues: dated pairs only, every pair, or none
  for (lk in list(TRUE, "all", FALSE)) {
    expect_s3_class(plot_tag_pred(fit, i = dtag, type = "osa", pred = pred,
                                  link = lk, plot_land = FALSE),
                    "data.frame")
  }

  expect_error(plot_tag_resid(fit, pred = pred, tag_type = "s"),
               "No predicted positions left")
})


test_that("add_tag_dist and plot_tag_dist accept tag ids (CTMC)", {

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

  ## the name survives the fit, while nll() saw the integer code
  expect_equal(fit$conf$engine, "ctmc")

  ids <- names(admove:::.split_tags(fit$dat$tags))
  fit <- suppressWarnings(suppressMessages(add_tag_dist(fit, i = ids[1])))

  expect_named(fit$tag_dist, "1")
  expect_equal(as.character(fit$tag_dist[[1]]$tag$id[1]), ids[1])

  pdf(NULL)
  on.exit(dev.off())

  ## selection by id, by index, and the alias for 'select'
  expect_silent(plot_tag_dist(fit, i = ids[1]))
  expect_silent(plot_tag_dist(fit, select = 1))
  expect_error(plot_tag_dist(fit, select = 1, i = ids[1]), "not both")
  expect_error(plot_tag_dist(fit, i = "nope"), "no precomputed distribution")
})
