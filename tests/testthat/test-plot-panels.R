## require(admove); require(testthat)


## Minimal stand-in carrying only what plot_fit() inspects before it starts
## drawing: the covariates, the parameter arrays, the map and the config.
make_panel_fit <- function(use_advection = FALSE, gamma = 0) {

  x <- list(
    dat = list(cov = list(a = NULL, b = NULL)),
    conf = list(use_advection = use_advection),
    par = list(alpha = array(0, c(3, 2, 2)),
               beta = array(0, c(1, 2, 2)),
               gamma = array(gamma, c(2, 2, 2))),
    map = list(alpha = factor(rep(c(NA, 1, 2), 4)),
               beta = factor(c(1, NA, 1, NA)),
               gamma = factor(rep(NA, 8)))
  )

  admove:::.add_class(x, "admove")
}


test_that(".adv_active is FALSE without advection or with all-zero gamma", {

  expect_false(admove:::.adv_active(make_panel_fit(use_advection = FALSE, gamma = 0)))
  expect_false(admove:::.adv_active(make_panel_fit(use_advection = TRUE, gamma = 0)))

  ## fitted with advection off but non-zero gamma left in par: still no field
  expect_false(admove:::.adv_active(make_panel_fit(use_advection = FALSE, gamma = 2)))
})


test_that(".adv_active is TRUE for fixed but non-zero gamma", {

  ## gamma is entirely mapped NA, so it is fixed rather than estimated -- it
  ## still produces a real advection field and must not be dropped
  fit <- make_panel_fit(use_advection = TRUE, gamma = 1.5)

  expect_true(all(is.na(fit$map$gamma)))
  expect_true(admove:::.adv_active(fit))
})


test_that(".adv_active prefers the estimates over the starting values", {

  fit <- make_panel_fit(use_advection = TRUE, gamma = 0)
  fit$pl <- list(gamma = array(0.3, c(2, 2, 2)))

  expect_true(admove:::.adv_active(fit))
})


test_that("plot_fit drops the advection panels for a model without advection", {

  fit <- make_panel_fit(use_advection = FALSE, gamma = 0)

  expect_warning(res <- plot_fit(fit, quantity = "advection"),
                 "Nothing to plot")
  expect_null(res)
})


test_that(".is_constant_field detects fields without spatial variation", {

  expect_true(admove:::.is_constant_field(rep(5.4, 100)))
  expect_true(admove:::.is_constant_field(rep(0, 100)))
  expect_true(admove:::.is_constant_field(1))
  expect_true(admove:::.is_constant_field(numeric(0)))

  expect_false(admove:::.is_constant_field(c(1, 2, 3)))
  expect_false(admove:::.is_constant_field(c(0, 0, 1e-6)))
})


test_that(".is_constant_field tolerates rounding noise and ignores non-finite values", {

  z <- rep(239.36, 500) + c(0, rep(1e-12, 499))
  expect_true(admove:::.is_constant_field(z))

  expect_true(admove:::.is_constant_field(c(3, 3, NA, NaN, Inf)))
  expect_false(admove:::.is_constant_field(c(3, 4, NA)))
})


## The preference plots mark the knots on the curve. Because the smooth
## interpolates its knots, coefficient i IS the curve's value at knot i, so the
## markers must land exactly on the line. They did not when do_sdreport was
## FALSE: the plotting code reshaped opt$par by hand, prepending a zero row on
## the assumption that the first coefficient is fixed. That holds for alpha but
## not for beta -- .make_beta_map() frees the intercept of one covariate and
## fixes the others -- so the diffusion knot was drawn at 0 instead of its
## fitted value, and it fails for any hand-edited map too.

test_that(".fitted_par returns the stored parameter list when present", {

  fit <- small_fit()

  expect_false(is.null(fit$pl))
  expect_identical(admove:::.fitted_par(fit), fit$pl)
})


test_that(".fitted_par rebuilds the parameter list when pl is absent", {

  fit <- small_fit()
  pl <- fit$pl
  fit$pl <- NULL

  expect_equal(admove:::.fitted_par(fit), pl)
})


test_that("preference knots lie on the preference curve", {

  for (sdrep in c(FALSE, TRUE)) {

    fit <- small_fit(sdreport = sdrep, report = TRUE)
    pe <- admove:::.fitted_par(fit)

    for (ty in c("taxis", "diffusion")) {
      co <- if (ty == "taxis") pe$alpha[, 1, 1] else pe$beta[, 1, 1]
      kn <- if (ty == "taxis") fit$dat$knots_tax[, 1] else fit$dat$knots_dif[, 1]
      f <- admove:::.poly_fun(kn, co, method = fit$conf$smooth_method)

      expect_equal(as.numeric(sapply(kn, f)), as.numeric(co),
                   tolerance = 1e-8,
                   info = paste(ty, "with do_sdreport =", sdrep))
    }

    ## the diffusion intercept is a free parameter, not a fixed zero: the old
    ## reconstruction assumed the opposite and drew the marker at the origin
    expect_true(any(abs(pe$beta) > 1e-6))
  }
})


test_that("plot_tag_dist reports an empty tag_dist instead of failing on mfrow", {

  ## add_tag_dist() warns and skips a tag it cannot build a prediction for, so
  ## fit$tag_dist can be an empty list rather than NULL. That used to reach
  ## par(mfrow = c(0, -Inf)) via max(integer(0)) and die with an opaque
  ## "invalid value specified for graphical parameter" message.
  fit <- list(tag_dist = list())
  class(fit) <- "admove"
  expect_error(plot_tag_dist(fit), "skipped every tag")

  ## an absent store keeps its own, different message
  fit$tag_dist <- NULL
  expect_error(plot_tag_dist(fit), "Run add_tag_dist\\(\\) first")

  ## selecting nothing out of a populated store is also reported, not drawn
  fit$tag_dist <- list("1" = list(engine = 1L, i = 1L,
                                  tag = data.frame(t = c(0, 1), x = c(0, 1),
                                                   y = c(0, 1))))
  expect_error(plot_tag_dist(fit, n_tags = 0), "No tags left to plot")
})


test_that("plot_data reserves one panel per covariate and per tag type", {

  dat <- skjepo$sim$dat
  dat$cov <- suppressMessages(make_x_y_cov(dat$grid))
  expect_setequal(unique(dat$tags$tag_type), c("d", "c"))

  ## grid + 2 covariates + d + c tags = 5 panels; too few panels spill onto a
  ## second page
  dir <- withr::local_tempdir()
  grDevices::png(file.path(dir, "p%03d.png"))
  plot(dat)
  grDevices::dev.off()
  expect_length(list.files(dir, pattern = "\\.png$"), 1L)
})


test_that("plot_sim reserves exactly as many panels as it draws", {

  sim <- skjepo$sim
  n_types <- length(unique(sim$tags$tag_type))
  ncov <- if (inherits(sim$cov, "list")) length(sim$cov) else 1L

  n_seen <- NULL
  local_mocked_bindings(
    n2mfrow = function(nr.plots, ...) {
      ## the first call is plot_sim's own layout; lower-level plots may call it too
      if (is.null(n_seen)) n_seen <<- nr.plots
      grDevices::n2mfrow(nr.plots, ...)
    },
    .package = "admove"
  )

  .count <- function(...) {
    n_seen <<- NULL
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off())
    plot(sim, ...)
    n_seen
  }

  expect_equal(.count(by_tag_type = TRUE), 2 + 2 * ncov + n_types)
  expect_equal(.count(by_tag_type = FALSE), 2 + 2 * ncov + 1)

  sim$tags <- NULL
  expect_equal(.count(), 2 + 2 * ncov)
})
