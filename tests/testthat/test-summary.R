## require(admove); require(testthat)


test_that(".select_estimated_par drops fixed elements", {

  pl <- list(alpha = array(0, c(3, 1, 1)), logKappa = 0)
  map <- list(alpha = factor(c(NA, 1, 2)), logKappa = factor(NA))

  sel <- admove:::.select_estimated_par(pl, map)

  expect_equal(sel$keep, c(FALSE, TRUE, TRUE, FALSE))
  expect_equal(sel$labels[sel$keep], c("alpha2", "alpha3"))
})


test_that(".select_estimated_par shows coupled elements once", {

  pl <- list(gamma = array(0, c(2, 3, 1)))
  map <- list(gamma = factor(c(NA, NA, 1, NA, NA, 1)))

  sel <- admove:::.select_estimated_par(pl, map)

  expect_equal(sum(sel$keep), 1L)
  expect_equal(sel$labels[sel$keep], "gamma3,6")

  ## the kept element is the first of the level, so it carries the estimate
  expect_equal(which(sel$keep), 3L)
})


test_that(".select_estimated_par keeps unmapped parameters", {

  pl <- list(alpha = array(0, c(2, 1, 1)), logKappa = 0)
  map <- list(alpha = factor(c(NA, 1)))

  sel <- admove:::.select_estimated_par(pl, map)

  expect_equal(sel$keep, c(FALSE, TRUE, TRUE))
  expect_equal(sel$labels[sel$keep], c("alpha2", "logKappa"))
})


test_that(".select_estimated_par drops empty map entries", {

  pl <- list(gamma = array(0, c(2, 1, 1)))
  map <- list(gamma = factor(character(0)))

  sel <- admove:::.select_estimated_par(pl, map)

  expect_false(any(sel$keep))
})


test_that(".par_display_labels maps element names to summary row names", {

  fit <- list(pl = list(gamma = array(0, c(2, 3, 1))),
              map = list(gamma = factor(c(NA, NA, 1, NA, NA, 1))))
  class(fit) <- "admove"

  expect_equal(admove:::.par_display_labels(fit, "gamma3"), "gamma3,6")

  ## keys without a label are passed through unchanged
  expect_equal(admove:::.par_display_labels(fit, c("gamma3", "beta1")),
               c("gamma3,6", "beta1"))
})


test_that("summary prints one row per estimated parameter", {

  sim <- skjepo$sim
  conf <- sim$conf
  conf$use_advection <- TRUE
  map <- default_map(sim$dat, conf, sim$par)

  ## default_map couples the x- and y-direction advection coefficients
  expect_equal(nlevels(map$gamma), 1L)
  expect_equal(sum(!is.na(map$gamma)), 2L)

  fit <- suppressWarnings(
    admove(sim, conf = conf, map = map, do_predictions = FALSE,
           do_report = FALSE, verbose = FALSE)
  )

  out <- capture.output(summary(fit))

  ## one row per entry in opt$par, i.e. per logLik degree of freedom
  expect_equal(sum(grepl("^ (alpha|beta|gamma|logSdO)", out)),
               length(fit$opt$par))
  expect_equal(sum(grepl("gamma", out)), 1L)
  expect_true(any(grepl("gamma1,2", out, fixed = TRUE)))

  ## the parameter plot labels its x-axis with the same names
  sel <- admove:::.select_estimated_par(fit$pl, fit$map)
  keys <- names(unlist(fit$pl))[sel$keep]

  expect_equal(admove:::.par_display_labels(fit, keys), sel$labels[sel$keep])
  expect_true(all(vapply(sel$labels[sel$keep],
                         function(l) any(grepl(l, out, fixed = TRUE)),
                         logical(1))))
})
