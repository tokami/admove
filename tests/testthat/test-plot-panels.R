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
