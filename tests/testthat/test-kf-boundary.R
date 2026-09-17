## KF: predicted positions that leave the covariate field (conf$kf_boundary)

test_that(".clamp_ad clamps values and is the identity inside the bounds", {

  v <- c(-3, -1, 0, 0.5, 1, 2)
  expect_equal(admove:::.clamp_ad(v, -1, 1), c(-1, -1, 0, 0.5, 1, 1))

  ## on the tape: derivative 1 inside, 0 outside
  f <- RTMB::MakeTape(function(p) admove:::.clamp_ad(p, -1, 1), c(0.5, 2, -3))
  expect_equal(as.numeric(f(c(0.5, 2, -3))), c(0.5, 1, -1))
  expect_equal(diag(f$jacobian(c(0.5, 2, -3))), c(1, 0, 0))
})


test_that("check_conf validates kf_boundary", {

  dat <- small_sim()$dat
  expect_equal(default_conf(dat, verbose = FALSE)$kf_boundary, "clamp")
  expect_error(check_conf(list(kf_boundary = "reflect"), dat, verbose = FALSE),
               "kf_boundary")
})


## One tape per setting, built at taxis coefficients 10 times the simulated
## ones, which push some predicted tracks out of the covariate field.
.boundary_objs <- function() {
  .cached("boundary_objs", {
    sim <- small_sim()
    par <- sim$par
    par$alpha[] <- 10 * sim$par_true$alpha
    lapply(c(none = "none", clamp = "clamp"), function(b) {
      conf <- sim$conf
      conf$kf_boundary <- b
      suppressWarnings(suppressMessages(
        admove(sim$dat, conf, par, sim$map, run = FALSE, verbose = FALSE)
      ))$obj
    })
  })
}


test_that("clamping does not change the likelihood while tracks stay inside", {

  objs <- .boundary_objs()
  sim <- small_sim()

  ## the simulated taxis keeps all predicted tracks inside the field
  p <- objs$clamp$par
  p[names(p) == "alpha"] <- p[names(p) == "alpha"] / 10

  expect_true(is.finite(objs$none$fn(p)))
  expect_equal(objs$clamp$fn(p), objs$none$fn(p))
  expect_equal(objs$clamp$gr(p), objs$none$gr(p))

  ids <- names(split(sim$dat$tags, sim$dat$tags$id))
  expect_length(admove:::.boundary_tags(objs$clamp, p, ids), 0)
})


test_that("clamping keeps the likelihood finite when tracks leave the field", {

  objs <- .boundary_objs()
  sim <- small_sim()
  p <- objs$clamp$par

  expect_true(is.nan(objs$none$fn(p)))
  expect_true(is.finite(objs$clamp$fn(p)))
  expect_true(all(is.finite(objs$clamp$gr(p))))

  ## the tracks held at the edge are reported, and only with clamping
  ids <- names(split(sim$dat$tags, sim$dat$tags$id))
  held <- admove:::.boundary_tags(objs$clamp, p, ids)
  expect_gt(length(held), 0)
  expect_true(all(held %in% ids))
  expect_length(admove:::.boundary_tags(objs$none, p, ids), 0)

  expect_match(admove:::.boundary_warning(held), "held at its edge")
})


test_that("admove stores the tags held at the edge of the field", {

  expect_identical(small_fit()$boundary_tags, character(0))

  ## from the report when there is one
  ids <- c("a", "b", "c")
  expect_identical(admove:::.boundary_ids(c(0, 0.5, NaN), ids), "b")
  expect_identical(admove:::.boundary_ids(NULL, ids), character(0))
})
