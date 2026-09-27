## aggregate_cov(): block means of covariate fields

## The block-average used in the dolphinfish analysis before aggregate_cov()
## existed (admove_dolphinfish_v2.R), as the reference.
script_aggregate <- function(x, res = 1) {
  bx <- floor(as.numeric(dimnames(x)[[1]]) / res)
  by <- floor(as.numeric(dimnames(x)[[2]]) / res)
  ux <- sort(unique(bx))
  uy <- sort(unique(by))
  ix <- match(bx, ux)
  iy <- match(by, uy)
  out <- array(NA_real_, c(length(ux), length(uy), dim(x)[3]))
  for (t in seq_len(dim(x)[3])) {
    M <- x[, , t]
    s <- t(rowsum(t(rowsum(ifelse(is.na(M), 0, M), ix)), iy))
    n <- t(rowsum(t(rowsum((!is.na(M)) * 1, ix)), iy))
    out[, , t] <- ifelse(n > 0, s / n, NA_real_)
  }
  out
}

make_fine <- function() {
  xc <- seq(-99.875, -90.125, by = 0.25)
  yc <- seq(10.125, 14.875, by = 0.25)
  a <- withr::with_seed(1, array(rnorm(length(xc) * length(yc) * 3, 25),
                                 c(length(xc), length(yc), 3),
                                 dimnames = list(xc, yc, c("0", "1", "2"))))
  a[1:6, 1:4, ] <- NA
  a[10, 10, 2] <- NA
  a
}


test_that("cellsize blocks match the reference block average", {

  a <- make_fine()
  out <- aggregate_cov(a, cellsize = 1, verbose = FALSE)

  expect_equal(unname(unclass(out)), script_aggregate(a, 1))
  expect_equal(as.numeric(dimnames(out)[[1]]), seq(-99.5, -90.5, by = 1))
  expect_equal(as.numeric(dimnames(out)[[2]]), seq(10.5, 14.5, by = 1))
  expect_equal(dimnames(out)[[3]], dimnames(a)[[3]])
})


test_that("factor gives the same blocks when they align", {

  a <- make_fine()
  expect_equal(aggregate_cov(a, factor = 4, verbose = FALSE),
               aggregate_cov(a, cellsize = 1, verbose = FALSE))
})


test_that("min_frac sets mostly missing blocks to NA", {

  a <- make_fine()
  ## block (2, 1) holds cells 5:8 x 1:4, of which 5:6 x 1:4 are NA (8 of 16)
  blk <- a[5:8, 1:4, 1]
  expect_equal(aggregate_cov(a, cellsize = 1, verbose = FALSE)[2, 1, 1],
               mean(blk, na.rm = TRUE))
  expect_true(is.na(aggregate_cov(a, cellsize = 1, min_frac = 0.6,
                                  verbose = FALSE)[2, 1, 1]))
  expect_false(is.na(aggregate_cov(a, cellsize = 1, min_frac = 0.5,
                                   verbose = FALSE)[2, 1, 1]))
  ## block (1, 1) is entirely missing whatever min_frac is
  expect_true(is.na(aggregate_cov(a, cellsize = 1, verbose = FALSE)[1, 1, 1]))
})


test_that("time_factor averages consecutive slices, labelled by the first", {

  a <- make_fine()
  out <- aggregate_cov(a, time_factor = 2, verbose = FALSE)

  expect_equal(dim(out), c(dim(a)[1:2], 2L))
  expect_equal(dimnames(out)[[3]], c("0", "2"))
  expect_equal(out[20, 5, 1], mean(a[20, 5, 1:2]))
  ## a cell missing in one slice takes the other
  expect_equal(out[10, 10, 1], a[10, 10, 1])
  expect_equal(out[20, 5, 2], a[20, 5, 3])
})


test_that("admove_cov objects and lists keep their class and references", {

  cov <- skjepo$sim$cov
  out <- aggregate_cov(cov, factor = 2, time_factor = 2, verbose = FALSE)

  expect_s3_class(out, "admove_cov")
  expect_equal(sref(out), sref(cov))
  expect_equal(tref(out), tref(cov))
  expect_equal(dim(out), c(ceiling(dim(cov)[1:2] / 2), dim(cov)[3] / 2))
  ## the new cells are a regular lattice, twice as coarse
  dx <- diff(as.numeric(dimnames(cov)[[1]]))[1]
  expect_equal(diff(as.numeric(dimnames(out)[[1]])),
               rep(2 * dx, dim(out)[1] - 1))

  lst <- .make_cov_list(cov)
  names(lst) <- "sst"
  out_l <- aggregate_cov(lst, factor = 2, verbose = FALSE)
  expect_s3_class(out_l, "admove_cov_list")
  expect_equal(unclass(out_l[[1]]), unclass(aggregate_cov(cov, factor = 2,
                                                          verbose = FALSE)))

  ## the result can be used as data
  expect_s3_class(suppressWarnings(suppressMessages(
    setup_data(cov = out, verbose = FALSE))), "admove_data")
})


test_that("aggregate_cov drops the fill marks and handles matrices", {

  a <- make_fine()
  f <- fill_cov(a, verbose = FALSE)
  expect_null(attr(aggregate_cov(f, factor = 2, verbose = FALSE), "filled"))

  m <- a[, , 1]
  out <- aggregate_cov(m, cellsize = 1, verbose = FALSE)
  expect_equal(dim(out), c(10L, 5L))
  expect_equal(unname(out), script_aggregate(a, 1)[, , 1])
})


test_that("labels rounded for printing give the same blocks, silently", {

  ## as in the dolphinfish script: sprintf("%.2f") turns -99.875 into -99.88,
  ## so the label steps alternate 0.26 / 0.24
  a <- make_fine()
  r <- a
  dimnames(r)[1:2] <- lapply(dimnames(a)[1:2],
                             function(v) sprintf("%.2f", as.numeric(v)))

  expect_silent(out <- aggregate_cov(r, cellsize = 1, verbose = FALSE))
  expect_equal(unclass(out), unclass(aggregate_cov(a, cellsize = 1,
                                                   verbose = FALSE)))
  expect_message(aggregate_cov(r, cellsize = 1), "cells of 0.25 x 0.25")
})


test_that("aggregate_cov checks its arguments", {

  a <- make_fine()
  expect_error(aggregate_cov(a, verbose = FALSE), "Nothing to aggregate")
  expect_error(aggregate_cov(a, cellsize = 1, factor = 2), "not both")
  expect_error(aggregate_cov(a, cellsize = 0.1), "only makes cells coarser")
  expect_error(aggregate_cov(a, factor = 1.5), "whole numbers")
  expect_error(aggregate_cov(a, factor = 2, min_frac = 2), "between 0 and 1")
  expect_warning(aggregate_cov(a, cellsize = c(0.6, 1), verbose = FALSE),
                 "not a multiple of the x cell size")
  expect_message(aggregate_cov(a, cellsize = 1),
                 "40 x 20 cells of 0.25 x 0.25 -> 10 x 5 of 1 x 1")
})
