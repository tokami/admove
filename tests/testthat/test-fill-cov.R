## fill_cov(): ring-wise filling of NA covariate cells, and setup_data(fill_na)

## 9 x 9 x 2 field, NA in the block [3:9, 3:9] (a "continent" in one corner)
toy_field <- function(dx = 1, dy = 1) {
  a <- array(outer(1:9, 1:9, "+"), c(9, 9, 2),
             dimnames = list((1:9) * dx, (1:9) * dy, c(0, 1)))
  a[, , 2] <- a[, , 2] + 10
  a[3:9, 3:9, ] <- NA
  a
}

## NA cells whose 8-neighbourhood touches a non-NA cell
touching <- function(a) {
  na <- is.na(a)
  nb <- array(FALSE, dim(a))
  for (di in -1:1) for (dj in -1:1) {
    if (di == 0 && dj == 0) next
    si <- seq_len(dim(a)[1]) + di
    sj <- seq_len(dim(a)[2]) + dj
    vi <- si >= 1 & si <= dim(a)[1]
    vj <- sj >= 1 & sj <= dim(a)[2]
    s <- array(FALSE, dim(a))
    s[which(vi), which(vj), ] <- !na[si[vi], sj[vj], , drop = FALSE]
    nb <- nb | s
  }
  na & nb
}


test_that("one ring fills exactly the cells next to data", {

  a <- toy_field()
  f1 <- fill_cov(a, verbose = FALSE)

  ring1 <- touching(a)
  expect_equal(sum(ring1), 2 * 13)
  expect_false(anyNA(f1[ring1]))
  expect_true(all(is.na(f1[is.na(a) & !ring1])))
  expect_equal(f1[!is.na(a)], a[!is.na(a)])
  expect_equal(attr(f1, "filled"), ring1)

  ## extrapolated from the neighbours: within the range of the data per slice
  for (t in 1:2) {
    r <- range(a[, , t], na.rm = TRUE)
    v <- f1[, , t][ring1[, , t]]
    expect_true(all(v >= r[1] & v <= r[2]))
  }

  ## the second ring is the next layer, grown from the first
  f2 <- fill_cov(a, n_rings = 2, verbose = FALSE)
  ring2 <- touching(unclass(f1))
  expect_equal(attr(f2, "filled"), ring1 | ring2)
  expect_equal(f2[ring1], f1[ring1])

  ## repeated calls accumulate the mask and equal filling two rings at once
  f11 <- fill_cov(f1, verbose = FALSE)
  expect_equal(attr(f11, "filled"), attr(f2, "filled"))
  expect_equal(as.numeric(f11), as.numeric(f2))
})


test_that("n_rings = Inf fills everything, 0 is the identity", {

  a <- toy_field()
  expect_false(anyNA(fill_cov(a, n_rings = Inf, verbose = FALSE)))
  expect_identical(fill_cov(a, n_rings = 0), a)

  ## a slice without data stays NA
  b <- a
  b[, , 2] <- NA
  expect_message(fb <- fill_cov(b, n_rings = Inf), "without any data")
  expect_true(all(is.na(fb[, , 2])))
  expect_false(anyNA(fb[, , 1]))
})


test_that("weights follow the cell sizes, not the cell indices", {

  ## one NA cell between a low column (x = 1) and a high one (x = 3), with
  ## y rows of value 5 around; with y cells 10 times wider than x cells, the
  ## x neighbours dominate
  a <- array(5, c(3, 3, 1), dimnames = list(1:3, (1:3) * 10, 0))
  a[1, , 1] <- 0
  a[3, , 1] <- 20
  a[2, 2, 1] <- NA
  sq <- fill_cov(unname(a), verbose = FALSE)[2, 2, 1]
  wide <- fill_cov(a, verbose = FALSE)[2, 2, 1]
  expect_equal(wide, 10, tolerance = 1e-6)
  expect_false(isTRUE(all.equal(sq, 10)))
})


test_that("class, references and list structure are kept", {

  cov <- small_sim()$dat$cov
  expect_s3_class(cov, "admove_cov_list")
  f <- fill_cov(cov, verbose = FALSE)

  expect_s3_class(f, "admove_cov_list")
  expect_s3_class(f[[1]], "admove_cov")
  expect_equal(sref(f), sref(cov))
  expect_equal(tref(f), tref(cov))
  expect_equal(dimnames(f[[1]]), dimnames(cov[[1]]))
  expect_lt(sum(is.na(f[[1]])), sum(is.na(cov[[1]])))

  expect_error(fill_cov(cov, n_rings = 1.5), "n_rings")
  expect_error(fill_cov(cov, sd = 0), "sd")
})


test_that("setup_data(fill_na = 1) keeps coastal grid cells", {

  ## 4 x 4 grid on a 4 x 4 covariate with a 2 x 2 NA corner: without filling
  ## the 4 grid cells centred on NA cells are dropped; one ring fills 3 of
  ## them, the corner cell touches no data
  grid <- create_grid(cellsize = 0.25, verbose = FALSE)
  cov <- suppressMessages(sim_cov(grid, nt = 2))
  cov[3:4, 3:4, ] <- NA
  build <- function(...) suppressWarnings(suppressMessages(
    setup_data(grid = grid, cov = cov, trange = c(0, 1), verbose = FALSE, ...)
  ))

  d0 <- build()
  d1 <- build(fill_na = 1)
  expect_equal(sum(!is.na(d0$grid$celltable)), 12)
  expect_equal(sum(!is.na(d1$grid$celltable)), 15)
  expect_equal(sum(!is.na(build(fill_na = Inf)$grid$celltable)), 16)

  expect_identical(d0$cov, build(fill_na = 0)$cov)
  expect_error(build(fill_na = -1), "fill_na")
})


test_that("fill_na = 'grid' fills until every grid cell is kept", {

  ## covariate cells 4 times finer than the grid cells, NA over the upper
  ## right grid cell: its centre lies between the second and third covariate
  ## cell inside the gap, so it needs two rings
  grid <- create_grid(cellsize = 0.5, verbose = FALSE)
  fine <- create_grid(cellsize = 0.125, verbose = FALSE)
  cov <- suppressMessages(sim_cov(fine, nt = 2))
  cov[5:8, 5:8, ] <- NA
  build <- function(...) suppressWarnings(suppressMessages(
    setup_data(grid = grid, cov = cov, trange = c(0, 1), verbose = FALSE, ...)
  ))

  expect_equal(sum(!is.na(build(fill_na = 1)$grid$celltable)), 3)
  dg <- build(fill_na = "grid")
  expect_equal(sum(!is.na(dg$grid$celltable)), 4)
  expect_equal(sum(is.na(dg$cov[[1]])), sum(is.na(build(fill_na = 2)$cov[[1]])))
  expect_error(build(fill_na = "all"), "fill_na")

  ## tags in a removed cell: the message names the covariate and the fix
  tags <- prep_tags(data.frame(id = c(1, 1, 2, 2), t = c(0.2, 0.6, 0.2, 0.6),
                               x = c(0.9, 0.3, 0.2, 0.3),
                               y = c(0.9, 0.3, 0.2, 0.3)),
                    tag_type = "d",
                    names = c(t = "t", x = "x", y = "y", id = "id"),
                    verbose = FALSE)
  expect_message(
    suppressWarnings(setup_data(grid = grid, cov = cov, tags = tags,
                                trange = c(0, 1), fill_na = 1)),
    "fill_na = \"grid\""
  )
})
