## prep_cov() keeps the cell centres at full precision: they are the positions
## the likelihood interpolates on. Two-decimal labels put 1/12 degree cells at
## steps of 0.08 / 0.09 and shifted the whole field.

twelfth <- function() {
  xc <- -60 + (seq_len(24) - 0.5) / 12
  yc <- 10 + (seq_len(12) - 0.5) / 12
  list(xc = xc, yc = yc)
}


test_that("data.frame input keeps full-precision cell centres", {

  g <- twelfth()
  df <- expand.grid(x = g$xc, y = g$yc)
  df$sst <- 20 + df$x / 100 + df$y / 10

  cov <- suppressMessages(prep_cov(df, verbose = FALSE))
  cov1 <- if (inherits(cov, "admove_cov_list")) cov[[1]] else cov

  expect_equal(as.numeric(dimnames(cov1)[[1]]), g$xc)
  expect_equal(as.numeric(dimnames(cov1)[[2]]), g$yc)
  expect_equal(diff(range(diff(as.numeric(dimnames(cov1)[[1]])))), 0,
               tolerance = 1e-12)
})


test_that("a data.frame with one covariate column gives a named list", {

  ## the single column used to become the time label of one slice
  df <- expand.grid(x = 1:3 + 0.5, y = 1:2 + 0.5)
  df$sst <- 1:6
  cov <- prep_cov(df, verbose = FALSE)

  expect_s3_class(cov, "admove_cov_list")
  expect_named(cov, "sst")
  expect_equal(dimnames(cov[[1]])[[3]], "0")
  expect_equal(unname(cov[[1]][, , 1]), matrix(1:6, 3, 2))
  expect_equal(dimnames(prep_cov(df, times = 5, verbose = FALSE)[[1]])[[3]],
               "5")

  ## as with several columns
  df$chl <- 6:1
  expect_named(prep_cov(df, verbose = FALSE), c("sst", "chl"))
})


test_that("SpatRaster input keeps full-precision cell centres", {

  skip_if_not_installed("terra")
  g <- twelfth()
  r <- terra::rast(nrows = length(g$yc), ncols = length(g$xc),
                   xmin = -60, xmax = -58, ymin = 10, ymax = 11,
                   crs = "EPSG:4326")
  terra::values(r) <- seq_len(terra::ncell(r))

  cov <- suppressMessages(prep_cov(r, verbose = FALSE))
  cov1 <- if (inherits(cov, "admove_cov_list")) cov[[1]] else cov

  expect_equal(as.numeric(dimnames(cov1)[[1]]), g$xc)
  expect_equal(as.numeric(dimnames(cov1)[[2]]), g$yc)
})


test_that("sim_cov time labels are the exact slice starts", {

  grid <- create_grid(cellsize = 0.25, verbose = FALSE)
  cov <- suppressMessages(sim_cov(grid, nt = 4, trange = c(0, 1)))

  expect_equal(as.numeric(dimnames(cov)[[3]]), seq(0, 1, length.out = 4))
})
