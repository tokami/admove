## require(admove); require(testthat)

## Smoke tests for plot_cov(): layouts, colour bar, titles, zlim


.with_null_device <- function(expr) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  force(expr)
}


test_that("plot_cov plots all layers of a covariate list in one call", {

  grid <- create_grid(xrange = c(0, 10), yrange = c(0, 10), cellsize = 1,
                      verbose = FALSE)
  geo <- suppressMessages(make_x_y_cov(grid))

  expect_s3_class(geo, "admove_cov_list")
  expect_equal(names(geo), c("x", "y"))

  expect_silent(.with_null_device(plot(geo, i = 1:2)))
  expect_silent(.with_null_device(
    plot(geo, i = 1:2, legend = FALSE, plot_contour = FALSE,
         titles = c("East", "North"), col = hcl.colors(20, "Blues"))
  ))
  expect_error(.with_null_device(plot(geo, i = 1:2, titles = "only one")),
               "one entry per covariate")
})


test_that("plot_cov handles multiple time steps, zlim and constant fields", {

  expect_silent(.with_null_device(plot_cov(skjepo$cov, select = 1:2)))
  expect_silent(.with_null_device(plot_cov(skjepo$cov, select = 1:2,
                                           zlim = c(22, 26))))
  expect_error(.with_null_device(plot_cov(skjepo$cov, select = 1:2,
                                          titles = "a")),
               "one entry per plotted time step")

  grid <- create_grid(xrange = c(0, 5), yrange = c(0, 5), cellsize = 1,
                      verbose = FALSE)
  flat <- suppressMessages(prep_cov(matrix(3, 5, 5),
                                    x_centers = x_centers(grid),
                                    y_centers = y_centers(grid), times = 0))
  expect_silent(.with_null_device(plot_cov(flat, plot_contour = FALSE)))
})
