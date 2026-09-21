## require(admove); require(testthat)

## create_grid(): xrange / yrange with and without 'x', and which units they are in

## metre-based projection, coordinates stored in km (as in the IOTC analysis)
crs_aeqd_m <- paste("+proj=aeqd +lat_0=-5 +lon_0=55 +x_0=0 +y_0=0",
                    "+datum=WGS84 +units=m +no_defs")

.km_tags <- function() {
  d <- data.frame(t = 1:4,
                  x = c(-3000, 7000, 0, 4000),
                  y = c(-4000, 4000, 1000, -1000),
                  id = 1)
  prep_dtags(d, names = c(t = "t", x = "x", y = "y", id = "id"),
             sref = list(crs = crs_aeqd_m, units = "km"), verbose = FALSE)
}


test_that("xrange / yrange take precedence over the extent of tags", {
  skip_if_not_installed("sf")

  tags <- suppressMessages(.km_tags())

  expect_message(
    grid <- create_grid(tags, xrange = c(-400, 4700), yrange = c(-2000, 3700),
                        cellsize = 500),
    "2 of 4 tag positions lie outside")

  expect_equal(grid$xrange, c(-400, 5100))
  expect_equal(grid$yrange, c(-2000, 4000))
  expect_equal(units_space(grid), "km")
  expect_equal(crs_scale(grid), 0.001)

  ## one axis only: the other still comes from the tags
  grid <- create_grid(tags, xrange = c(-400, 4700), cellsize = 500,
                      verbose = FALSE)
  expect_equal(grid$xrange, c(-400, 5100))
  expect_lt(grid$yrange[1], min(tags$y))
  expect_gt(grid$yrange[2], max(tags$y))
})


test_that("without ranges the grid still encloses every tag with a margin", {
  skip_if_not_installed("sf")

  tags <- suppressMessages(.km_tags())
  grid <- create_grid(tags, cellsize = 500, verbose = FALSE)

  expect_lt(grid$xrange[1], min(tags$x))
  expect_gt(grid$xrange[2], max(tags$x))
  expect_lt(grid$yrange[1], min(tags$y))
  expect_gt(grid$yrange[2], max(tags$y))
})


test_that("without x, ranges are in 'units', whatever form the crs takes", {
  skip_if_not_installed("sf")

  make <- function(crs) {
    create_grid(xrange = c(-400, 4700), yrange = c(-2000, 3700),
                cellsize = 500, crs = crs, units = "km", verbose = FALSE)
  }
  g_str <- make(crs_aeqd_m)
  g_sf <- make(sf::st_crs(crs_aeqd_m))
  g_wkt <- make(sf::st_crs(crs_aeqd_m)$wkt)

  ## an sf::crs used to build the grid in metres and shrink it by 1000
  expect_equal(g_sf$xrange, c(-400, 5100))
  expect_equal(g_sf$cellsize, c(500, 500))
  expect_equal(crs_scale(g_sf), 0.001)
  expect_equal(g_str$xygrid, g_sf$xygrid)
  expect_equal(g_wkt$xygrid, g_sf$xygrid)
  expect_equal(sref(g_str), sref(g_sf))
})


test_that("without x, crs_scale also declares the units of the ranges", {
  skip_if_not_installed("sf")

  grid <- create_grid(xrange = c(2, 5), yrange = c(1.5, 5), cellsize = 0.5,
                      crs = sf::st_crs(3035), crs_scale = 1e-6, verbose = FALSE)

  expect_equal(grid$xrange, c(2, 5))
  expect_equal(crs_scale(grid), 1e-6)
  expect_equal(units_space(grid), "metre_x_1e-06")
})


test_that("without units the ranges are in the unit of the crs", {
  skip_if_not_installed("sf")

  grid <- create_grid(xrange = c(2e6, 5e6), yrange = c(1.5e6, 5e6),
                      cellsize = 5e5, crs = sf::st_crs(3035), verbose = FALSE)

  expect_equal(grid$xrange, c(2e6, 5e6))
  expect_equal(crs_scale(grid), 1)
  expect_equal(units_space(grid), "metre")
})
