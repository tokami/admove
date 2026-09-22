
## crs, units_space and crs_scale describe one relation together:
##   stored_xy = T(lonlat, crs) * crs_scale
## so crs() must keep returning the CRS of the UNSCALED coordinates.
## See dev/code_notes.org, "Spatial reference: CRS units vs stored units".

crs_aeqd_m <- paste("+proj=aeqd +lat_0=40 +lon_0=-6 +x_0=0 +y_0=0",
                    "+datum=WGS84 +units=m +no_defs")
crs_aeqd_km <- sub("units=m ", "units=km ", crs_aeqd_m)


test_that("a metre CRS with km storage keeps the metre CRS", {
  skip_if_not_installed("sf")

  sp <- create_sref(crs = crs_aeqd_m, units = "km")

  expect_equal(units_space(sp), "km")
  expect_equal(crs_scale(sp), 0.001)
  ## crs() does not adopt the stored unit
  expect_equal(sf::st_crs(crs(sp))$units_gdal, "metre")
})


test_that("a km-native CRS with km storage needs no scaling", {
  skip_if_not_installed("sf")

  sp <- create_sref(crs = crs_aeqd_km, units = "km")

  expect_equal(units_space(sp), "km")
  expect_equal(crs_scale(sp), 1)
  expect_equal(sf::st_crs(crs(sp))$units_gdal, "kilometre")
})


test_that("create_sref() derives units and crs_scale from each other", {
  skip_if_not_installed("sf")

  ## units given -> crs_scale derived
  expect_equal(crs_scale(create_sref(32631, units = "km")), 0.001)

  ## crs_scale given -> units derived
  expect_equal(units_space(create_sref(32631, crs_scale = 0.001)), "km")

  ## neither given -> both from the CRS
  sp <- create_sref(32631)
  expect_equal(units_space(sp), "metre")
  expect_equal(crs_scale(sp), 1)

  ## lon/lat has no metric scaling
  sp <- create_sref(4326)
  expect_equal(units_space(sp), "degree")
  expect_equal(crs_scale(sp), 1)

  ## neither crs nor units
  sp <- create_sref()
  expect_true(is.na(crs(sp)))
  expect_true(is.na(crs_scale(sp)))
})


test_that("create_sref() keeps an explicitly supplied crs_scale", {
  skip_if_not_installed("sf")

  expect_warning(sp <- create_sref(32631, units = "m", crs_scale = 0.001),
                 "does not match")
  expect_equal(crs_scale(sp), 0.001)
  expect_equal(units_space(sp), "m")

  ## consistent pair: no warning
  expect_silent(create_sref(32631, units = "km", crs_scale = 0.001))
})


test_that("scalings without a unit name round-trip through the composite label", {
  skip_if_not_installed("sf")

  sp <- create_sref(32631, crs_scale = 1e-6)

  expect_equal(units_space(sp), "metre_x_1e-06")
  expect_equal(crs_scale(sp), 1e-6)
  ## .in_m() must read back what .units_from_scale() wrote
  expect_equal(.in_m(units_space(sp)), 1e6)
  expect_equal(.infer_sref_scale(crs(sp), units_space(sp)), 1e-6)
})


test_that(".units_from_scale() inverts .in_m()", {

  expect_equal(.units_from_scale("metre", 1e-3), "km")
  expect_equal(.units_from_scale("kilometre", 1e3), "m")
  expect_equal(.units_from_scale("metre", 1), "metre")
  expect_equal(.units_from_scale("metre", 0.5), "metre_x_5e-01")

  for (scale in c(1e-6, 1e-3, 0.5, 1, 2, 1e3)) {
    u <- .units_from_scale("metre", scale)
    expect_equal(.in_m("metre") / .in_m(u), scale, tolerance = 1e-12)
  }

  expect_true(is.na(.units_from_scale(NA_character_, 1)))
  expect_true(is.na(.units_from_scale("metre", NA_real_)))
})


test_that("sref setters recompute the dependent fields", {
  skip_if_not_installed("sf")

  sp <- create_sref(crs = 32631)

  units_space(sp) <- "km"
  expect_equal(crs_scale(sp), 0.001)

  ## same native unit: the stored label survives the CRS change
  crs(sp) <- 3035
  expect_equal(units_space(sp), "km")
  expect_equal(crs_scale(sp), 0.001)

  ## different native unit: re-inferred
  crs(sp) <- 4326
  expect_equal(units_space(sp), "degree")
  expect_equal(crs_scale(sp), 1)

  ## setting the scale relabels the units
  sp <- create_sref(crs = 32631)
  crs_scale(sp) <- 0.001
  expect_equal(units_space(sp), "km")

  ## crs<- normalises to WKT, like create_sref()
  sp <- create_sref()
  crs(sp) <- 32631
  expect_equal(crs(sp), sf::st_crs(32631)$wkt)
})


test_that("sref setters refuse objects that carry coordinates", {
  skip_if_not_installed("sf")

  grid <- create_grid(xrange = c(0, 5e4), yrange = c(0, 5e4),
                      cellsize = 5e3, crs = 32631, verbose = FALSE)

  expect_error(units_space(grid) <- "km", "scale_sref")
  expect_error(crs_scale(grid) <- 0.001, "scale_sref")
  expect_error(crs(grid) <- 4326, "transform_sref")
})


test_that("scale_sref() moves coordinates, crs_scale and label together", {
  skip_if_not_installed("sf")

  grid <- create_grid(xrange = c(0, 5e4), yrange = c(0, 5e4),
                      cellsize = 5e3, crs = 32631, verbose = FALSE)

  km <- scale_sref(grid, scale = 0.001, verbose = FALSE)
  expect_equal(units_space(km), "km")
  expect_equal(crs_scale(km), 0.001)
  expect_equal(km$xrange, grid$xrange * 0.001)
  expect_equal(km$cellsize, grid$cellsize * 0.001)

  back <- scale_sref(km, scale = 1000, verbose = FALSE)
  expect_equal(units_space(back), "m")
  expect_equal(crs_scale(back), 1)
  expect_equal(as.matrix(back$xygrid), as.matrix(grid$xygrid))

  ## the units= form goes through .guess_crs_scale()
  km2 <- scale_sref(grid, units = "km", verbose = FALSE)
  expect_equal(crs_scale(km2), 0.001)
  expect_equal(km2$xrange, grid$xrange * 0.001)
})


test_that("add_sref() keeps an explicit crs_scale and rescales the grid", {
  skip_if_not_installed("sf")

  grid <- create_grid(xrange = c(0, 5e4), yrange = c(0, 5e4),
                      cellsize = 5e3, crs = 32631, verbose = FALSE)

  res <- add_sref(grid, list(crs = 32631, units = "km", crs_scale = 0.001),
                  verbose = FALSE)

  expect_equal(units_space(res), "km")
  expect_equal(crs_scale(res), 0.001)
  expect_equal(res$xrange, grid$xrange * 0.001)
})


test_that("validate_sref() rejects an unusable crs_scale", {
  sp <- create_sref()
  sp$crs_scale <- -1
  expect_error(validate_sref(sp), "positive finite")
})


test_that("print.admove_sref shows the CRS unit next to the stored unit", {
  skip_if_not_installed("sf")

  expect_snapshot(print(create_sref(crs = crs_aeqd_m, units = "km")))
})
