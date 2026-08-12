## require(admove); require(testthat)


test_that(".lonlat_box follows the box edges in small steps", {

  skip_if_not_installed("sf")

  box <- admove:::.lonlat_box(-20, 20, 10, 30, crs = sf::st_crs(4326), step = 1)
  xy <- sf::st_coordinates(box)[, c("X", "Y")]

  ## corners are present and the ring is closed
  expect_equal(range(xy[, 1]), c(-20, 20))
  expect_equal(range(xy[, 2]), c(10, 30))
  expect_equal(xy[1, ], xy[nrow(xy), ])

  ## a corner-to-corner rectangle would have 5 points; this one steps along
  expect_gt(nrow(xy), 100)
  expect_true(all(sf::st_is_valid(box)))
})


test_that(".crop_land_to_window keeps only the land near a projected window", {

  skip_if_not_installed("sf")

  land <- sf::st_geometry(admove:::.get_land(FALSE))
  crs <- sf::st_crs(paste("+proj=aea +lat_0=50 +lon_0=-154 +lat_1=55 +lat_2=65",
                          "+x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs"))
  usr <- c(-960, -180, 550, 1260)   ## eastern Bering Sea, in km

  out <- admove:::.crop_land_to_window(land, usr, crs, 1)
  bb <- sf::st_bbox(out)

  ## the whole world went in; only the north Pacific comes out
  expect_lt(bb[["xmax"]] - bb[["xmin"]], 60)
  expect_gt(bb[["ymin"]], 40)
})


test_that("projecting the cropped land leaves the sea as sea", {

  skip_if_not_installed("sf")

  land <- sf::st_geometry(admove:::.get_land(FALSE))
  crs <- sf::st_crs(paste("+proj=aea +lat_0=50 +lon_0=-154 +lat_1=55 +lat_2=65",
                          "+x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs"))
  usr <- c(-960, -180, 550, 1260)

  ## without the crop, the rings that cross the projection's antimeridian are
  ## projected into a shape that swallows the window: everything reads as land
  whole <- sf::st_make_valid(sf::st_transform(land, crs))
  cropped <- sf::st_make_valid(
    sf::st_transform(admove:::.crop_land_to_window(land, usr, crs, 1), crs))

  sea <- sf::st_sfc(sf::st_point(c(-600, 800)), crs = crs)   ## middle of the shelf
  shore <- sf::st_sfc(sf::st_point(c(-250, 1150)), crs = crs) ## mainland Alaska

  expect_true(lengths(sf::st_intersects(sea, whole)) > 0)     ## the bug
  expect_true(lengths(sf::st_intersects(sea, cropped)) == 0)  ## the fix
  expect_true(lengths(sf::st_intersects(shore, cropped)) > 0) ## land still there
})


test_that("polygons crossing the antimeridian do not become bands of land", {

  skip_if_not_installed("sf")

  land <- sf::st_geometry(admove:::.get_land(FALSE))
  ## a window straddling +/-180: Fiji is stored as a polygon crossing it, which
  ## in the [-180, 180] frame is a sliver reaching right across the map
  crs <- sf::st_crs("+proj=aeqd +lat_0=-10 +lon_0=175 +units=km")

  out <- admove:::.crop_land_to_window(land, c(-1500, 1500, -1200, 1200), crs, 1)
  polys <- suppressWarnings(sf::st_cast(out, "POLYGON"))

  lon_span <- vapply(seq_along(polys), function(i) {
    b <- sf::st_bbox(polys[i])
    b[["xmax"]] - b[["xmin"]]
  }, numeric(1))

  ## every island in that window is a few degrees across at most
  expect_true(all(lon_span < 15))
})


test_that(".crop_land_to_window returns nothing for a window of open ocean", {

  skip_if_not_installed("sf")

  land <- sf::st_geometry(admove:::.get_land(FALSE))
  crs <- sf::st_crs("+proj=aeqd +lat_0=-30 +lon_0=-140 +units=km")

  out <- admove:::.crop_land_to_window(land, c(-400, 400, -400, 400), crs, 1)

  expect_length(out, 0L)
})
