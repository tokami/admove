## require(admove); require(testthat)

## Land/water mask in create_grid()

## tags spanning the west coast of Africa (Gulf of Guinea): part land, part sea
.mask_tags <- function(crs = 4326) {
  d <- data.frame(t = 1:4,
                  x = c(-15, 5, -10, 10),
                  y = c(-5, 15, 10, -2),
                  id = 1)
  prep_dtags(d, names = c(t = "t", x = "x", y = "y", id = "id"),
             sref = list(crs = crs), verbose = FALSE)
}

.on_land <- function(grid) {
  s2_old <- suppressMessages(sf::sf_use_s2(FALSE))
  on.exit(suppressMessages(sf::sf_use_s2(s2_old)), add = TRUE)
  land <- suppressMessages(sf::st_union(sf::st_make_valid(sf::st_geometry(admove:::.get_land()))))
  pts <- sf::st_as_sf(grid$xygrid, coords = c("x", "y"), crs = 4326)
  suppressMessages(lengths(sf::st_intersects(pts, land)) > 0)
}


test_that("mask = 'water' drops land cells and 'land' keeps their complement", {
  skip_if_not_installed("sf")

  tags <- .mask_tags()
  g_all <- create_grid(tags, cellsize = 2, verbose = FALSE)
  g_water <- create_grid(tags, cellsize = 2, mask = "water", verbose = FALSE)
  g_land <- create_grid(tags, cellsize = 2, mask = "land", verbose = FALSE)

  expect_s3_class(g_water, "admove_grid")
  expect_lt(nrow(g_water$xygrid), nrow(g_all$xygrid))
  expect_gt(nrow(g_water$xygrid), 0)
  expect_gt(nrow(g_land$xygrid), 0)
  expect_equal(nrow(g_water$xygrid) + nrow(g_land$xygrid), nrow(g_all$xygrid))

  expect_false(any(.on_land(g_water)))
  expect_true(all(.on_land(g_land)))

  ## celltable indexes the kept cells consecutively
  expect_equal(sort(g_water$celltable[!is.na(g_water$celltable)]),
               seq_len(nrow(g_water$xygrid)))
})


test_that("mask works with a projected sref (aeqd, km)", {
  skip_if_not_installed("sf")

  crs_aeqd <- "+proj=aeqd +lat_0=5 +lon_0=0 +units=km +datum=WGS84"
  d <- data.frame(t = 1:4, x = c(-15, 5, -10, 10), y = c(-5, 15, 10, -2))
  xy <- sf::st_coordinates(sf::st_transform(
    sf::st_as_sf(d, coords = c("x", "y"), crs = 4326), crs_aeqd))
  d$x <- xy[, 1]
  d$y <- xy[, 2]
  d$id <- 1
  tags <- prep_dtags(d, names = c(t = "t", x = "x", y = "y", id = "id"),
                     sref = list(crs = crs_aeqd), verbose = FALSE)

  g_all <- create_grid(tags, cellsize = 200, verbose = FALSE)
  g_water <- create_grid(tags, cellsize = 200, mask = "water", verbose = FALSE)

  expect_lt(nrow(g_water$xygrid), nrow(g_all$xygrid))
  expect_gt(nrow(g_water$xygrid), 0)

  ## kept cell centres, back in lon/lat, are all off land
  ll <- sf::st_coordinates(sf::st_transform(
    sf::st_as_sf(g_water$xygrid, coords = c("x", "y"), crs = crs_aeqd), 4326))
  expect_false(any(.on_land(list(xygrid = data.frame(x = ll[, 1], y = ll[, 2])))))
})


test_that("mask requires a CRS and a valid value", {
  skip_if_not_installed("sf")

  d <- data.frame(t = 1:2, x = c(0, 10), y = c(0, 10), id = 1)
  tags <- prep_dtags(d, names = c(t = "t", x = "x", y = "y", id = "id"),
                     verbose = FALSE)
  expect_error(create_grid(tags, cellsize = 2, mask = "water", verbose = FALSE),
               "requires a CRS")

  expect_error(create_grid(.mask_tags(), cellsize = 2, mask = "sea", verbose = FALSE))
})


test_that("mask reports tag positions in removed cells", {
  skip_if_not_installed("sf")

  ## (5, 15) is on land in Nigeria/Niger
  expect_message(create_grid(.mask_tags(), cellsize = 2, mask = "water"),
                 "tag position")
})
