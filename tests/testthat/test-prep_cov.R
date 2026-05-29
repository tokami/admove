

## data ---------------------------------------------

set.seed(1)

grid <- create_grid(cellsize = c(0.1,0.2))

cov_array <- sim_cov(grid)

cov_list <- lapply(seq_len(dim(cov_array)[3]),
                   function(k) cov_array[,,k, drop = FALSE][,,1])

nx <- dim(cov_array)[1]
ny <- dim(cov_array)[2]

## raster expects: [y,x]
m <- t(cov_array[,,1])
## raster expects: row 1 = max y (top)
m <- m[ny:1, , drop = FALSE]

xmin <- grid$xrange[1]
xmax <- grid$xrange[2]
ymin <- grid$yrange[1]
ymax <- grid$yrange[2]

cov_rast <- raster::raster(
  m,
  xmn = xmin, xmx = xmax,
  ymn = ymin, ymx = ymax
)


nx <- dim(cov_array)[1]
ny <- dim(cov_array)[2]
nl <- dim(cov_array)[3]

## raster expects: [y, x, layer]
a_yx <- aperm(cov_array, c(2, 1, 3))
## raster expects: row 1 is top (max y)
a_yx <- a_yx[ny:1, , , drop = FALSE]

cov_brick <- raster::brick(
  a_yx,
  xmn = xmin, xmx = xmax,
  ymn = ymin, ymx = ymax
)
names(cov_brick) <- seq_len(nl)

cov_stack <- raster::stack(cov_brick)
names(cov_stack) <- format(as.Date("2002-01-01") + seq_len(nl), "%Y-%m-%d")



## tests ---------------------------------------------

test_that("prep_cov: cov as 3D-array", {

  out <- prep_cov(cov_array)

  expect_s3_class(out, "array")
  expect_s3_class(out, "admove_cov")
  expect_equal(dim(out), dim(cov_array))

})


test_that("prep_cov: cov without dimnames", {
  cov_array2 <- cov_array
  attributes(cov_array2) <- NULL
  expect_error(prep_cov(cov_array2))
})



test_that("prep_cov: cov as list", {

  out <- prep_cov(cov_list)

  expect_s3_class(out, "array")
  expect_s3_class(out, "admove_cov")
  expect_equal(dim(out), c(dim(cov_list[[1]]),length(cov_list)))

})


test_that("prep_cov: cov as raster", {

  out <- prep_cov(cov_rast)

  expect_s3_class(out, "array")
  expect_s3_class(out, "admove_cov")
  ## raster: first y, then x
  expect_equal(dim(out)[c(1,2)], dim(cov_rast)[c(2,1)])

})


## test_that("prep_cov: cov as rasterBrick", {

##   ## warning because of invalid dates in dimnames[[3]]
##   expect_warning({
##     out <- prep_cov(cov_brick)
##   })

##   expect_s3_class(out, "array")
##   expect_s3_class(out, "admove_cov")
##   ## raster: first y, then x
##   expect_equal(dim(out)[c(1,2)], dim(cov_brick)[c(2,1)])

## })


## test_that("prep_cov: cov as rasterStack", {

##   out <- prep_cov(cov_stack,
##                   date_format = "%Y.%m.%d")

##   expect_s3_class(out, "array")
##   expect_s3_class(out, "admove_cov")
##   expect_equal(dim(out)[c(1,2)], dim(cov_stack)[c(2,1)])

## })
