
set.seed(1)

df_wide <- data.frame(time0 = 1:10 + 0.03,
                      time1 = 2:11 + 0.3,
                      lon_rel = rnorm(10),
                      lat_rel = rnorm(10),
                      lon_rec = rnorm(10),
                      lat_rec = rnorm(10))

df_long <- data.frame(time = 1:20 + 0.543,
                      lon = rnorm(20),
                      lat = rnorm(20),
                      id = rep(1:10, each = 2),
                      other = sample(letters, 20, TRUE))

df_list <- split(df_long, df_long$id)
df_list <- lapply(df_list, function(x) x[,-which(colnames(x) == "id")])


test_that("prep_ctags: wide format", {

  out <- prep_ctags(df_wide, names = c(t0 = "time0",
                                       y0 = "lat_rel",
                                       x0 = "lon_rel",
                                       x1 = "lon_rec",
                                       t1 = "time1",
                                       y1 = "lat_rec"))

  expect_s3_class(out, "data.frame")
  expect_s3_class(out, "admove_tags")
  expect_equal(nrow(out), nrow(df_wide)*2)
  expect_true(all(c("t","x","y") %in% names(out)))
  expect_type(out$t, "double")
  expect_type(out$x, "double")
  expect_type(out$y, "double")
  expect_true(is.numeric(out$x))
  expect_true(is.numeric(out$y))
  expect_true(is.numeric(out$t))

})


test_that("prep_ctags: required id not provided", {

  expect_error(
    prep_ctags(df_long, names = c(y = "lat",
                                  x = "lon",
                                  t = "time"))
    )

})


test_that("prep_ctags: missing id with list is okay", {

  ## no error if id is missing if split into list
  out <- prep_ctags(df_list,
                    names = c(t = "time",
                              y = "lat",
                              x = "lon"))
  ## add use
  expect_true("use" %in% names(out))
  expect_true(any(out$use %in% c(FALSE,1,TRUE,0)))
  ## keep existing columns
  expect_true("other" %in% names(out))

})
