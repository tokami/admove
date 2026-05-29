
set.seed(1)

df_long <- data.frame(time = 1:50 + 0.543,
                      id = rep(1:10, each = 5),
                      lat = rnorm(50),
                      lon = rnorm(50),
                      other = sample(letters, 50, TRUE))

df_list <- split(df_long, df_long$id)
df_list <- lapply(df_list, function(x) x[,-which(colnames(x) == "id")])



test_that("prep_stags: long format with id", {

  out <- prep_stags(df_long,
                    names = c(t = "time",
                              y = "lat",
                              id = "id",
                              x = "lon"))

  expect_s3_class(out, "data.frame")
  expect_s3_class(out, "admove_tags")
  expect_equal(nrow(out), nrow(df_long))
  expect_true(all(c("t","x","y") %in% names(out)))
  expect_type(out$t, "double")
  expect_type(out$x, "double")
  expect_type(out$y, "double")
  expect_true(is.numeric(out$x))
  expect_true(is.numeric(out$y))
  expect_true(is.numeric(out$t))

})



test_that("prep_stags: required id not provided", {

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
