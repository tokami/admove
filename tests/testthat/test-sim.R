
## require(admove)
## require(testthat)

## set-up ---------------------------------------

grid <- create_grid()



## tests ----------------------------------------

test_that("sim_release_events", {

  out <- sim_release_events(grid)

  expect_true(ncol(out) == 3)

  ## no grid
  expect_error(sim_release_events())

  ## incorrect class
  expect_error(sim_release_events(out))
})



test_that("sim_tags", {

  out <- sim_tags("dtags")
  ## class check
  expect_s3_class(out, "admove_sim")
  expect_s3_class(out$tags, "admove_tags")

  out <- sim_tags("ctags")
  ## class check
  expect_s3_class(out, "admove_sim")
  expect_s3_class(out$tags, "admove_tags")

  ## out <- sim_tags("stags")
  ## ## class check
  ## expect_s3_class(out, "admove_sim")
  ## expect_s3_class(out$tags, "admove_tags")

  ## tags not specified
  expect_error(sim_tags())
})



test_that("sim_cov", {

  out <- sim_cov()

  ## class check
  expect_s3_class(out, "admove_cov")

  out <- admove:::.make_cov_list(out)

  ## class check
  expect_s3_class(out, "admove_cov_list")
})



## test_that("sim_data", {

##   out <- sim_data()

##   ## class check
##   expect_s3_class(out, "admove_sim")

##   ## plot check
##   skip_on_cran()
##   skip_if_not_installed("grDevices")

##   tmp <- tempfile(fileext = ".png")
##   grDevices::png(tmp, width = 800, height = 600)
##   on.exit(grDevices::dev.off(), add = TRUE)

##   expect_error(plot(out), NA)

## })
