## require(admove); require(testthat)

## Guards against a time axis that is silently wrong: a tref attached to an
## object whose own tref fields were NA (values re-interpreted, not converted)
## and tag times that t2index() clamps onto the last covariate slice.


test_that("add_tref warns when it can only relabel an NA tref", {

  dtags <- prep_tags(skjepo$dtags, tag_type = "d",
                     names = c(t = "time", x = "mptlon", y = "mptlat"),
                     verbose = FALSE)

  ## numeric time with no date information -> origin and units stay NA
  expect_true(admove:::.is_na_scalar(tref(dtags)$origin))

  expect_warning(
    add_tref(dtags, create_tref("2000-01-01 UTC", units = "months")),
    "re-interpreted, not converted"
  )
})


test_that("add_tref stays silent when relabelling with verbose = FALSE", {

  dtags <- prep_tags(skjepo$dtags, tag_type = "d",
                     names = c(t = "time", x = "mptlon", y = "mptlat"),
                     verbose = FALSE)

  expect_silent(
    add_tref(dtags, create_tref("2000-01-01 UTC", units = "months"),
             verbose = FALSE)
  )
})


test_that("add_tref does not warn when the object has a real tref", {

  ctags <- prep_tags(skjepo$ctags, tag_type = "c",
                     names = c(t0 = "date_time", t1 = "date_caught",
                               x0 = "rel_lon", x1 = "recap_lon",
                               y0 = "rel_lat", y1 = "recap_lat"),
                     date_origin = "1899-12-30", verbose = FALSE)

  expect_no_warning(
    add_tref(ctags, create_tref("2000-01-01 UTC", units = "months"),
             verbose = FALSE, shift_origin = TRUE)
  )
})


test_that("setup_data warns when tag times are clamped to the last cov slice", {

  ctags <- prep_tags(skjepo$ctags, tag_type = "c",
                     names = c(t0 = "date_time", t1 = "date_caught",
                               x0 = "rel_lon", x1 = "recap_lon",
                               y0 = "rel_lat", y1 = "recap_lat"),
                     date_origin = "1899-12-30", verbose = FALSE)
  grid <- create_grid(x = skjepo$grid, cellsize = c(10, 10), plot = FALSE,
                      verbose = FALSE)
  cov <- prep_cov(skjepo$cov, verbose = FALSE)

  ## in range (a tag in the last slice may sit up to one slice spacing beyond
  ## the last label, which is normal and must not warn)
  expect_no_warning(
    setup_data(grid = grid, cov = cov, tags = ctags, transform_sref = TRUE,
               shift_tref = TRUE, verbose = FALSE)
  )

  ## well past the end of the covariate
  ctags$t <- ctags$t + 50

  expect_warning(
    setup_data(grid = grid, cov = cov, tags = ctags, transform_sref = TRUE,
               shift_tref = TRUE, verbose = FALSE),
    "lie beyond the last time slice"
  )
})
