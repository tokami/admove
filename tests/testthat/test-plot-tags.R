## require(admove); require(testthat)


test_that("plot_tags() draws subsets of elements via 'show'", {

  tags <- skjepo$sim$tags

  pdf(NULL)
  on.exit(dev.off())

  expect_no_error(plot_tags(tags))
  expect_no_error(plot_tags(tags, show = "release"))
  expect_no_error(plot_tags(tags, show = "recovery"))
  expect_no_error(plot_tags(tags, show = c("release", "recovery")))
  expect_no_error(plot_tags(tags, show = "path", by_tag_type = FALSE))
  expect_no_error(plot_tags(tags, legend = FALSE))
  expect_no_error(plot_tags(tags, show = "release", col = "red"))
  expect_no_error(plot_tags(tags, show = c("path", "recovery"), lty = 2))
})


test_that("plot_tags() rejects unknown 'show' elements", {

  pdf(NULL)
  on.exit(dev.off())

  expect_error(plot_tags(skjepo$sim$tags, show = "releases_only"))
})
