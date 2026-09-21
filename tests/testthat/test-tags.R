## require(admove); require(testthat)



make_tags <- function() {
  data.frame(
    id = c("c_ok", "c_ok",
           "c_reversed", "c_reversed",
           "c_zero", "c_zero",
           rep("d_unsorted", 4),
           rep("d_dup", 4)),
    t = c(0, 10,
          10, 2,
          5, 5,
          3, 1, 4, 2,
          1, 2, 2, 3),
    x = c(1, 2, 1, 2, 1, 2, 30, 10, 40, 20, 1, 2, 99, 3),
    y = c(1, 2, 1, 2, 1, 2, 30, 10, 40, 20, 1, 2, 99, 3),
    tag_type = c(rep("c", 6), rep("d", 8))
  )
}


test_that(".format_ids truncates long id lists", {

  expect_equal(.format_ids(1:3), "1, 2, 3")
  expect_equal(.format_ids(1:5), "1, 2, 3, 4, 5")
  expect_equal(.format_ids(1:8), "1, 2, 3, 4, 5, ... and 3 more")
  expect_equal(.format_ids(1:8, max_show = 2), "1, 2, ... and 6 more")
  expect_equal(.format_ids(character(0)), "")
})


test_that("check_tags removes recaptures at or before release", {

  out <- check_tags(make_tags(), verbose = FALSE)

  expect_false(any(out$id %in% c("c_reversed", "c_zero")))
  expect_true(all(c("c_ok", "d_unsorted", "d_dup") %in% out$id))
})


test_that("check_tags reorders archival tags without mixing up positions", {

  out <- check_tags(make_tags(), verbose = FALSE)
  tag <- out[out$id == "d_unsorted", ]

  expect_equal(tag$t, c(1, 2, 3, 4))
  expect_equal(tag$x, c(10, 20, 30, 40))
  expect_equal(tag$y, c(10, 20, 30, 40))
})


test_that("check_tags drops repeated times within a tag", {

  out <- check_tags(make_tags(), verbose = FALSE)
  tag <- out[out$id == "d_dup", ]

  expect_equal(tag$t, c(1, 2, 3))
  expect_false(any(out$x == 99))
})


test_that("check_tags leaves valid tags untouched", {

  tags <- make_tags()
  tags <- tags[tags$id %in% c("c_ok", "d_unsorted"), ]
  tags <- tags[order(tags$id, tags$t), ]
  rownames(tags) <- NULL

  out <- check_tags(tags, verbose = FALSE)

  expect_equal(nrow(out), nrow(tags))
  expect_equal(out$t, tags$t)
  expect_equal(out$x, tags$x)
})


test_that("check_tags keeps single-observation tags when asked to", {

  tags <- make_tags()
  tags <- rbind(tags, data.frame(id = "d_single", t = 1, x = 1, y = 1,
                                 tag_type = "d"))

  out <- check_tags(tags, remove_non_recovered_tags = FALSE, verbose = FALSE)

  expect_true("d_single" %in% out$id)
})


test_that("setup_data accepts tags combined with list() as well as c()", {

  ctags <- prep_tags(skjepo$ctags, tag_type = "c",
                     names = c(t0 = "date_time", t1 = "date_caught",
                               x0 = "rel_lon", x1 = "recap_lon",
                               y0 = "rel_lat", y1 = "recap_lat"),
                     date_origin = "1899-12-30", verbose = FALSE)
  dtags <- prep_tags(skjepo$dtags, tag_type = "d",
                     names = c(t = "time", x = "mptlon", y = "mptlat"),
                     date_origin = "1899-12-30", verbose = FALSE)

  .setup <- function(tags) {
    suppressMessages(suppressWarnings(
      setup_data(tags = tags, shift_tref = TRUE, verbose = FALSE)
    ))
  }

  dat_c <- .setup(c(ctags = ctags, dtags = dtags))
  dat_list <- .setup(list(ctags = ctags, dtags = dtags))

  expect_identical(dat_list, dat_c)
  expect_error(.setup(as.data.frame(dtags)), "must be an 'admove_tags' object")
})


test_that("get_recaptured_tags() splits recaptured and never-recaptured tags", {

  w <- data.frame(id = c("a", "b", "c"),
                  t0 = c(0, 0, 0), t1 = c(1, NA, 2),
                  x0 = c(0, 1, 2), x1 = c(1, NA, 3),
                  y0 = c(0, 1, 2), y1 = c(1, NA, 3))
  tags <- suppressMessages(
    prep_ctags(w, names = c(id = "id", t0 = "t0", t1 = "t1", x0 = "x0",
                            x1 = "x1", y0 = "y0", y1 = "y1"),
               verbose = FALSE))

  rec <- get_recaptured_tags(tags)
  expect_setequal(unique(rec$id), c("a", "c"))
  expect_equal(nrow(rec), 4L)
  expect_s3_class(rec, "admove_tags")
  expect_identical(sref(rec), sref(tags))
  expect_identical(tref(rec), tref(tags))

  non <- get_recaptured_tags(tags, invert = TRUE)
  expect_identical(unique(non$id), "b")
  expect_equal(nrow(non), 2L)
})
