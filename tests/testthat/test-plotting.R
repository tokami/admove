
local_plot_device <- function() {
  tf <- tempfile(fileext = ".pdf")
  grDevices::pdf(tf)
  withr::defer(grDevices::dev.off())
  invisible(tf)
}

test_that("plot_land() returns invisibly with missing CRS", {
  skip_if_not_installed("sf")

  local_plot_device()
  plot(0, 0, type = "n")

  expect_message(
    res <- plot_land(
      sref = list(crs = NULL),
      verbose = TRUE,
      warn_once = FALSE
    ),
    "No valid crs provided in sref"
  )
  expect_null(res)

})


test_that("plot_land() errors for invalid CRS", {

  skip_if_not_installed("sf")

  local_plot_device()
  plot(0, 0, type = "n")

  expect_error(
    plot_land(
      sref = list(crs = "definitely_not_a_valid_crs"),
      verbose = FALSE
    ),
    "invalid crs"
  )

})


test_that("plot_taxis() works for a minimal admove object", {

  local_plot_device()

  obj <- structure(
    list(
      dat = structure(
        list(
          grid = list(
            xrange = c(0, 1),
            yrange = c(0, 1)
          ),
          pred = list(
            time = c(1, 2),
            grid = list(
              xygrid = matrix(
                c(0.25, 0.25,
                  0.75, 0.25,
                  0.25, 0.75,
                  0.75, 0.75),
                ncol = 2,
                byrow = TRUE
              )
            )
          )
        ),
        class = "admove_data"
      ),
      pred = list(
        hTdx = matrix(c(0.1, 0.2, 0.0, 0.1,
                        0.1, 0.2, 0.0, 0.1), ncol = 2),
        hTdy = matrix(c(0.0, 0.1, 0.2, 0.1,
                        0.0, 0.1, 0.2, 0.1), ncol = 2)
      )
    ),
    class = "admove"
  )

  expect_invisible(
    plot_taxis(
      obj,
      average = TRUE,
      plot_land = FALSE
    )
  )

})


test_that("plot_taxis() errors for unsupported input", {
  expect_error(
    plot_taxis(list()),
    "Don't know how to plot taxis"
  )
})


test_that("plot_diffusion() works for a minimal admove object", {

  local_plot_device()

  obj <- structure(
    list(
      dat = structure(
        list(
          grid = list(
            xrange = c(0, 1),
            yrange = c(0, 1)
          ),
          pred = list(
            grid = list(
              xygrid = matrix(
                c(0.25, 0.25,
                  0.75, 0.25,
                  0.25, 0.75,
                  0.75, 0.75),
                ncol = 2,
                byrow = TRUE
              )
            )
          )
        ),
        class = "admove_data"
      ),
      pred = list(
        hD = matrix(log(c(0.1, 0.2, 0.15, 0.25,
                          0.1, 0.2, 0.15, 0.25)), ncol = 2)
      )
    ),
    class = "admove"
  )

  expect_invisible(
    plot_diffusion(
      obj,
      plot_land = FALSE
    )
  )

})


test_that("plot_diffusion() errors for admove_sim without par_sim", {

  local_plot_device()

  sim <- structure(
    list(
      grid = list(
        xrange = c(0, 1),
        yrange = c(0, 1)
      ),
      cov = array(1, dim = c(2, 2, 1)),
      dat = list(
        knots_tax = NULL,
        knots_dif = NULL,
        pred = list(
          grid = list(
            xygrid = matrix(c(0.5, 0.5), ncol = 2),
            igrid = 1
          )
        )
      ),
      par_sim = NULL
    ),
    class = "admove_sim"
  )

  expect_error(
    plot_diffusion(sim),
    "No parameters provided!"
  )

})
