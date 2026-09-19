## KF: predicted positions that leave the covariate field (conf$kf_boundary)

test_that(".clamp_ad clamps values and is the identity inside the bounds", {

  v <- c(-3, -1, 0, 0.5, 1, 2)
  expect_equal(admove:::.clamp_ad(v, -1, 1), c(-1, -1, 0, 0.5, 1, 1))

  ## on the tape: derivative 1 inside, 0 outside
  f <- RTMB::MakeTape(function(p) admove:::.clamp_ad(p, -1, 1), c(0.5, 2, -3))
  expect_equal(as.numeric(f(c(0.5, 2, -3))), c(0.5, 1, -1))
  expect_equal(diag(f$jacobian(c(0.5, 2, -3))), c(1, 0, 0))
})


test_that("check_conf validates kf_boundary", {

  dat <- small_sim()$dat
  expect_equal(default_conf(dat, verbose = FALSE)$kf_boundary, "clamp")
  expect_error(check_conf(list(kf_boundary = "reflect"), dat, verbose = FALSE),
               "kf_boundary")
})


## One tape per setting, built at taxis coefficients 10 times the simulated
## ones, which push some predicted tracks out of the covariate field.
##
## Two fixtures, and the difference between them is the point. small_sim()'s
## covariate carries the skjepo land mask (280 NA cells of 1248, which
## .dxfield() turns into 35 NA gradient cells of 156 per slice), so a predicted
## mean there can hit an *interior* gap as well as leave the rectangle. The
## clamp only ever covered the rectangle -- see dev/code_notes.org,
## "Bounding-box clamp" -> "Scope" -- so the guarantee it does provide can only
## be tested on a field without gaps.
.fill_cov_gaps <- function(x) {
  for (t in seq_len(dim(x)[3])) {
    M <- x[, , t]
    for (it in 1:50) {
      na <- is.na(M)
      if (!any(na)) break
      sh <- function(di, dj) {
        o <- matrix(NA_real_, nrow(M), ncol(M))
        si <- seq_len(nrow(M)) + di
        sj <- seq_len(ncol(M)) + dj
        vi <- si >= 1 & si <= nrow(M)
        vj <- sj >= 1 & sj <= ncol(M)
        o[which(vi), which(vj)] <- M[si[vi], sj[vj]]
        o
      }
      nb <- list(sh(1, 0), sh(-1, 0), sh(0, 1), sh(0, -1))
      num <- Reduce(`+`, lapply(nb, function(z) { z[is.na(z)] <- 0; z }))
      den <- Reduce(`+`, lapply(nb, function(z) !is.na(z)))
      upd <- na & den > 0
      M[upd] <- num[upd] / den[upd]
    }
    x[, , t] <- M
  }
  x
}

.boundary_objs_at <- function(key, fill) {
  .cached(key, {
    sim <- small_sim()
    dat <- sim$dat
    if (fill) dat$cov <- lapply(dat$cov, .fill_cov_gaps)
    par <- sim$par
    par$alpha[] <- 10 * sim$par_true$alpha
    lapply(c(none = "none", clamp = "clamp"), function(b) {
      conf <- sim$conf
      conf$kf_boundary <- b
      suppressWarnings(suppressMessages(
        admove(dat, conf, par, sim$map, run = FALSE, verbose = FALSE)
      ))$obj
    })
  })
}

## land mask kept: means can reach an interior gap
.boundary_objs <- function() .boundary_objs_at("boundary_objs", FALSE)

## gaps filled by neighbour averaging: the only way out is past the edge
.boundary_objs_nogap <- function() .boundary_objs_at("boundary_objs_nogap", TRUE)


test_that("clamping does not change the likelihood while tracks stay inside", {

  objs <- .boundary_objs()
  sim <- small_sim()

  ## the simulated taxis keeps all predicted tracks inside the field
  p <- objs$clamp$par
  p[names(p) == "alpha"] <- p[names(p) == "alpha"] / 10

  expect_true(is.finite(objs$none$fn(p)))
  expect_equal(objs$clamp$fn(p), objs$none$fn(p))
  expect_equal(objs$clamp$gr(p), objs$none$gr(p))

  ids <- names(split(sim$dat$tags, sim$dat$tags$id))
  expect_length(admove:::.boundary_tags(objs$clamp, p, ids), 0)
})


test_that("clamping keeps the likelihood finite when tracks leave the field", {

  ## no interior gaps: leaving the rectangle is the only way out, which is
  ## exactly what the clamp is for
  objs <- .boundary_objs_nogap()
  sim <- small_sim()
  p <- objs$clamp$par

  expect_true(is.nan(objs$none$fn(p)))
  expect_true(is.finite(objs$clamp$fn(p)))
  expect_true(all(is.finite(objs$clamp$gr(p))))

  ## the tracks held at the edge are reported, and only with clamping
  ids <- names(split(sim$dat$tags, sim$dat$tags$id))
  held <- admove:::.boundary_tags(objs$clamp, p, ids)
  expect_gt(length(held), 0)
  expect_true(all(held %in% ids))
  expect_length(admove:::.boundary_tags(objs$none, p, ids), 0)

  expect_match(admove:::.boundary_warning(held), "held at its edge")
})


test_that("clamping does not rescue a mean that reaches an interior gap", {

  ## .dxfield() leaves NA where the covariate is NA (it used to fill those with
  ## zero, which silently gave a mean over a gap no taxis drift at all). The
  ## NaN therefore appears in moveT, i.e. *before* clamp_xy() runs, and
  ## clamping a NaN gives a NaN. Only a boundary correction field (dat$bnd,
  ## make_bnd_field(); on the bound2 branch, not on dev) handles this case.
  gap <- .boundary_objs()
  p <- gap$clamp$par

  expect_true(is.nan(gap$none$fn(p)))
  expect_true(is.nan(gap$clamp$fn(p)))
  expect_length(admove:::.boundary_tags(gap$clamp, p,
                                        names(split(small_sim()$dat$tags,
                                                    small_sim()$dat$tags$id))), 0)

  ## the gaps are the whole difference: same tags, same parameters, same
  ## kf_boundary -- filling them makes the very same objective finite
  expect_true(is.finite(.boundary_objs_nogap()$clamp$fn(p)))
})


test_that("admove stores the tags held at the edge of the field", {

  expect_identical(small_fit()$boundary_tags, character(0))

  ## from the report when there is one
  ids <- c("a", "b", "c")
  expect_identical(admove:::.boundary_ids(c(0, 0.5, NaN), ids), "b")
  expect_identical(admove:::.boundary_ids(NULL, ids), character(0))
})
