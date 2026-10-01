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
## Two fixtures. small_sim()'s covariate carries the skjepo land mask, and the
## clamp only ever covered the rectangle -- see dev/code_notes.org,
## "Bounding-box clamp" -> "Scope" -- so the guarantee it does provide is
## tested on the field with its gaps filled. Whether a simulated track happens
## to cross a land gap depends on data/skjepo.rda, so the interior-gap case
## has its own deterministic fixture (.gap_obj() below).
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


## A field with a hole in its interior and one tag whose predicted mean the
## taxis drives into it (+x at ~1 per unit time from x = 0.3, steps of 0.1, the
## hole at x 0.5-0.7) well before it could reach the edge. Deterministic on
## purpose: this test used small_sim(), and when data/skjepo.rda was
## regenerated (86de79c, 2026-09-29) its simulated tracks stopped crossing a
## land gap, so the test failed without anything having changed in the clamp.
## Filled (hole = FALSE), the same mean runs on to the edge, which is the case
## the clamp does handle.
.gap_obj <- function(hole, boundary) {
  .cached(paste("gap_obj", hole, boundary), {
    grid <- create_grid(xrange = c(0, 1), yrange = c(0, 1), cellsize = 0.1,
                        verbose = FALSE)
    xc <- seq(0.05, 0.95, by = 0.1)
    m <- outer(xc, xc, function(x, y) 20 + 10 * x)
    if (hole) m[6:7, 4:7] <- NA
    cov <- prep_cov(m, x_centers = xc, y_centers = xc, times = 0, verbose = FALSE)
    tags <- prep_ctags(data.frame(id = "a", t = c(0, 1), x = c(0.3, 0.35), y = 0.5),
                       names = c(t = "t", x = "x", y = "y", id = "id"),
                       verbose = FALSE)
    dat <- suppressWarnings(suppressMessages(
      setup_data(grid = grid, cov = list(cov1 = cov), tags = tags,
                 trange = c(0, 1), knots_tax = matrix(c(21, 25, 29), ncol = 1),
                 knots_dif = matrix(25, ncol = 1), verbose = FALSE)))
    conf <- default_conf(dat, verbose = FALSE)
    conf$engine <- "kf"
    conf$kf_boundary <- boundary
    par <- default_par(dat, conf, verbose = FALSE)
    ## drift = kappa x d(pref)/dz x dz/dx = 1 x 0.1 x 10
    par$logKappa <- 0
    par$alpha[] <- c(0, 0.4, 0.8)
    par$beta[] <- log(0.001)
    map <- default_map(dat, conf, par)
    suppressWarnings(suppressMessages(
      admove(dat, conf, par, map, run = FALSE, verbose = FALSE)))$obj
  })
}


test_that("clamping does not rescue a mean that reaches an interior gap", {

  ## .dxfield() leaves NA where the covariate is NA (it used to fill those with
  ## zero, which silently gave a mean over a gap no taxis drift at all). The
  ## NaN therefore appears in moveT, i.e. *before* clamp_xy() runs, and
  ## clamping a NaN gives a NaN. Only a boundary correction field (dat$bnd,
  ## make_bnd_field(); on the bound2 branch, not on dev) handles this case.
  none <- .gap_obj(TRUE, "none")
  clamp <- .gap_obj(TRUE, "clamp")
  p <- clamp$par

  expect_true(is.nan(none$fn(p)))
  expect_true(is.nan(clamp$fn(p)))
  expect_length(admove:::.boundary_tags(clamp, p, "a"), 0)

  ## the gap is the whole difference: same tag, same parameters, same
  ## kf_boundary -- without the hole the mean reaches the edge instead, where
  ## the clamp holds it and the objective is finite
  filled <- .gap_obj(FALSE, "clamp")
  expect_true(is.nan(.gap_obj(FALSE, "none")$fn(p)))
  expect_true(is.finite(filled$fn(p)))
  expect_identical(admove:::.boundary_tags(filled, p, "a"), "a")
})


test_that("admove stores the tags held at the edge of the field", {

  expect_identical(small_fit()$boundary_tags, character(0))

  ## from the report when there is one
  ids <- c("a", "b", "c")
  expect_identical(admove:::.boundary_ids(c(0, 0.5, NaN), ids), "b")
  expect_identical(admove:::.boundary_ids(NULL, ids), character(0))
})
