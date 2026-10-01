## CTMC drift scheme "sg" (Scharfetter-Gummel). See dev/code_notes.org,
## "Generator drift scheme".

## Generator of a 4 x 4 grid with spatially constant drift and diffusion,
## through the likelihood's own .ctmc_generator().
.sg_generator <- function(scheme, v, D, taxis = TRUE) {
  grid <- create_grid(xrange = c(0, 1), yrange = c(0, 1), cellsize = 0.25,
                      verbose = FALSE)
  nextTo <- get_neighbours(grid)
  next_dist <- rep(0.25, 4)
  funcs <- list(tax = function(xy, t) matrix(v, nrow(xy), 2, byrow = TRUE),
                dif = function(xy, t) log(D),
                adv = function(xy, t) matrix(0, nrow(xy), 2))
  conf <- list(use_taxis = taxis, use_advection = !taxis, drift_scheme = scheme)
  if (!taxis) {
    funcs$adv <- funcs$tax
    funcs$tax <- function(xy, t) matrix(0, nrow(xy), 2)
  }
  ctx <- admove:::.sim_ctmc_ctx(conf, funcs, kappa = 1, grid$xygrid, nextTo,
                                next_dist, "expm")
  list(Q = as.matrix(admove:::.ctmc_generator(ctx, 0, NULL)), nextTo = nextTo,
       grid = grid)
}

## rate from cell i to its neighbour in column `dir` of nextTo, for all cells
## that have that neighbour
.sg_rates <- function(g, dir) {
  ind <- which(!is.na(g$nextTo[, dir]))
  g$Q[cbind(ind, g$nextTo[ind, dir])]
}


test_that("Scharfetter-Gummel rates are non-negative, mass-conserving and exact in the net drift", {

  h <- 0.25
  D <- 0.01
  for (v in list(c(0, 0), c(0.02, -0.01), c(0.4, 0.3), c(-3, 2))) {
    g <- .sg_generator("sg", v, D)
    off <- g$Q
    diag(off) <- NA
    expect_true(all(off >= 0, na.rm = TRUE))
    expect_equal(unname(rowSums(g$Q)), rep(0, nrow(g$Q)), tolerance = 1e-12)
    ## columns of nextTo: 2 top (+y), 3 down (-y), 4 left (-x), 5 right (+x);
    ## an interior cell has all four, and the net rate (with minus against the
    ## drift) is v / h, exactly as for upwind and central
    i <- which(rowSums(is.na(g$nextTo)) == 0)[1]
    r <- function(dir) g$Q[i, g$nextTo[i, dir]]
    expect_equal(r(5) - r(4), v[1] / h, tolerance = 1e-10)
    expect_equal(r(2) - r(3), v[2] / h, tolerance = 1e-10)
  }
})


test_that("Scharfetter-Gummel is central at small and upwind at large grid-Peclet", {

  ## Pe = v h / D = 1e-3: central to O(Pe^2)
  small_sg <- .sg_generator("sg", c(0.04, -0.04), 10)
  small_c <- .sg_generator("central", c(0.04, -0.04), 10)
  expect_equal(small_sg$Q, small_c$Q, tolerance = 1e-6)

  ## Pe = 400: upwind to O(1/Pe). Upwind adds the diffusion D/h^2 on both
  ## sides on top of v/h; in SG it is part of the flux, so the rate against
  ## the drift vanishes exponentially rather than staying at D/h^2.
  big_sg <- .sg_generator("sg", c(16, -16), 0.01)
  big_u <- .sg_generator("upwind", c(16, -16), 0.01)
  ## with the drift (right, v_x > 0; down, v_y < 0): upwind's v/h + D/h^2
  ## against SG's v/h, a relative difference of 1/Pe
  expect_equal(.sg_rates(big_sg, 5), .sg_rates(big_u, 5), tolerance = 1.5 / 400)
  expect_equal(.sg_rates(big_sg, 3), .sg_rates(big_u, 3), tolerance = 1.5 / 400)
  expect_equal(.sg_rates(big_sg, 5), rep(16 / 0.25, length(.sg_rates(big_sg, 5))),
               tolerance = 1e-10)
  ## against the drift: exponentially small instead of D/h^2
  expect_lt(max(.sg_rates(big_sg, 4), .sg_rates(big_sg, 2)), 1e-10)

  ## at no drift all three are the pure diffusion generator
  z <- lapply(c("sg", "central", "upwind"), function(s) .sg_generator(s, c(0, 0), 0.01)$Q)
  expect_equal(z[[1]], z[[2]], tolerance = 1e-12)
  expect_equal(z[[1]], z[[3]], tolerance = 1e-12)

  ## the advection drift goes through the same rates as the taxis drift
  adv <- .sg_generator("sg", c(0.4, 0.3), 0.01, taxis = FALSE)
  tax <- .sg_generator("sg", c(0.4, 0.3), 0.01, taxis = TRUE)
  expect_equal(adv$Q, tax$Q, tolerance = 1e-12)
})


test_that("calc_mstar() builds the same Scharfetter-Gummel generator as the likelihood", {

  v <- c(0.4, -0.3)
  D <- 0.02
  g <- .sg_generator("sg", v, D)
  nc <- nrow(g$Q)
  fit <- list(dat = list(pred = list(grid = g$grid, time = c(0, 1))),
              conf = list(drift_scheme = "sg", use_taxis = TRUE,
                          use_advection = FALSE),
              pred = list(hTdx = matrix(v[1], nc, 1), hTdy = matrix(v[2], nc, 1),
                          hD = matrix(log(D), nc, 1)))
  m <- as.matrix(admove:::calc_mstar(fit)[[1]])
  expect_equal(unname(m), unname(g$Q), tolerance = 1e-12)
})


test_that("the Scharfetter-Gummel likelihood is smooth: AD gradient = finite differences, also at zero drift", {

  grid <- create_grid(cellsize = 0.25, verbose = FALSE)
  cov <- list(cov1 = withr::with_seed(7, sim_cov(grid, nt = 2)))
  tags <- data.frame(id = c("a", "a", "b", "b", "c", "c"),
                     t = c(0, 0.5, 0, 0.6, 0, 0.8),
                     x = c(0.3, 0.6, 0.3, 0.4, 0.6, 0.8),
                     y = c(0.3, 0.4, 0.3, 0.8, 0.6, 0.3))
  tags <- suppressMessages(prep_ctags(tags, names = c(t = "t", x = "x", y = "y", id = "id"),
                                      sref = sref(grid), verbose = FALSE))
  dat <- suppressWarnings(suppressMessages(
    setup_data(grid = grid, cov = cov, tags = tags, trange = c(0, 1), verbose = FALSE)))
  conf <- default_conf(dat, verbose = FALSE)
  conf$engine <- "ctmc"
  conf$drift_scheme <- "sg"
  par <- default_par(dat, conf, verbose = FALSE)
  map <- default_map(dat, conf, par)
  obj <- suppressWarnings(suppressMessages(
    admove(dat, conf, par, map, run = FALSE, verbose = FALSE)))$obj

  fd <- function(p, h = 1e-5) {
    vapply(seq_along(p), function(j) {
      pp <- pm <- p
      pp[j] <- p[j] + h
      pm[j] <- p[j] - h
      (obj$fn(pp) - obj$fn(pm)) / (2 * h)
    }, numeric(1))
  }

  ## the default start: alpha = 0, so the drift is exactly 0 in every cell,
  ## where upwind has its kink (AD = 2 x FD there) and x / expm1(x) is 0/0
  p0 <- obj$par
  expect_true(is.finite(obj$fn(p0)))
  expect_equal(as.numeric(obj$gr(p0)), fd(p0), tolerance = 1e-5)

  ## and away from it
  p1 <- p0
  p1[names(p1) == "alpha"] <- c(0.4, -0.3)[seq_len(sum(names(p1) == "alpha"))]
  expect_equal(as.numeric(obj$gr(p1)), fd(p1), tolerance = 1e-5)
})


test_that("check_conf() accepts drift_scheme = \"sg\" and rejects unknown schemes", {

  dat <- skjepo$sim$dat
  conf <- default_conf(dat, verbose = FALSE)
  conf$drift_scheme <- "sg"
  expect_identical(check_conf(conf, dat, verbose = FALSE)$drift_scheme, "sg")
  conf$drift_scheme <- "exponential"
  expect_error(check_conf(conf, dat, verbose = FALSE), "drift_scheme")
})
