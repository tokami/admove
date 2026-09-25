## CTMC engine: slice-cached generators, aligned lattice and release events.
##
## The constants below were recorded on the implementation that preceded this
## work (branch ctmc_events at 25cfe44, per-tag loop, generator rebuilt every
## step). They pin the *model*, so every later step has to reproduce them under
## conf$ctmc_groups = "off". Two steps are allowed to move them and neither may
## move them silently:
##   - the generator cache: (D*dt)/h^2 and dt*(D/h^2) round differently, so the
##     agreement is ~1e-12 relative, not bit-identical. Never use identical().
##   - the aligned lattice ("align"/"auto"): a genuinely different, finer
##     discretisation. Only "off" is held to these numbers.
##
## They are recorded at the defaults of 25cfe44, which the fixtures set
## explicitly: beta = 0, and for ctmc_pin_obj() also logKappa and the knots.
## logKappa is fixed rather than estimated and the knots define the smooths, so
## both are part of the model. Later changes to the defaults (the diffusion
## start for beta, the kappa scaling, the default knot placement) must not move
## the pins, and would if the fixtures took them from default_par() and
## setup_data() as they are.


## ---------------------------------------------------------------- fixtures

## Mirrors cand_fit_obj() in test-candidates.R, which is the smallest CTMC tape
## in the suite. All three ids are conventional tags released in the same cell
## at the same time, so they also form one release event.
ctmc_pin_obj <- function(w = 0.3, t2 = 0.6, groups = "off") {

  grid <- create_grid(cellsize = 0.25, verbose = FALSE)
  ## seeded: at par = 0 the preference is flat, so fn does not depend on the
  ## covariate but d(nll)/d(alpha) does -- an unseeded field pins fn and leaves
  ## the gradient constants flaky
  cov <- list(cov1 = withr::with_seed(7, sim_cov(grid, nt = 2)))

  x0 <- 0.3; y0 <- 0.3
  x1 <- 0.6; y1 <- 0.4
  x2 <- 0.4; y2 <- 0.8

  tags <- data.frame(
    id    = c("a", "a", "b", "b", "m", "m", "m"),
    t     = c(0, 0.5, 0, t2, 0, 0.5, t2),
    x     = c(x0, x1, x0, x2, x0, x1, x2),
    y     = c(y0, y1, y0, y2, y0, y1, y2),
    event = c(1, 2, 1, 2, 1, 2, 2),
    prob  = c(1, 1, 1, 1, 1, w, 1 - w))
  nms <- c(t = "t", x = "x", y = "y", id = "id", event = "event", prob = "prob")
  tags <- suppressMessages(prep_ctags(tags, names = nms, sref = sref(grid),
                                      verbose = FALSE))

  dat <- suppressWarnings(suppressMessages(
    setup_data(grid = grid, cov = cov, tags = tags, trange = c(0, 1),
               verbose = FALSE)))
  dat$knots_tax[] <- c(21.1571676734481, 24.529971476953499, 27.526040460835198)
  dat$knots_dif[] <- 24.529971476953499
  conf <- default_conf(dat, verbose = FALSE)
  conf$engine <- 2
  conf$ctmc_groups <- groups
  par <- default_par(dat, conf, verbose = FALSE)
  par$beta[] <- 0
  ## the kappa scaling depends on the tag positions, so one value per t2
  log_kappa <- c("0.6" = -0.89422860017923178, "0.85" = -1.2425352944474477)
  par$logKappa[] <- log_kappa[[as.character(t2)]]
  map <- default_map(dat, conf, par)

  suppressWarnings(suppressMessages(
    admove(dat, conf, par, map, run = FALSE, verbose = FALSE)))$obj
}


## A simulated mixed-type set: 2 archival (updating) + 5 conventional (not),
## i.e. the fixture that has to exercise both passes at once. Same seed and
## shape as the CTMC fixtures in test-conf.R and test-tag-pred.R.
ctmc_pin_sim_obj <- function(method = "expav", groups = "off") {

  grid <- create_grid(xrange = c(0, 1), yrange = c(0, 1), cellsize = 0.25,
                      verbose = FALSE)
  sim <- withr::with_seed(1, suppressMessages(suppressWarnings(
    sim_data(grid = grid, n_dtags = 2, n_ctags = 5, trange = c(0, 1),
             verbose = FALSE)
  )))
  conf <- sim$conf
  conf$engine <- 2
  conf$ctmc_method <- method
  conf$ctmc_groups <- groups
  par <- sim$par
  par$beta[] <- 0

  suppressMessages(suppressWarnings(
    admove(sim$dat, conf, par, sim$map, run = FALSE, verbose = FALSE)))$obj
}


## ------------------------------------------------------------ pinned model

test_that("the CTMC likelihood is unchanged (candidates fixture)", {

  obj <- ctmc_pin_obj(w = 0.3, t2 = 0.6)
  p <- obj$par

  expect_equal(obj$fn(p), 8.32032525818662, tolerance = 1e-12)
  expect_equal(as.vector(obj$gr(p)),
               c(-0.620482580323745, 0.664255723232523, -0.0144507163943639),
               tolerance = 1e-10)
  expect_equal(obj$report(p)$loglik_tags,
               c(-2.77259581563871, -2.77408863315109, -2.77364055385788),
               tolerance = 1e-12)
})


test_that("the CTMC likelihood is unchanged (candidates at differing times)", {

  obj <- ctmc_pin_obj(w = 0.45, t2 = 0.85)
  p <- obj$par

  expect_equal(obj$fn(p), 8.31799955696118, tolerance = 1e-12)
  expect_equal(as.vector(obj$gr(p)),
               c(-0.464465647245129, 0.455685878393697, -0.00186065698509374),
               tolerance = 1e-10)
  expect_equal(obj$report(p)$loglik_tags,
               c(-2.77259581563871, -2.77273246053478, -2.7726709680209),
               tolerance = 1e-12)
})


test_that("the CTMC likelihood is unchanged (simulated mixed tag types)", {

  obj <- ctmc_pin_sim_obj()
  p <- obj$par

  expect_equal(obj$fn(p), 58.128093781311, tolerance = 1e-12)
  expect_equal(as.vector(obj$gr(p)),
               c(-2.57697505466118, 0.579053795510462,
                 0.166084137511498, 0.19929172432785),
               tolerance = 1e-10)
  ## element-wise, never the sum: a sum hides an index permutation completely
  expect_equal(obj$report(p)$loglik_tags,
               c(-2.77258878789954, -2.77226417149026, -2.77258878789955,
                 -2.77226417149026, -2.77180522175879,
                 -22.1374465248557, -22.1291352619326),
               tolerance = 1e-12)
})


## --------------------------------------------------- the generator cache key

test_that("habi$slice() changes exactly at the covariate and spline breaks", {

  grid <- create_grid(cellsize = 0.25, verbose = FALSE)
  cov <- list(cov1 = sim_cov(grid, nt = 4))
  dat <- suppressWarnings(suppressMessages(
    setup_data(grid = grid, cov = cov, trange = c(0, 1), verbose = FALSE)))
  conf <- default_conf(dat, verbose = FALSE)
  par <- default_par(dat, conf, verbose = FALSE)

  habi <- admove:::.build_habi(dat, conf, par, per = dat$period)$habi
  breaks <- sort(unique(c(unlist(dat$time_cov), unlist(dat$time_spline))))

  ## the key is constant between consecutive breaks and differs across one
  tt <- sort(c(breaks, breaks + 1e-8, breaks - 1e-8,
               seq(min(breaks), max(breaks), length.out = 37)))
  tt <- tt[tt >= min(breaks) & tt <= max(breaks)]
  keys <- vapply(tt, function(z) paste(habi$tax$slice(z), collapse = ","), "")
  ints <- findInterval(tt, breaks)

  ## one key per interval, one interval per key
  expect_equal(length(unique(paste(keys, ints))), length(unique(ints)))
  expect_equal(length(unique(keys)), length(unique(ints)))

  ## tax and dif are built with the same time arguments, so one of them serves
  ## the generator key; .ctmc_slice_key() relies on this
  for (z in tt) {
    expect_identical(habi$tax$slice(z), habi$dif$slice(z))
  }
})


## A current that changes at t = 0, 1, ..., 5 and a covariate that is constant
## in time, given either as 6 identical slices (its breaks happen to cover the
## current's) or as 1 slice (only the advection part of the key and the breaks
## can). Before advection entered .ctmc_slice_key() and .ctmc_breaks(), the
## single-slice version silently reused the t = 0 generator throughout.
ctmc_adv_obj <- function(nt_cov) {

  grid <- create_grid(cellsize = 0.25, verbose = FALSE)
  base <- withr::with_seed(7, sim_cov(grid, nt = 6))
  cov <- withr::with_seed(7, sim_cov(grid, nt = nt_cov))
  for (k in seq_len(nt_cov)) cov[, , k] <- base[, , 1]
  u <- withr::with_seed(3, sim_cov(grid, nt = 6))
  v <- withr::with_seed(4, sim_cov(grid, nt = 6))
  u[] <- u[] - mean(u)
  v[] <- v[] - mean(v)

  tags <- data.frame(
    id = c("a", "a", "b", "b", "c", "c"),
    t  = c(0, 4.6, 0, 3.4, 0, 4.2),
    x  = c(0.3, 0.6, 0.3, 0.4, 0.3, 0.8),
    y  = c(0.3, 0.4, 0.3, 0.8, 0.3, 0.6))
  nms <- c(t = "t", x = "x", y = "y", id = "id")
  tags <- suppressMessages(prep_ctags(tags, names = nms, sref = sref(grid),
                                      verbose = FALSE))
  dat <- suppressWarnings(suppressMessages(
    setup_data(grid = grid, cov = list(cov1 = cov), tags = tags,
               adv = list(cur = prep_adv(u, v, units = NULL)),
               trange = c(0, 5), verbose = FALSE)))
  dat$min_dt <- 0.4
  conf <- default_conf(dat, verbose = FALSE)
  conf$engine <- 2
  par <- default_par(dat, conf, verbose = FALSE)
  map <- default_map(dat, conf, par)

  suppressWarnings(suppressMessages(
    admove(dat, conf, par, map, run = FALSE, verbose = FALSE)))$obj
}


test_that("the CTMC generator follows the advection time slices", {

  a <- ctmc_adv_obj(6)
  b <- ctmc_adv_obj(1)

  ## away from par = 0, so that the current (gamma) and the taxis act
  p <- a$par
  p[] <- seq(0.3, 1.2, length.out = length(p))

  expect_equal(b$fn(p), a$fn(p), tolerance = 1e-8)
  expect_equal(b$report(p)$ctmc_ngen, a$report(p)$ctmc_ngen)
})


test_that(".ctmc_breaks() includes advection field times and season starts", {

  d <- list(tags = data.frame(t = c(0, 2.5)),
            time_cov = list(0), time_spline = list(0),
            use_advection = TRUE, time_adv = list(c(0, 1.5), c(0, 1.5)),
            n_seasons_adv = 2L, period = 1)

  expect_equal(admove:::.ctmc_breaks(d), seq(0, 3.5, by = 0.5))

  d$n_seasons_adv <- 1L
  expect_equal(admove:::.ctmc_breaks(d), c(0, 1.5))

  d$use_advection <- FALSE
  expect_equal(admove:::.ctmc_breaks(d), 0)
})


## ------------------------------------------------------- the aligned lattice

## Tags spanning several covariate slices, so the slice boundaries fall inside
## their time range and alignment actually changes the lattice. A fixture whose
## tags sit inside one slice makes every test below vacuous.
ctmc_span_obj <- function(groups, min_dt) {

  grid <- create_grid(cellsize = 0.25, verbose = FALSE)
  cov <- list(cov1 = withr::with_seed(7, sim_cov(grid, nt = 6)))
  tags <- data.frame(
    id = c("a", "a", "b", "b"),
    t  = c(0, 3.6, 0, 2.4),
    x  = c(0.3, 0.6, 0.3, 0.4),
    y  = c(0.3, 0.4, 0.3, 0.8))
  nms <- c(t = "t", x = "x", y = "y", id = "id")
  tags <- suppressMessages(prep_ctags(tags, names = nms, sref = sref(grid),
                                      verbose = FALSE))
  dat <- suppressWarnings(suppressMessages(
    setup_data(grid = grid, cov = cov, tags = tags, trange = c(0, 5),
               verbose = FALSE)))
  dat$min_dt <- min_dt
  conf <- default_conf(dat, verbose = FALSE)
  conf$engine <- 2
  conf$ctmc_groups <- groups
  par <- default_par(dat, conf, verbose = FALSE)
  map <- default_map(dat, conf, par)

  list(obj = suppressWarnings(suppressMessages(
         admove(dat, conf, par, map, run = FALSE, verbose = FALSE)))$obj,
       dat = dat)
}


test_that(".ctmc_breaks() returns the covariate and spline slice boundaries", {

  z <- ctmc_span_obj("align", 0.4)
  expect_equal(admove:::.ctmc_breaks(z$dat),
               sort(unique(c(unlist(z$dat$time_cov),
                             unlist(z$dat$time_spline)))))

  ## callable on both tag representations: nll() is handed the split() list,
  ## setup_data() stores a data.frame
  d2 <- z$dat
  d2$tags <- split(d2$tags, d2$tags$id)
  expect_equal(admove:::.ctmc_breaks(d2), admove:::.ctmc_breaks(z$dat))
})


test_that("the aligned lattice is a refinement, not a different model", {

  ## Evaluate AWAY from par = 0. At the origin there is no taxis and a single
  ## diffusion knot gives D = exp(0) = 1, so the generator is the same matrix in
  ## every slice: alignment provably cannot change the answer and the test would
  ## pass on a gap of pure expAv truncation noise (~1e-8).
  p <- c(1.2, -0.8, 0.4)

  gaps <- vapply(c(0.8, 0.4, 0.2), function(md) {
    a <- ctmc_span_obj("align", md)
    b <- ctmc_span_obj("off", md)
    abs(a$obj$fn(p) - b$obj$fn(p))
  }, numeric(1))

  ## the two lattices agree in the limit: halving min_dt must shrink the gap.
  ## "the difference is small" would not distinguish a refinement from a bug
  expect_true(all(diff(gaps) < 0))
  expect_lt(gaps[3], gaps[1] / 100)
})


test_that("conf$ctmc_groups is validated and defaults to auto", {

  dat <- small_sim()$dat
  expect_equal(default_conf(dat, verbose = FALSE)$ctmc_groups, "auto")
  expect_error(check_conf(list(ctmc_groups = "yes"), dat, verbose = FALSE),
               "ctmc_groups")

  ## NULL-safe: a conf built before the setting existed reaches nll() without
  ## it, and admove() never calls check_conf()
  expect_equal(admove:::.ctmc_mode(list()), "auto")
  expect_equal(admove:::.ctmc_mode(list(ctmc_groups = "off")), "off")
})


## ----------------------------------------------------------- the partition

## Tags that genuinely share a release cell and time. Built by hand on purpose:
## sim_data() scatters releases within a cell, so a simulated set forms no
## group at all and every assertion below would hold vacuously.
ctmc_event_obj <- function(n = 4, groups = "auto", do_update = NULL,
                           tag_types = NULL, obs_var = NULL, t2 = 0.6) {

  grid <- create_grid(cellsize = 0.25, verbose = FALSE)
  cov <- list(cov1 = withr::with_seed(7, sim_cov(grid, nt = 2)))

  tags <- data.frame(
    id = rep(paste0("t", seq_len(n)), each = 2),
    t  = rep(c(0, t2), n),
    x  = as.numeric(rbind(0.3, seq(0.4, 0.7, length.out = n))),
    y  = as.numeric(rbind(0.3, seq(0.4, 0.8, length.out = n))))
  nms <- c(t = "t", x = "x", y = "y", id = "id")
  tags <- suppressMessages(prep_ctags(tags, names = nms, sref = sref(grid),
                                      verbose = FALSE))

  dat <- suppressWarnings(suppressMessages(
    setup_data(grid = grid, cov = cov, tags = tags, trange = c(0, 1),
               verbose = FALSE)))
  conf <- default_conf(dat, verbose = FALSE)
  conf$engine <- 2
  conf$ctmc_groups <- groups
  if (!is.null(do_update)) conf$do_update <- do_update
  if (!is.null(obs_var)) conf$obs_var_type <- obs_var
  ## via the constructor, not a hand-rolled list: admove() validates the class,
  ## and a fixture that bypasses it would not exercise the real path
  if (!is.null(tag_types)) {
    dat <- set_release_events(dat, tag_types = tag_types, verbose = FALSE)
  }
  par <- default_par(dat, conf, verbose = FALSE)
  map <- default_map(dat, conf, par)

  suppressWarnings(suppressMessages(
    admove(dat, conf, par, map, run = FALSE, verbose = FALSE)))$obj
}


test_that("the fixture really does form a release event", {

  ## guards every partition and equivalence test below: with ngrouped == 0 they
  ## would compare the per-tag pass with itself and pass forever
  o <- ctmc_event_obj(n = 4)
  r <- o$report(o$par)
  expect_equal(r$ctmc_ngroup, 1)
  expect_equal(r$ctmc_ngrouped, 4)
})


test_that("grouping follows conf$do_update, not the tags$update column", {

  ## check_tags() writes tags$update and nll() has never read it. A rule built
  ## on that column would ignore the configuration and group updating tags.
  expect_equal(ctmc_event_obj(do_update = c(TRUE, TRUE, FALSE))$
                 report(ctmc_event_obj()$par)$ctmc_ngrouped, 4)

  o_upd <- ctmc_event_obj(do_update = c(TRUE, TRUE, TRUE))
  expect_equal(o_upd$report(o_upd$par)$ctmc_ngrouped, 0)

  o_all <- ctmc_event_obj(do_update = c(FALSE, FALSE, FALSE))
  expect_equal(o_all$report(o_all$par)$ctmc_ngrouped, 4)
})


test_that("grouping is gated by the mode and by the configured tag types", {

  for (g in c("align", "off")) {
    o <- ctmc_event_obj(groups = g)
    expect_equal(o$report(o$par)$ctmc_ngrouped, 0)
  }

  o_c <- ctmc_event_obj(tag_types = "c")
  expect_equal(o_c$report(o_c$par)$ctmc_ngrouped, 4)

  o_d <- ctmc_event_obj(tag_types = "d")
  expect_equal(o_d$report(o_d$par)$ctmc_ngrouped, 0)
})


## ------------------------------------------------- grouped == ungrouped

## A release event shares one forward pass; the per-tag pass is the reference.
## Compare against "align", not "off": both then use the same lattice, so the
## comparison isolates the grouping from the lattice change.
##
## Perturbations are built from length(par): obs_var_type != "none" frees
## logSdO, so a fixed-length offset would recycle and perturb something other
## than what it names.
.ctmc_pars <- function(obj) {
  p <- obj$par
  n <- length(p)
  list(p,
       p + rep_len(c(0.9, -0.6, 0.35), n),
       p + rep_len(c(-1.1, 0.8, -0.4), n))
}

.expect_group_equiv <- function(lab, ...) {
  a <- ctmc_event_obj(groups = "auto", ...)
  b <- ctmc_event_obj(groups = "align", ...)

  ra <- a$report(a$par)
  ## without this the comparison degenerates to the per-tag pass against
  ## itself, which passes for the wrong reason and keeps passing
  expect_gt(ra$ctmc_ngrouped, 0)

  ## one pass instead of one per member: deterministic, unlike wall time
  expect_lt(ra$ctmc_nstep, b$report(b$par)$ctmc_nstep)

  for (p in .ctmc_pars(a)) {
    expect_equal(a$fn(p), b$fn(p), tolerance = 1e-10, info = lab)
    expect_equal(as.vector(a$gr(p)), as.vector(b$gr(p)),
                 tolerance = 1e-8, info = lab)
    ## element-wise: a sum would hide a term landing on the wrong tag
    expect_equal(a$report(p)$loglik_tags, b$report(p)$loglik_tags,
                 tolerance = 1e-10, info = lab)
  }

  ## a cache bound to the wrong tape is right at the taping point and wrong
  ## away from it, so this has to be checked after a retape as well
  a$retape()
  p <- .ctmc_pars(a)[[2]]
  expect_equal(a$fn(p), b$fn(p), tolerance = 1e-10, info = lab)
}


test_that("a release event gives the same likelihood as one pass per tag", {
  .expect_group_equiv("plain", n = 4)
  .expect_group_equiv("six members", n = 6)
})


test_that("grouping preserves the observation-error variants", {
  .expect_group_equiv("all", obs_var = c("none", "none", "all"))
  .expect_group_equiv("all_but_last", obs_var = c("none", "none", "all_but_last"))
})


test_that("an ambiguous final event stays a per-member mixture when grouped", {

  ## Three ids sharing a release, the third with two candidate recaptures.
  ## The exact-mixture identity is the strong form of the check: it fails if
  ## the mixture is pooled across members or a candidate is counted twice,
  ## which equality-to-the-other-path alone might not catch.
  w <- 0.3
  obj <- ctmc_pin_obj(w = w, t2 = 0.6, groups = "auto")
  r <- obj$report(obj$par)
  expect_gt(r$ctmc_ngrouped, 0)

  ll <- r$loglik_tags
  expect_length(ll, 3L)
  expect_equal(ll[3], log(w * exp(ll[1]) + (1 - w) * exp(ll[2])),
               tolerance = 1e-10)
})


## ------------------------------------------------------ set_release_events

test_that("set_release_events stores a policy and validates its arguments", {

  dat <- small_sim()$dat

  re <- set_release_events(dat, verbose = FALSE)$release_events
  expect_s3_class(re, "admove_release_events")
  expect_true(re$exact)
  expect_null(re$keys)            ## exact grouping needs no map
  expect_equal(re$tag_types, "c")

  expect_error(set_release_events(dat, tag_types = "x"), "tag_types")
  expect_error(set_release_events(dat, t_tol = -1), "t_tol")
  expect_error(set_release_events(dat, dist_tol = NA), "dist_tol")
})


test_that("a tolerance merges tags, warns, and is honoured by the engine", {

  ## two releases a short distance and time apart: exact grouping keeps them
  ## separate, a tolerance merges them
  grid <- create_grid(cellsize = 0.25, verbose = FALSE)
  cov <- list(cov1 = withr::with_seed(7, sim_cov(grid, nt = 2)))
  tags <- data.frame(
    id = rep(c("a", "b"), each = 2),
    t  = c(0, 0.6, 0.02, 0.6),
    x  = c(0.30, 0.55, 0.32, 0.65),
    y  = c(0.30, 0.45, 0.31, 0.70))
  nms <- c(t = "t", x = "x", y = "y", id = "id")
  tags <- suppressMessages(prep_ctags(tags, names = nms, sref = sref(grid),
                                      verbose = FALSE))
  dat <- suppressWarnings(suppressMessages(
    setup_data(grid = grid, cov = cov, tags = tags, trange = c(0, 1),
               verbose = FALSE)))

  ## different release times, so exact grouping forms nothing
  fit_obj <- function(d) {
    conf <- default_conf(d, verbose = FALSE)
    conf$engine <- 2
    par <- default_par(d, conf, verbose = FALSE)
    map <- default_map(d, conf, par)
    suppressWarnings(suppressMessages(
      admove(d, conf, par, map, run = FALSE, verbose = FALSE)))$obj
  }
  o0 <- fit_obj(dat)
  expect_equal(o0$report(o0$par)$ctmc_ngrouped, 0)

  ## the approximation is announced, with its size
  expect_warning(dat2 <- set_release_events(dat, t_tol = 0.1, dist_tol = 0.1,
                                            verbose = FALSE),
                 "merged 2 tags")
  expect_false(dat2$release_events$exact)
  expect_length(dat2$release_events$keys, 2L)

  o1 <- fit_obj(dat2)
  expect_equal(o1$report(o1$par)$ctmc_ngrouped, 2)

  ## and it really is an approximation: the likelihood moves
  expect_false(isTRUE(all.equal(o0$fn(o0$par), o1$fn(o1$par))))
})


test_that("admove rejects a malformed release_events object", {

  dat <- small_sim()$dat
  dat$release_events <- list(tag_types = "c")
  conf <- default_conf(dat, verbose = FALSE)
  conf$engine <- 2

  expect_error(
    suppressMessages(admove(dat, conf, run = FALSE, verbose = FALSE)),
    "admove_release_events")
})
