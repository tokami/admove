## require(admove); require(testthat)

## Ambiguous final observations: a tag whose last position is known only up to a
## set of candidates with probabilities (e.g. the fishing sets of the vessel
## that recovered it).


make_cand_tags <- function() {
  data.frame(
    id = c("plain", "plain",
           rep("amb", 4)),
    t = c(0, 10,
          0, 10, 10, 12),
    x = c(1, 2,
          1, 2, 5, 7),
    y = c(1, 2,
          1, 2, 5, 7),
    event = c(1, 2,
              1, 2, 2, 2),
    prob = c(1, 1,
             1, 0.5, 0.3, 0.2),
    tag_type = rep("c", 6)
  )
}


## ---------------------------------------------------------------- input routes

test_that("wide input with candidate suffixes reaches the long form", {

  w <- data.frame(
    fish = c("a", "b"),
    rel_t = c(0, 0), rel_x = c(1, 2), rel_y = c(3, 4),
    date1 = c(10, 10), lon1 = c(5, 6), lat1 = c(7, 8), per1 = c(0.4, 1),
    date2 = c(12, NA), lon2 = c(9, NA), lat2 = c(11, NA), per2 = c(0.6, NA)
  )

  out <- prep_ctags(
    w,
    names = c(id = "fish", t0 = "rel_t", x0 = "rel_x", y0 = "rel_y",
              t1 = "date", x1 = "lon", y1 = "lat", p1 = "per"),
    candidates = 1:2, verbose = FALSE)

  ## tag "a" keeps both candidates, tag "b" only the one that is not NA-padded
  expect_equal(nrow(out), 5L)
  expect_equal(out$event, c(1, 2, 2, 1, 2))
  expect_equal(out$prob, c(1, 0.4, 0.6, 1, 1))

  ## the explicit list form must give exactly the same object
  out2 <- prep_ctags(
    w,
    names = list(id = "fish", t0 = "rel_t", x0 = "rel_x", y0 = "rel_y",
                 t1 = paste0("date", 1:2), x1 = paste0("lon", 1:2),
                 y1 = paste0("lat", 1:2), p1 = paste0("per", 1:2)),
    verbose = FALSE)
  expect_equal(out, out2)
})


test_that("a classic wide input is unchanged by the candidate machinery", {

  w <- data.frame(id = c("a", "b"),
                  t0 = c(0, 0), t1 = c(10, 11),
                  x0 = c(1, 2), y0 = c(3, 4),
                  x1 = c(5, 6), y1 = c(7, 8))

  out <- prep_ctags(w, names = c(id = "id", t0 = "t0", t1 = "t1", x0 = "x0",
                                 y0 = "y0", x1 = "x1", y1 = "y1"),
                    verbose = FALSE)

  expect_equal(nrow(out), 4L)
  expect_false("event" %in% colnames(out))
  expect_false("prob" %in% colnames(out))
})


test_that("prep_tags keeps the class of Date time columns through wide2long", {

  w <- data.frame(id = c("a", "b"),
                  t0 = as.Date("2020-01-01") + c(0, 1),
                  t1 = as.Date("2020-03-01") + c(0, 1),
                  x0 = c(1, 2), y0 = c(3, 4), x1 = c(5, 6), y1 = c(7, 8))

  out <- prep_ctags(w, names = c(id = "id", t0 = "t0", t1 = "t1", x0 = "x0",
                                 y0 = "y0", x1 = "x1", y1 = "y1"),
                    date_format = "%Y-%m-%d", verbose = FALSE)

  expect_equal(nrow(out), 4L)
  expect_true(is.numeric(out$t))
})


test_that("add_candidates replaces the final observation", {

  tags <- prep_ctags(
    data.frame(id = c("a", "b"), t0 = c(0, 0), t1 = c(10, 10),
               x0 = c(1, 1), y0 = c(1, 1), x1 = c(2, 2), y1 = c(2, 2)),
    names = c(id = "id", t0 = "t0", t1 = "t1", x0 = "x0", y0 = "y0",
              x1 = "x1", y1 = "y1"), verbose = FALSE)

  sets <- data.frame(key = c("a", "a", "a"),
                     st = c(9, 10, 11), sx = c(2, 3, 4), sy = c(2, 3, 4),
                     sp = c(0.2, 0.5, 0.3))

  out <- add_candidates(tags, sets, by = c(id = "key"),
                        names = c(t = "st", x = "sx", y = "sy", prob = "sp"),
                        verbose = FALSE)

  expect_equal(sum(out$id == "a"), 4L)   ## release + 3 candidates
  expect_equal(sum(out$id == "b"), 2L)   ## untouched
  expect_equal(out$prob[out$id == "a"], c(1, 0.2, 0.5, 0.3))
  expect_s3_class(out, "admove_tags")
})


## ----------------------------------------------------------------- validation

test_that("check_tags keeps candidates that share a time", {

  out <- check_tags(make_cand_tags(), verbose = FALSE)

  ## the two candidates at t = 10 must both survive the zero-time-step filter
  expect_equal(sum(out$id == "amb"), 4L)
  expect_equal(out$event[out$id == "amb"], c(1, 2, 2, 2))
})


test_that("check_tags adds event and prob when they are absent", {

  tg <- make_cand_tags()
  tg <- tg[tg$id == "plain", setdiff(colnames(tg), c("event", "prob"))]

  out <- check_tags(tg, verbose = FALSE)

  expect_equal(out$event, c(1, 2))
  expect_equal(out$prob, c(1, 1))
})


test_that("check_tags rejects candidate probabilities that do not sum to 1", {

  tg <- make_cand_tags()
  tg$prob[tg$id == "amb"][4] <- 0.9

  expect_error(check_tags(tg, verbose = FALSE), "sum to 1")
})


test_that("check_tags rejects ambiguity that is not the final observation", {

  tg <- make_cand_tags()
  tg$id <- "amb"
  tg <- data.frame(
    id = rep("x", 4),
    t = c(0, 5, 5, 10),
    x = c(1, 2, 3, 4), y = c(1, 2, 3, 4),
    event = c(1, 2, 2, 3),
    prob = c(1, 0.5, 0.5, 1),
    tag_type = rep("s", 4))

  expect_error(check_tags(tg, verbose = FALSE), "not the last one")
})


test_that("check_tags rejects an ambiguous release", {

  tg <- data.frame(
    id = rep("x", 3),
    t = c(0, 0, 10),
    x = c(1, 2, 3), y = c(1, 2, 3),
    event = c(1, 1, 2),
    prob = c(0.5, 0.5, 1),
    tag_type = rep("c", 3))

  expect_error(check_tags(tg, verbose = FALSE), "FIRST observation")
})


test_that("check_tags warns about a repeated candidate position", {

  tg <- make_cand_tags()
  tg$x[tg$id == "amb"][3] <- tg$x[tg$id == "amb"][2]
  tg$y[tg$id == "amb"][3] <- tg$y[tg$id == "amb"][2]

  expect_warning(check_tags(tg, verbose = TRUE), "same candidate position twice")
})


test_that("check_tags counts events, not rows, for duplicated ctag ids", {

  ## 3 candidates = 2 events -> no warning
  expect_silent(check_tags(make_cand_tags(), verbose = TRUE))

  ## three genuine observations = 3 events -> warning
  tg <- data.frame(id = rep("x", 3), t = c(0, 5, 10),
                   x = c(1, 2, 3), y = c(1, 2, 3), tag_type = rep("c", 3))
  expect_warning(check_tags(tg, verbose = TRUE), "more than two observations")
})


## ------------------------------------------------------------- the likelihood

## Build an admove_tags object from a plain long data frame, attaching the
## spatial reference of the grid. Doubles as a check that the long input route
## (explicit event / prob columns) survives prep_tags().
as_cand_tags <- function(df, grid) {
  nms <- c(t = "t", x = "x", y = "y", id = "id")
  if ("event" %in% colnames(df)) nms <- c(nms, event = "event", prob = "prob")
  suppressMessages(prep_ctags(df, names = nms, sref = sref(grid),
                              verbose = FALSE))
}


## Three tags evaluated together so they share one parameter vector:
##   "a" : release + candidate 1 only (prob 1)
##   "b" : release + candidate 2 only (prob 1)
##   "m" : release + both candidates, weights w and 1 - w
## The mixture tag must satisfy  L_m = w * L_a + (1 - w) * L_b  exactly.
cand_fit_obj <- function(engine, w = 0.3, t2 = 0.6, cellsize = 0.25) {

  grid <- create_grid(cellsize = cellsize, verbose = FALSE)
  cov <- list(cov1 = sim_cov(grid, nt = 2))

  x0 <- 0.3; y0 <- 0.3
  x1 <- 0.6; y1 <- 0.4      ## candidate 1
  x2 <- 0.4; y2 <- 0.8      ## candidate 2

  tags <- data.frame(
    id    = c("a", "a", "b", "b", "m", "m", "m"),
    t     = c(0, 0.5, 0, t2, 0, 0.5, t2),
    x     = c(x0, x1, x0, x2, x0, x1, x2),
    y     = c(y0, y1, y0, y2, y0, y1, y2),
    event = c(1, 2, 1, 2, 1, 2, 2),
    prob  = c(1, 1, 1, 1, 1, w, 1 - w))
  tags <- as_cand_tags(tags, grid)

  dat <- suppressWarnings(suppressMessages(
    setup_data(grid = grid, cov = cov, tags = tags, trange = c(0, 1),
               verbose = FALSE)))
  conf <- default_conf(dat, verbose = FALSE)
  conf$engine <- engine
  par <- default_par(dat, conf, verbose = FALSE)
  map <- default_map(dat, conf, par)

  suppressWarnings(suppressMessages(
    admove(dat, conf, par, map, run = FALSE, verbose = FALSE)))$obj
}


test_that("the mixture is the exact weighted sum of its components (KF)", {

  w <- 0.3
  obj <- cand_fit_obj(engine = 1, w = w)
  invisible(obj$fn(obj$par))
  ll <- obj$report(obj$par)$loglik_tags   ## ordered by id: a, b, m

  expect_length(ll, 3L)
  expect_equal(ll[3], log(w * exp(ll[1]) + (1 - w) * exp(ll[2])),
               tolerance = 1e-10)
})


test_that("the mixture is the exact weighted sum of its components (CTMC)", {

  w <- 0.3
  obj <- cand_fit_obj(engine = 2, w = w)
  invisible(obj$fn(obj$par))
  ll <- obj$report(obj$par)$loglik_tags

  expect_length(ll, 3L)
  expect_equal(ll[3], log(w * exp(ll[1]) + (1 - w) * exp(ll[2])),
               tolerance = 1e-10)
})


test_that("candidates at different times are each evaluated at their own time", {

  ## candidate 2 sits at a different time from candidate 1; the identity must
  ## still hold, which it only can if the state is propagated to each in turn
  for (eng in 1:2) {
    w <- 0.45
    obj <- cand_fit_obj(engine = eng, w = w, t2 = 0.85)
    invisible(obj$fn(obj$par))
    ll <- obj$report(obj$par)$loglik_tags
    expect_equal(ll[3], log(w * exp(ll[1]) + (1 - w) * exp(ll[2])),
                 tolerance = 1e-10)
  }
})


test_that("a single candidate at prob 1 reproduces a tag with no event columns", {

  mk <- function(with_cols, engine) {
    grid <- create_grid(cellsize = 0.25, verbose = FALSE)
    cov <- list(cov1 = sim_cov(grid, nt = 2))
    tags <- data.frame(id = c("a", "a", "b", "b"),
                       t = c(0, 0.5, 0, 0.7),
                       x = c(0.3, 0.6, 0.2, 0.8),
                       y = c(0.3, 0.4, 0.7, 0.2))
    if (with_cols) {
      tags$event <- c(1, 2, 1, 2)
      tags$prob <- c(1, 1, 1, 1)
    }
    tags <- as_cand_tags(tags, grid)
    dat <- suppressWarnings(suppressMessages(
      setup_data(grid = grid, cov = cov, tags = tags, trange = c(0, 1),
                 verbose = FALSE)))
    conf <- default_conf(dat, verbose = FALSE)
    conf$engine <- engine
    par <- default_par(dat, conf, verbose = FALSE)
    map <- default_map(dat, conf, par)
    obj <- suppressWarnings(suppressMessages(
      admove(dat, conf, par, map, run = FALSE, verbose = FALSE)))$obj
    list(fn = as.numeric(obj$fn(obj$par)), gr = as.numeric(obj$gr(obj$par)))
  }

  for (eng in 1:2) {
    set.seed(1); a <- mk(TRUE, eng)
    set.seed(1); b <- mk(FALSE, eng)
    expect_equal(a$fn, b$fn, tolerance = 1e-12)
    expect_equal(a$gr, b$gr, tolerance = 1e-12)
  }
})


test_that("a far-away candidate does not underflow the objective or gradient", {

  grid <- create_grid(cellsize = 0.25, verbose = FALSE)
  cov <- list(cov1 = sim_cov(grid, nt = 2))

  ## candidate 2 is on the far corner: its density is vanishingly small, which
  ## a naive log(sum(exp(.))) would turn into -Inf
  tags <- data.frame(id = rep("m", 3),
                     t = c(0, 0.02, 0.02),
                     x = c(0.15, 0.2, 0.85),
                     y = c(0.15, 0.2, 0.85),
                     event = c(1, 2, 2),
                     prob = c(1, 0.5, 0.5))
  tags <- as_cand_tags(tags, grid)

  dat <- suppressWarnings(suppressMessages(
    setup_data(grid = grid, cov = cov, tags = tags, trange = c(0, 1),
               verbose = FALSE)))
  conf <- default_conf(dat, verbose = FALSE)
  par <- default_par(dat, conf, verbose = FALSE)
  map <- default_map(dat, conf, par)
  obj <- suppressWarnings(suppressMessages(
    admove(dat, conf, par, map, run = FALSE, verbose = FALSE)))$obj

  expect_true(is.finite(obj$fn(obj$par)))
  expect_true(all(is.finite(obj$gr(obj$par))))
})


## -------------------------------------------------------- dropping candidates

test_that("probabilities are rescaled when a candidate is removed", {

  grid <- create_grid(cellsize = 0.25, verbose = FALSE)

  ## the third candidate lies outside the grid and is dropped; the remaining
  ## two must be rescaled from 0.2 / 0.3 back to summing to 1
  tags <- data.frame(id = rep("m", 4),
                     t = c(0, 0.5, 0.5, 0.5),
                     x = c(0.3, 0.6, 0.4, 9),
                     y = c(0.3, 0.4, 0.8, 9),
                     event = c(1, 2, 2, 2),
                     prob = c(1, 0.2, 0.3, 0.5),
                     tag_type = rep("c", 4))

  out <- suppressMessages(check_tags(tags, grid, verbose = FALSE))

  expect_equal(sum(out$id == "m"), 3L)
  expect_equal(sum(out$prob[out$event == 2]), 1)
  expect_equal(out$prob[out$event == 2], c(0.4, 0.6))
})


test_that("a tag that loses its release is dropped", {

  grid <- create_grid(cellsize = 0.25, verbose = FALSE)

  tags <- data.frame(id = rep("m", 3),
                     t = c(0, 0.5, 0.5),
                     x = c(9, 0.6, 0.4),      ## release outside the grid
                     y = c(9, 0.4, 0.8),
                     event = c(1, 2, 2),
                     prob = c(1, 0.5, 0.5),
                     tag_type = rep("c", 3))

  expect_error(suppressMessages(check_tags(tags, grid, verbose = FALSE)),
               "No tags passed the checks")
})


## ------------------------------------------------------------------- recovery

test_that("simulated candidate sets build and obey the mixture bound", {

  skip_on_cran()

  set.seed(42)

  grid <- create_grid(cellsize = 0.2, verbose = FALSE)
  sim <- suppressWarnings(suppressMessages(
    sim_data(grid = grid, n_ctags = 150, n_dtags = 0,
             n_candidates = 3, verbose = FALSE)))

  tags <- sim$tags
  expect_true("prob" %in% colnames(tags))
  expect_true(any(tags$event == 2 & duplicated(paste(tags$id, tags$event))))

  fit_ll <- function(tg) {
    dat <- sim$dat
    dat$tags <- tg
    conf <- sim$conf
    par <- sim$par
    map <- sim$map
    obj <- suppressWarnings(suppressMessages(
      admove(dat, conf, par, map, run = FALSE, verbose = FALSE)))$obj
    obj
  }

  ## picking one candidate at random throws away the weights
  w_kept <- numeric(0)
  picked <- do.call(rbind, lapply(split(tags, tags$id), function(z) {
    k <- which(z$event == max(z$event))
    j <- if (length(k) == 1L) k else sample(k, 1)
    w_kept <<- c(w_kept, z$prob[j])
    z <- z[c(setdiff(seq_len(nrow(z)), k), j), , drop = FALSE]
    z$event <- seq_len(nrow(z))
    z$prob <- 1
    z
  }))
  picked <- .add_class(picked, "admove_tags")
  sref(picked) <- sref(tags)
  tref(picked) <- tref(tags)

  o_mix <- fit_ll(tags)
  o_pick <- fit_ll(picked)

  ## both must build and evaluate to finite values with finite gradients
  expect_true(is.finite(o_mix$fn(o_mix$par)))
  expect_true(all(is.finite(o_mix$gr(o_mix$par))))
  expect_true(is.finite(o_pick$fn(o_pick$par)))

  ## sum_k w_k L_k >= w_j L_j for the candidate j that was kept, so in negative
  ## log space the mixture can be worse than a single candidate by at most
  ## -log(w_j). (It is NOT bounded by the single candidate itself: a mixture
  ## that spreads weight away from a very good candidate scores below it.)
  expect_lte(as.numeric(o_mix$fn(o_mix$par)),
             as.numeric(o_pick$fn(o_pick$par)) - sum(log(w_kept)) + 1e-8)
})


test_that("use_release_events keeps merged recoveries sequential, not ambiguous", {

  ## Three physical tags sharing a release event are pooled into one tag. Their
  ## recoveries are SUCCESSIVE observations of the pooled tag; if they kept the
  ## event id they carried individually (all 2) the likelihood would read them
  ## as mutually exclusive candidate positions for one recovery.
  grid <- create_grid(cellsize = 0.25, verbose = FALSE)

  tg <- data.frame(id = c("a", "a", "b", "b", "c", "c"),
                   t = c(0, 0.4, 0, 0.5, 0, 0.6),
                   x = c(0.2, 0.5, 0.2, 0.6, 0.2, 0.7),
                   y = c(0.2, 0.4, 0.2, 0.5, 0.2, 0.3),
                   tag_type = rep("c", 6))

  tags <- suppressMessages(check_tags(tg, grid, verbose = FALSE))
  expect_equal(tags$event, c(1, 2, 1, 2, 1, 2))

  out <- use_release_events(tags, grid, seq(0, 1, by = 0.1))

  expect_equal(nrow(out), 4L)
  expect_equal(out$event, 1:4)
  expect_equal(out$prob, rep(1, 4))
})
