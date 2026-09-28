##' Predict where released animals are over time
##'
##' @description
##' `release_predictions()` releases animals at one or more locations and times
##' and predicts, from a fitted model, the probability that an animal is in
##' each grid cell at given times after its release. Multiplied by the number
##' of animals released, `n`, this is the expected number of animals in each
##' cell. [plot_release_pred()] draws the result.
##'
##' @param fit A fitted object of class `admove`, as returned by [admove()].
##' @param x0,y0 Release locations, in the coordinates of the model (see
##'   [units_space()]). Several releases can be given as vectors; `x0`, `y0`,
##'   `t0` and `n` are recycled to a common length.
##' @param t0 Release times, in model time (as the tag times; see
##'   [date_2_time()] to convert dates).
##' @param times Times after release at which to predict, in model time units
##'   (e.g. months), the same for every release. Default: `1:12`.
##' @param n Number of animals released at each release. Only scales the
##'   predicted distribution into expected numbers of animals. Default: `1`.
##'
##' @details
##' Animals released together at the same place and time move independently
##' under the model, so all of them have the same probability distribution over
##' the grid cells at any later time: one prediction serves any number of
##' animals, and `n` only scales it. The number of animals in a cell at one
##' time is binomial with `n` and the cell probability. The prediction says
##' nothing about the path of individual animals (e.g. how many animals that
##' are in one area in month 3 reach another area in month 6); that needs
##' individually simulated tracks.
##'
##' The prediction uses the continuous-time Markov chain (CTMC) of the model on
##' the model grid (`fit$dat$grid`) with the estimated taxis, diffusion and
##' advection, stepped exactly as the CTMC likelihood steps a tag (see
##' [tag_predictions()]): the same time steps (at most `fit$dat$min_dt`, aligned
##' to the covariate slices), the same generator and the same matrix
##' exponential (`conf$ctmc_method`). A release at a tag's release location
##' and time therefore reproduces that tag's forecast. It also works for a
##' model fitted with the Kalman filter, which estimates the same movement
##' model in continuous space; the CTMC is then its discretisation on the grid.
##'
##' Movement is only defined where the covariates are: releases and
##' prediction times outside the covariate time range give a warning, as in
##' [add_predictions()].
##'
##' @return An object of class `admove_release_pred`: a list with one element
##'   per release (`x0`, `y0`, `t0`, `n`, the release cell `ic`, the
##'   prediction times after release `times` and in model time `t`, and `prob`,
##'   a matrix of cell probabilities with one column per prediction time), and
##'   the grid and reference systems needed to plot it.
##'
##' @seealso [plot_release_pred()], [tag_predictions()], [sim_tags()] for
##'   individually simulated tracks
##'
##' @examples
##' \dontrun{
##' rp <- release_predictions(fit, x0 = 500, y0 = -800, t0 = 48, n = 1000)
##' plot_release_pred(rp, type = "count", plot_land = TRUE)
##' }
##'
##' @export
release_predictions <- function(fit, x0, y0, t0, times = 1:12, n = 1) {

  .check_class(fit, "admove")

  dat <- fit$dat
  grid <- dat$grid

  nrel <- max(length(x0), length(y0), length(t0), length(n))
  if (min(length(x0), length(y0), length(t0), length(n)) == 0L) {
    stop("'x0', 'y0', 't0' and 'n' must not be empty.", call. = FALSE)
  }
  rel <- data.frame(x0 = rep_len(as.numeric(x0), nrel),
                    y0 = rep_len(as.numeric(y0), nrel),
                    t0 = rep_len(as.numeric(t0), nrel),
                    n = rep_len(as.numeric(n), nrel))
  if (any(!is.finite(as.matrix(rel)))) {
    stop("'x0', 'y0', 't0' and 'n' must be finite numbers.", call. = FALSE)
  }
  if (any(rel$n <= 0)) stop("'n' must be positive.", call. = FALSE)

  times <- sort(unique(as.numeric(times)))
  if (length(times) == 0L || any(!is.finite(times)) || any(times <= 0)) {
    stop("'times' must be positive times after release.", call. = FALSE)
  }

  if (is.null(grid)) {
    stop("The fit has no grid (fit$dat$grid), which the prediction needs.",
         call. = FALSE)
  }
  rel$ic <- grid$celltable[cbind(cut(rel$x0, grid$xgr), cut(rel$y0, grid$ygr))]
  bad <- which(is.na(rel$ic))
  if (length(bad) > 0L) {
    stop("Release location(s) ", .format_ids(bad), " are not in a grid cell ",
         "of the model grid.", call. = FALSE)
  }

  ## same warnings as for predictions outside the covariate time range
  dcheck <- dat
  dcheck$pred$time <- sort(unique(c(rel$t0, outer(times, rel$t0, "+"))))
  .check_pred_time_coverage(dcheck, fit$conf)

  ## the CTMC context of the likelihood, and its lattice; the slice
  ## boundaries are taken over the release windows rather than the tags
  tmb_all <- .nll_data(dat, fit$conf)
  cc <- .ctmc_ctx_from_fit(fit, tmb_all)
  ctx <- cc$ctx
  gen <- cc$gen
  dbr <- tmb_all
  dbr$tags <- data.frame(t = range(dcheck$pred$time))
  lat <- .ctmc_lattice(dbr)

  res <- vector("list", nrel)
  for (r in seq_len(nrel)) {

    t_obs <- c(rel$t0[r], rel$t0[r] + times)
    ## the time step rule of .ctmc_tag_pass(), so that a release at a tag's
    ## release, predicted at the tag's times, reproduces its forecast
    dt_min <- min(tmb_all$min_dt, stats::median(diff(t_obs)), na.rm = TRUE)
    out <- .build_time_breaks(t_obs,
                              mode = lat$time_mode,
                              dt_min = dt_min,
                              dt = tmb_all$dt,
                              eps = 0.1,
                              breaks = lat$breaks)

    prob <- matrix(0, ctx$nc, length(times))
    last_dist <- rep(0, ctx$nc)
    last_dist[rel$ic[r]] <- 1
    for (k in 2:out$nts) {
      dt <- out$dts[k - 1L]
      last_dist <- .ctmc_step(ctx, gen$Q(out$ts[k - 1L], dt), dt, last_dist)
      prob[, out$observed == k] <- last_dist
    }

    res[[r]] <- list(x0 = rel$x0[r],
                     y0 = rel$y0[r],
                     t0 = rel$t0[r],
                     n = rel$n[r],
                     ic = rel$ic[r],
                     times = times,
                     t = rel$t0[r] + times,
                     prob = prob)
  }

  structure(list(releases = res,
                 grid = grid,
                 xg = x_centers(grid),
                 yg = y_centers(grid),
                 sref = sref(dat),
                 tref = tref(dat),
                 units_space = tryCatch(units_space(dat), error = function(e) NULL),
                 units_time = tryCatch(units_time(dat), error = function(e) NULL)),
            class = "admove_release_pred")
}


##' @rdname release_predictions
##' @param x An object of class `admove_release_pred`.
##' @param ... Unused.
##' @export
print.admove_release_pred <- function(x, ...) {
  rel <- x$releases
  cat("Predicted distribution of released animals\n")
  cat("  releases:", length(rel), "\n")
  for (r in seq_along(rel)) {
    z <- rel[[r]]
    cat(sprintf("  %d: x0 = %s, y0 = %s, t0 = %s, n = %s\n", r,
                .fmt_num(z$x0), .fmt_num(z$y0), .fmt_num(z$t0), .fmt_num(z$n)))
  }
  tt <- rel[[1L]]$times
  un <- x$units_time
  cat("  times after release:", paste(.fmt_num(range(tt)), collapse = " to "),
      if (length(un) == 1L && !is.na(un)) un, paste0("(", length(tt), " times)"),
      "\n")
  invisible(x)
}


##' Plot the predicted distribution of released animals
##'
##' @description
##' Draws the predictions of [release_predictions()]: the probability that a
##' released animal is in each grid cell, or the expected number of animals,
##' at each prediction time. With one release, the times are spread over a
##' grid of panels; with several releases, each release is a row (its number
##' leads the panel titles) and each time a column.
##'
##' @param x An object of class `admove_release_pred`, as returned by
##'   [release_predictions()].
##' @param select Prediction times to show, as indices into the times of `x`.
##'   Default: `NULL`, at most `n_time_steps` times evenly spread.
##' @param select_release Releases to show, as indices. Default: `NULL`, all.
##' @param n_time_steps Maximum number of times shown when `select` is `NULL`.
##'   Default: `6`.
##' @param type `"prob"` (default) for the probability per cell, `"count"` for
##'   the expected number of animals per cell (the probability times `n`).
##' @param min_prob Smallest cell probability drawn; smaller ones are left
##'   blank. The colours are on a log10 scale, shared by all panels. Default:
##'   `1e-4`.
##' @param col Colour palette. Default: `NULL`, the light teal palette for
##'   estimated quantities shared with [plot_taxis()] and the other plots of
##'   estimates.
##' @param legend Logical; if `TRUE` (default), one colour bar for all panels is
##'   drawn in the right margin.
##' @param plot_land Logical; if `TRUE`, land is added: filled underneath the
##'   prediction, coastline on top. Default: `FALSE`.
##' @param land_col,land_border Fill and coastline colour of the land.
##' @param main Main title, drawn above the panels. Default: `NULL`, none.
##' @param xlab,ylab Axis labels. If `NULL` (default), `"x"` and `"y"` with the
##'   spatial units in brackets.
##' @param asp Target aspect ratio passed to [grDevices::n2mfrow()] for the
##'   panel layout of a single release. Default: `1`.
##'
##' @details
##' Each panel is titled with its date (at the resolution of the time units,
##' e.g. `"Mar 2007"` for monthly units) or, without a time origin, with its
##' model time; the time after release is added in brackets. The release
##' location is marked with a cross.
##'
##' @return Invisibly `NULL`.
##'
##' @seealso [release_predictions()]
##'
##' @export
plot_release_pred <- function(x,
                              select = NULL,
                              select_release = NULL,
                              n_time_steps = 6L,
                              type = c("prob", "count"),
                              min_prob = 1e-4,
                              col = NULL,
                              legend = TRUE,
                              plot_land = FALSE,
                              land_col = grey(0.85),
                              land_border = grey(0.3),
                              main = NULL,
                              xlab = NULL,
                              ylab = NULL,
                              asp = 1) {

  .check_class(x, "admove_release_pred")
  type <- match.arg(type)
  if (is.null(col)) col <- .est_col()
  if (!is.numeric(min_prob) || length(min_prob) != 1L || !(min_prob > 0) ||
        min_prob >= 1)
    stop("'min_prob' must be a single number in (0, 1).", call. = FALSE)

  un <- x$units_space
  if (is.null(xlab)) xlab <- .axis_lab("x", un)
  if (is.null(ylab)) ylab <- .axis_lab("y", un)

  rel <- x$releases
  if (is.null(select_release)) select_release <- seq_along(rel)
  if (any(!select_release %in% seq_along(rel))) {
    stop("'select_release' must be indices of the releases (1 to ",
         length(rel), ").", call. = FALSE)
  }
  rel <- rel[select_release]
  nt <- length(rel[[1L]]$times)
  if (is.null(select)) {
    select <- if (nt <= n_time_steps) seq_len(nt) else
      unique(round(seq(1L, nt, length.out = n_time_steps)))
  }
  if (any(!select %in% seq_len(nt))) {
    stop("'select' must be indices of the prediction times (1 to ", nt, ").",
         call. = FALSE)
  }

  ## values drawn and their log10 colour scale, shared by all panels
  val <- function(z, k) {
    p <- z$prob[, k]
    if (type == "count") p * z$n else p
  }
  low <- min(vapply(rel, function(z) if (type == "count") min_prob * z$n else
    min_prob, numeric(1)))
  high <- max(vapply(rel, function(z) max(val(z, select)), numeric(1)))
  zlim <- if (high > low) log10(c(low, high)) else NULL
  use_legend <- legend && !is.null(zlim)

  grid <- x$grid
  on_grid <- function(p) {
    ct <- grid$celltable
    ok <- !is.na(ct)
    ct[ok] <- p[ct[ok]]
    ct
  }
  xrange <- grid$xrange + c(-0.05, 0.05) * diff(grid$xrange)
  yrange <- grid$yrange + c(-0.05, 0.05) * diff(grid$yrange)
  sr <- if (plot_land) x$sref

  opar <- par(no.readonly = TRUE)
  on.exit(suppressWarnings(graphics::par(opar)))

  nrel <- length(rel)
  ns <- length(select)
  mfrow <- if (nrel == 1L) n2mfrow(ns, asp = asp) else c(nrel, ns)
  has_main <- !is.null(main) && nzchar(main)
  par(mfrow = mfrow, mar = c(0.2, 0.2, 1.4, 0.2),
      oma = c(3.5, 4, if (has_main) 2 else 0.5, if (use_legend) 5 else 1),
      mgp = c(2, 0.5, 0), tcl = -0.3)

  n_panel <- nrel * ns
  nc <- mfrow[2L]
  p <- 0L
  for (r in seq_len(nrel)) {
    z <- rel[[r]]
    lab_t <- .time_labels(z$t, x$tref)
    lab_rel <- paste0("+", .fmt_num(z$times))
    for (k in select) {
      p <- p + 1L
      plot(NA, NA, xlim = xrange, ylim = yrange, asp = 1,
           xaxt = "n", yaxt = "n", xlab = "", ylab = "")
      col_p <- ((p - 1L) %% nc) + 1L
      if (p + nc > n_panel) axis(1)
      if (col_p == 1L) axis(2)
      .prob_image(x$xg, x$yg, on_grid(val(z, k)), zlim, col, sr,
                  land_col, land_border)
      points(z$x0, z$y0, pch = 4, lwd = 2, cex = 1.2, col = "grey10")
      title(main = paste0(if (nrel > 1L) paste0(select_release[r], ": "),
                          lab_t[k], " (", lab_rel[k], ")"),
            line = 0.3, font.main = 1, cex.main = 0.9)
      box(lwd = 1.5)
    }
  }
  ## pad the grid of a single release
  if (nrel == 1L) for (q in seq_len(prod(mfrow) - n_panel)) plot.new()

  mtext(xlab, 1, 2, outer = TRUE)
  mtext(ylab, 2, 2.5, outer = TRUE)
  if (has_main) mtext(main, 3, 0.5, outer = TRUE, font = 2)

  if (use_legend) {
    .log_color_bar(col, zlim, if (type == "count") "animals" else "probability")
  }

  invisible(NULL)
}
