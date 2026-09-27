## Main functions ---------------------------------------------------------------------

##' Predicted and observed tag positions
##'
##' @description
##' Returns, for every tag observation, the position the model predicts and the
##' position that was observed, together with the prediction uncertainty and
##' standardised residuals. This is the table behind [plot_tag_pred()],
##' [plot_tag_resid()] and [summarise_tag_pred()].
##'
##' The predictions are the ones the likelihood itself evaluates, so they
##' include the estimated observation error, the per-tag-type update rules
##' (`conf$do_update`) and the boundary treatment (`conf$kf_boundary`). With
##' the Kalman filter they are reported by `nll()`; with the CTMC engine they
##' come from the likelihood's own forward pass, which also gives the full
##' probability distribution over the grid cells, drawn by [plot_tag_pred()].
##'
##' @param fit A fitted object of class `admove`, as returned by [admove()].
##' @param type Which prediction to return:
##'   \describe{
##'     \item{`"osa"`}{One step ahead (the default): each observation is
##'       predicted from the one before it, as in the likelihood. For
##'       mark-recapture tags, which are not updated, this is the prediction
##'       from release to recapture.}
##'     \item{`"forecast"`}{From release: the track is never updated on the
##'       observations, so the prediction for an archival tag is a genuine
##'       forecast over the whole time at liberty. Needs a second model build
##'       and is therefore slower.}
##'   }
##' @param i Optional tag indices (into the tags split by id) or tag ids to
##'   restrict the table to. A small subset is much faster for
##'   `type = "forecast"` with the Kalman filter, where only those tags are
##'   built, and always with the CTMC engine, where every tag is a forward pass.
##' @param seed CTMC engine only: seed of the randomisation in the quantile
##'   residuals of observations without observation error (see Details), set
##'   locally so the caller's random numbers are not affected. `NULL` uses the
##'   current random stream. Default `1`, so repeated calls give the same
##'   residuals; try a few seeds to check that conclusions do not depend on it.
##' @param verbose Logical; if `TRUE` (default), messages about progress are
##'   printed.
##'
##' @details
##' With the Kalman filter (`conf$engine = "kf"`) the prediction is Gaussian
##' with independent x and y, so `z_x` and `z_y` should behave like standard
##' normal draws and `d2 = z_x^2 + z_y^2` like a chi-squared variable with 2
##' degrees of freedom. `inside_50` and `inside_95` compare `d2` with its
##' quantiles, i.e. they say whether the observation falls inside the 50% and
##' 95% prediction ellipse.
##'
##' With the CTMC engine (`conf$engine = "ctmc"`) the prediction is a
##' probability for every grid cell: the distribution each observation is
##' scored against in the likelihood, before its update. `pred_x`, `pred_y`
##' are its mean over the cell centres, and `sd_x`, `sd_y` its spread,
##' including the position within a cell (a cell of width \eqn{h} adds
##' \eqn{h^2/12}, as the CTMC does not resolve positions within a cell) and
##' the observation error where the likelihood applies one; `res_x`, `res_y`
##' and `dist` are measured from that mean. These describe the prediction but
##' are not standardised residuals, since the distribution can be skewed or
##' multimodal (e.g. along a coast).
##'
##' `z_x` and `z_y` are instead quantile residuals of the whole distribution,
##' by the Rosenblatt transform: `z_x` from the distribution of the grid
##' column (x) of the observation, `z_y` from the distribution of its row (y)
##' within that column. Under the model they are independent standard normal
##' whatever the shape of the distribution, so `d2`, `inside_50` and
##' `inside_95` are exact checks. For an observation without observation error
##' the likelihood only uses its grid cell, a discrete outcome, and the
##' residual is randomised within the probability of that cell (randomised
##' quantile residuals, Dunn and Smyth 1996; as one-step-ahead residuals of
##' state-space models, Thygesen et al. 2017; see `seed`). With observation
##' error the position is continuous (uniform within a cell plus the normal
##' error, as in the likelihood) and no randomisation is needed.
##'
##' `logdens` is the likelihood term of the observation divided by the cell
##' area, a density comparable with the Kalman filter's. The distributions
##' themselves are kept in `attr(, "ctmc")` for [plot_tag_pred()]; subsetting
##' the rows drops them.
##'
##' @references
##' Dunn, P. K. and Smyth, G. K. (1996). Randomized quantile residuals.
##' *Journal of Computational and Graphical Statistics* 5: 236--244.
##'
##' Rosenblatt, M. (1952). Remarks on a multivariate transformation. *Annals
##' of Mathematical Statistics* 23: 470--472.
##'
##' Thygesen, U. H., Albertsen, C. M., Berg, C. W., Kristensen, K. and
##' Nielsen, A. (2017). Validation of ecological state space models using the
##' Laplace approximation. *Environmental and Ecological Statistics* 24:
##' 317--339.
##'
##' The first observation of each tag is its release: it is conditioned on, not
##' predicted, and therefore not in the table.
##'
##' @return
##' A `data.frame` with one row per predicted observation and the columns
##' \describe{
##'   \item{`id`, `tag_type`}{Tag id and type.}
##'   \item{`obs`}{Index of the observation within the tag.}
##'   \item{`t`, `date`}{Time of the observation, in model units and as a date.}
##'   \item{`horizon`}{Time since the previous observation (`"osa"`) or since
##'     release (`"forecast"`).}
##'   \item{`x`, `y`}{Observed position.}
##'   \item{`x_from`, `y_from`}{Position the prediction starts from: the previous
##'     observation (`"osa"`) or the release (`"forecast"`).}
##'   \item{`pred_x`, `pred_y`}{Predicted position.}
##'   \item{`sd_x`, `sd_y`}{Prediction standard deviations.}
##'   \item{`res_x`, `res_y`, `dist`}{Observed minus predicted, and the distance
##'     between the two positions.}
##'   \item{`z_x`, `z_y`, `d2`}{Standardised residuals and their squared sum.}
##'   \item{`inside_50`, `inside_95`}{Whether the observation falls inside the
##'     50% / 95% prediction ellipse.}
##'   \item{`logdens`}{Log predictive density of the observation.}
##'   \item{`type`}{The `type` argument, so tables can be combined.}
##' }
##'
##' @examples
##' \donttest{
##' fit <- admove(skjepo$sim)
##' pred <- tag_predictions(fit)
##' head(pred)
##' summarise_tag_pred(pred)
##' }
##'
##' @seealso [plot_tag_pred()], [plot_tag_resid()], [summarise_tag_pred()]
##'
##' @export
tag_predictions <- function(fit,
                            type = c("osa", "forecast"),
                            i = NULL,
                            seed = 1,
                            verbose = TRUE) {

  .check_class(fit, "admove")
  type <- match.arg(type)

  tags <- .split_tags(fit$dat$tags)
  ids <- .resolve_tag_ids(i, names(tags))

  if (identical(.get_engine_integer(.get_engine_name(fit$conf$engine)), 2L)) {
    return(.ctmc_tag_pred(fit, ids, type, seed))
  }

  if (identical(type, "osa")) {
    rep <- .pred_report(fit)
  } else {
    if (verbose) {
      message("Building the forecast model for ", length(ids), " tag(s).")
    }
    rep <- .forecast_report(fit, ids)
    tags <- tags[ids]
  }

  out <- .tag_pred_table(tags, rep, tref(fit), type)

  out <- out[out$id %in% ids, , drop = FALSE]
  rownames(out) <- NULL

  if (nrow(out) == 0L) {
    warning("No predicted positions available for the selected tag(s).",
            call. = FALSE)
  }

  out
}


##' Prediction skill of one or more fitted models
##'
##' @description
##' Summarises how well predicted tag positions match the observed ones, from
##' the tables returned by [tag_predictions()]. Several tables (or fitted
##' models) can be passed at once to compare models, for example a model with
##' habitat preference against one with diffusion only.
##'
##' @param ... One or more tables from [tag_predictions()], or fitted `admove`
##'   objects (for which `tag_predictions(type = type)` is called). Names are
##'   used to label the rows.
##' @param by Optional column to group by, typically `"tag_type"` (the default)
##'   or `NULL` for one row per model.
##' @param type Passed to [tag_predictions()] when fitted models are supplied.
##'
##' @details
##' `rmse` is the root mean squared distance between the predicted and the
##' observed position, in the spatial units of the model. `cover_50` and
##' `cover_95` are the shares of observations inside the 50% and 95% prediction
##' ellipses: a well calibrated model is close to 0.5 and 0.95. Lower `rmse` and
##' higher `logdens` mean better predictions, while `sd_z` above 1 means the
##' model is more confident than it should be.
##'
##' @return
##' A `data.frame` with one row per model (and group) and the columns `model`,
##' the grouping column, `n`, `rmse`, `bias_x`, `bias_y`, `cover_50`,
##' `cover_95`, `sd_z` and `logdens`.
##'
##' @examples
##' \donttest{
##' fit <- admove(skjepo$sim)
##' conf0 <- fit$conf
##' conf0$use_taxis <- FALSE
##' fit0 <- admove(fit$dat, conf0)
##' summarise_tag_pred(preference = fit, diffusion_only = fit0)
##' }
##'
##' @seealso [tag_predictions()]
##'
##' @export
summarise_tag_pred <- function(..., by = "tag_type", type = c("osa", "forecast")) {

  type <- match.arg(type)
  inputs <- list(...)

  if (length(inputs) == 0L) stop("No prediction tables supplied.", call. = FALSE)

  nms <- names(inputs)
  if (is.null(nms)) nms <- rep("", length(inputs))
  nms[!nzchar(nms)] <- paste0("model", seq_along(inputs))[!nzchar(nms)]

  preds <- lapply(inputs, function(x) {
    if (inherits(x, "admove")) tag_predictions(x, type = type, verbose = FALSE) else x
  })

  for (k in seq_along(preds)) {
    if (!is.data.frame(preds[[k]]) || !all(c("res_x", "d2") %in% names(preds[[k]]))) {
      stop("Input ", k, " is neither an 'admove' fit nor a table from ",
           "tag_predictions().", call. = FALSE)
    }
  }

  if (!is.null(by) && !all(vapply(preds, function(p) by %in% names(p), logical(1)))) {
    stop("Column '", by, "' is missing from at least one prediction table.",
         call. = FALSE)
  }

  out <- do.call(rbind, lapply(seq_along(preds), function(k) {
    p <- preds[[k]]
    grp <- if (is.null(by)) rep("all", nrow(p)) else as.character(p[[by]])
    do.call(rbind, lapply(sort(unique(grp)), function(g) {
      pg <- p[grp == g, , drop = FALSE]
      data.frame(model = nms[k],
                 group = g,
                 n = nrow(pg),
                 rmse = sqrt(mean(pg$res_x^2 + pg$res_y^2, na.rm = TRUE)),
                 bias_x = mean(pg$res_x, na.rm = TRUE),
                 bias_y = mean(pg$res_y, na.rm = TRUE),
                 cover_50 = mean(pg$inside_50, na.rm = TRUE),
                 cover_95 = mean(pg$inside_95, na.rm = TRUE),
                 sd_z = stats::sd(c(pg$z_x, pg$z_y), na.rm = TRUE),
                 logdens = mean(pg$logdens, na.rm = TRUE),
                 stringsAsFactors = FALSE)
    }))
  }))

  names(out)[names(out) == "group"] <- if (is.null(by)) "group" else by
  rownames(out) <- NULL

  out
}


## Internal functions ---------------------------------------------------------------

## Tags as a list of data frames, in the same order as nll() sees them.
.split_tags <- function(tags) {
  if (inherits(tags, "list")) return(tags)
  split(as.data.frame(tags), tags$id)
}


## Tag ids selected by index or by id.
.resolve_tag_ids <- function(i, ids) {

  if (is.null(i)) return(ids)

  if (is.character(i)) {
    unknown <- setdiff(i, ids)
    if (length(unknown) > 0L) {
      stop("Unknown tag id(s): ", .format_ids(unknown), ".", call. = FALSE)
    }
    return(i)
  }

  i <- as.integer(i)
  if (any(is.na(i) | i < 1L | i > length(ids))) {
    stop("Tag index 'i' contains values outside [1, ", length(ids), "].",
         call. = FALSE)
  }

  ids[i]
}


## The reported predictions of a fit, from the stored report when there is one
## and otherwise by evaluating obj$report() at the estimates.
.pred_report <- function(fit) {

  if (!is.null(fit$rep$pred_set)) return(fit$rep)

  if (is.null(fit$obj) || !is.function(fit$obj$report)) {
    stop("This fit carries neither a report nor a usable model object. ",
         "Add the report with fit <- add_report(fit).", call. = FALSE)
  }

  pp <- fit$obj$env$last.par.best
  rep <- tryCatch(if (is.null(pp)) fit$obj$report() else fit$obj$report(pp),
                  error = function(e) NULL)

  if (is.null(rep$pred_set)) {
    stop("Could not evaluate the model report. Objects read back from disk ",
         "lose the pointers of their RTMB model, so refit or add the report ",
         "before saving: fit <- add_report(fit).", call. = FALSE)
  }

  rep
}


## Predictions with every update switched off, i.e. forecasts from release.
## Built for the selected tags only, because the model build dominates the cost.
.forecast_report <- function(fit, ids) {

  dat <- fit$dat
  tags <- dat$tags

  keep <- as.character(tags$id) %in% ids
  if (!any(keep)) stop("No tag rows left for the forecast.", call. = FALSE)
  dat$tags <- .subset_tags(tags, keep)

  conf <- fit$conf
  conf$do_update <- rep(FALSE, length(conf$do_update))

  obj <- suppressWarnings(suppressMessages(
    admove(dat, conf, fit$pl, fit$map, run = FALSE, verbose = FALSE)
  ))$obj

  obj$report(obj$par)
}


## Subset tag rows, keeping the class and the spatial/temporal reference.
.subset_tags <- function(tags, keep) {
  out <- as.data.frame(tags)[keep, , drop = FALSE]
  out <- .add_class(out, "admove_tags")
  out <- add_sref(out, sref(tags), verbose = FALSE)
  add_tref(out, tref(tags), verbose = FALSE)
}


## Assemble the prediction table from the tag list and a model report.
.tag_pred_table <- function(tags, rep, tr, type) {

  nobs <- vapply(tags, nrow, integer(1))
  offset <- c(0L, cumsum(nobs))

  out <- do.call(rbind, lapply(seq_along(tags), function(k) {

    tag <- tags[[k]]
    ind <- offset[k] + seq_len(nrow(tag))
    set <- rep$pred_set[ind] == 1

    if (!any(set)) return(NULL)

    t_obs <- as.numeric(tag$t)
    ## reference time of the prediction: the previous observation one step
    ## ahead, the release for a forecast from release
    x_obs <- as.numeric(tag$x)
    y_obs <- as.numeric(tag$y)
    if (identical(type, "forecast")) {
      t_from <- rep(t_obs[1], length(t_obs))
      x_from <- rep(x_obs[1], length(x_obs))
      y_from <- rep(y_obs[1], length(y_obs))
    } else {
      t_from <- c(t_obs[1], t_obs[-length(t_obs)])
      x_from <- c(x_obs[1], x_obs[-length(x_obs)])
      y_from <- c(y_obs[1], y_obs[-length(y_obs)])
    }

    data.frame(id = as.character(tag$id[1]),
               tag_type = as.character(tag$tag_type)[set],
               obs = which(set),
               t = t_obs[set],
               horizon = (t_obs - t_from)[set],
               x = x_obs[set],
               y = y_obs[set],
               x_from = x_from[set],
               y_from = y_from[set],
               pred_x = rep$pred_x[ind][set],
               pred_y = rep$pred_y[ind][set],
               sd_x = sqrt(rep$pred_var_x[ind][set]),
               sd_y = sqrt(rep$pred_var_y[ind][set]),
               stringsAsFactors = FALSE)
  }))

  if (is.null(out)) {
    out <- data.frame(id = character(0), tag_type = character(0),
                      obs = integer(0), t = numeric(0), horizon = numeric(0),
                      x = numeric(0), y = numeric(0),
                      x_from = numeric(0), y_from = numeric(0),
                      pred_x = numeric(0), pred_y = numeric(0),
                      sd_x = numeric(0), sd_y = numeric(0),
                      stringsAsFactors = FALSE)
  }

  out$res_x <- out$x - out$pred_x
  out$res_y <- out$y - out$pred_y
  out$dist <- sqrt(out$res_x^2 + out$res_y^2)
  out$z_x <- out$res_x / out$sd_x
  out$z_y <- out$res_y / out$sd_y
  ## per-row values of a flat report vector, in the order of the table rows
  set_rows <- function(v) unlist(lapply(seq_along(tags), function(k) {
    ind <- offset[k] + seq_len(nrow(tags[[k]]))
    v[ind][rep$pred_set[ind] == 1]
  }))
  ## CTMC: quantile residuals of the cell distribution instead
  if (!is.null(rep$pred_zx) && nrow(out) > 0L) {
    out$z_x <- set_rows(rep$pred_zx)
    out$z_y <- set_rows(rep$pred_zy)
  }
  out$d2 <- out$z_x^2 + out$z_y^2
  out$inside_50 <- out$d2 <= stats::qchisq(0.5, df = 2)
  out$inside_95 <- out$d2 <= stats::qchisq(0.95, df = 2)
  out$logdens <- stats::dnorm(out$res_x, 0, out$sd_x, log = TRUE) +
    stats::dnorm(out$res_y, 0, out$sd_y, log = TRUE)
  if (!is.null(rep$pred_logdens) && nrow(out) > 0L) {
    out$logdens <- set_rows(rep$pred_logdens)
  }
  out$type <- type

  date <- tryCatch(time_2_date(out$t, tref = tr), error = function(e) NULL)
  out$date <- if (is.null(date)) as.POSIXct(rep(NA_real_, nrow(out))) else date

  out[, c("id", "tag_type", "obs", "t", "date", "horizon",
          "x", "y", "x_from", "y_from", "pred_x", "pred_y", "sd_x", "sd_y",
          "res_x", "res_y", "dist", "z_x", "z_y", "d2",
          "inside_50", "inside_95", "logdens", "type")]
}


## CTMC counterpart of the reported KF predictions: the cell distribution each
## observation is scored against, recorded off the likelihood's own pass
## (.ctmc_tag_pass(record = )), summarised into the columns of the KF table.
## "forecast" switches the updates off, as .forecast_report() does for the KF.
## The distributions go into attr(, "ctmc") for plot_tag_pred().
## See dev/code_notes.org, "Tag location distributions".
.ctmc_tag_pred <- function(fit, ids, type, seed = 1) {

  ## the quantile residuals of exact positions are randomised: a local seed
  ## makes them reproducible without touching the caller's random stream
  if (!is.null(seed)) {
    return(.with_seed(seed, .ctmc_tag_pred(fit, ids, type, seed = NULL)))
  }

  dat <- fit$dat
  grid <- dat$grid
  tmb_all <- .nll_data(dat, fit$conf)
  if (identical(type, "forecast")) tmb_all$do_update[] <- FALSE
  tags <- tmb_all$tags

  cc <- .ctmc_ctx_from_fit(fit, tmb_all)
  ctx <- cc$ctx
  lat <- .ctmc_lattice(tmb_all)

  xc <- grid$xygrid[, 1]
  yc <- grid$xygrid[, 2]
  cs <- grid$cellsize
  nobs <- vapply(tags, nrow, integer(1))
  nall <- sum(nobs)
  offset <- c(0L, cumsum(nobs))
  rep <- list(pred_x = numeric(nall), pred_y = numeric(nall),
              pred_var_x = numeric(nall), pred_var_y = numeric(nall),
              pred_logdens = numeric(nall), pred_set = numeric(nall),
              pred_zx = numeric(nall), pred_zy = numeric(nall))
  dist <- list()

  for (idx in match(ids, names(tags))) {

    tag <- tags[[idx]]
    rec <- new.env(parent = emptyenv())
    .ctmc_tag_pass(ctx, idx, rep(0, length(tags)), cc$gen,
                   lat$time_mode, lat$breaks, record = rec)
    if (is.null(rec$dist)) next

    ev <- .tag_events(tag)
    last_ev <- max(ev)
    tt <- as.integer(tag$tag_type)
    p0 <- rep(0, ctx$nc)
    p0[tag$ic[1]] <- 1
    prob <- matrix(p0, ctx$nc, nrow(tag))

    for (j in seq_len(nrow(tag))[-1]) {
      p <- rec$dist[[j]]
      if (is.null(p)) next
      prob[, j] <- p
      if (is.na(tag$ic[j]) || tag$use[j] == 0) next

      mx <- sum(p * xc)
      my <- sum(p * yc)
      ## spread over the cells, the position within a cell (uniform: the CTMC
      ## does not resolve it), and the observation error where the likelihood
      ## applies one to this row
      ## same condition as in .ctmc_obs_dist()
      ovt <- tmb_all$obs_var_type[tt[j]]
      obs_err <- (ovt == 1 && ev[j] != last_ev) || ovt == 2 || ovt == 3
      this_dist <- .ctmc_obs_dist(ctx, tag, j, tt[j], ev, last_ev)
      so <- if (!obs_err) c(0, 0) else if (ovt == 3) {
        c(tag$sdx[j], tag$sdy[j])
      } else ctx$sdO[, tt[j]]

      k <- offset[idx] + j
      rep$pred_x[k] <- mx
      rep$pred_y[k] <- my
      rep$pred_var_x[k] <- sum(p * (xc - mx)^2) + cs[1]^2 / 12 + so[1]^2
      rep$pred_var_y[k] <- sum(p * (yc - my)^2) + cs[2]^2 / 12 + so[2]^2
      ## the likelihood term of the row, as a density per unit area so that it
      ## compares with the KF's Gaussian density
      rep$pred_logdens[k] <- log(sum(p * this_dist)) - log(cs[1] * cs[2])
      z <- .ctmc_quantile_resid(p, grid, tag$ic[j], tag$x[j], tag$y[j], so)
      rep$pred_zx[k] <- z[1]
      rep$pred_zy[k] <- z[2]
      rep$pred_set[k] <- 1
    }

    dist[[names(tags)[idx]]] <- list(prob = prob, tag = tag)
  }

  out <- .tag_pred_table(tags, rep, tref(fit), type)
  out <- out[out$id %in% ids, , drop = FALSE]
  rownames(out) <- NULL

  if (nrow(out) == 0L) {
    warning("No predicted positions available for the selected tag(s).",
            call. = FALSE)
  }

  attr(out, "ctmc") <- list(dist = dist, grid = grid,
                            xg = x_centers(grid), yg = y_centers(grid),
                            type = type)
  out
}


## Quantile residuals of one observation under a CTMC cell distribution `p`:
## the Rosenblatt transform (x from its marginal, then y given x), exactly
## N(0, 1) and independent under the model whatever the shape of p. See
## dev/code_notes.org, "Tag location distributions and release predictions".
##
## Without observation error (`sd` 0) the likelihood only uses the cell `ic`
## of the observation, so the residual is the randomised quantile residual of
## that discrete outcome (Dunn & Smyth 1996): uniform between the CDF values
## at the cell's edges, x over grid columns, y over the cells of the observed
## column. With observation error the likelihood treats the position as
## uniform within a cell plus N(0, sd^2) (the cell-integrated normal of
## .ctmc_obs_dist()), which is continuous: its CDF at the observed position,
## no randomisation.
.ctmc_quantile_resid <- function(p, grid, ic, x_obs, y_obs, sd = c(0, 0),
                                 eps = 1e-10) {

  ix <- grid$igrid$idx
  iy <- grid$igrid$idy

  if (all(sd == 0)) {

    i0 <- ix[ic]
    col <- ix == i0
    lo_x <- sum(p[ix < i0])
    u_x <- stats::runif(1, lo_x, lo_x + sum(p[col]))
    pcol <- sum(p[col])
    u_y <- if (pcol > 0) {
      lo_y <- sum(p[col & iy < iy[ic]]) / pcol
      stats::runif(1, lo_y, lo_y + p[ic] / pcol)
    } else NA_real_

  } else {

    ## Uniform(a, b) + N(0, s^2): CDF and density
    G <- function(v, a, b, s) {
      psi <- function(z) z * stats::pnorm(z) + stats::dnorm(z)
      s / (b - a) * (psi((v - a) / s) - psi((v - b) / s))
    }
    g <- function(v, a, b, s) {
      (stats::pnorm((v - a) / s) - stats::pnorm((v - b) / s)) / (b - a)
    }
    ax <- grid$xgr[ix]
    bx <- grid$xgr[ix + 1L]
    ay <- grid$ygr[iy]
    by <- grid$ygr[iy + 1L]
    u_x <- sum(p * G(x_obs, ax, bx, sd[1]))
    w <- p * g(x_obs, ax, bx, sd[1])
    u_y <- if (sum(w) > 0) sum(w * G(y_obs, ay, by, sd[2])) / sum(w) else NA_real_
  }

  stats::qnorm(pmin(pmax(c(u_x, u_y), eps), 1 - eps))
}
