## Main functions ---------------------------------------------------------------------

##' Predicted and observed tag positions
##'
##' @description
##' Returns, for every tag observation, the position the model predicts and the
##' position that was observed, together with the prediction uncertainty and
##' standardised residuals. This is the table behind [plot_tag_pred()],
##' [plot_tag_resid()] and [summarise_tag_pred()].
##'
##' The predictions are the ones the likelihood itself evaluates, reported by
##' `nll()` rather than recomputed, so they include the estimated observation
##' error, the per-tag-type update rules (`conf$do_update`) and the boundary
##' treatment (`conf$kf_boundary`).
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
##' @param i Optional tag indices (as in [add_tag_dist()]) or tag ids to restrict
##'   the table to. For `type = "forecast"` a small subset is much faster,
##'   because only those tags are built.
##' @param verbose Logical; if `TRUE` (default), messages about progress are
##'   printed.
##'
##' @details
##' Predicted positions are only available for the Kalman filter engine
##' (`conf$engine = 1`), which tracks a mean and a variance per tag. The CTMC
##' engine carries a distribution over grid cells instead; use [add_tag_dist()]
##' and [plot_tag_dist()] for it.
##'
##' The prediction is Gaussian with independent x and y, so `z_x` and `z_y`
##' should behave like standard normal draws and `d2 = z_x^2 + z_y^2` like a
##' chi-squared variable with 2 degrees of freedom. `inside_50` and `inside_95`
##' compare `d2` with its quantiles, i.e. they say whether the observation falls
##' inside the 50% and 95% prediction ellipse.
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
                            verbose = TRUE) {

  .check_class(fit, "admove")
  type <- match.arg(type)

  if (!identical(.get_engine_integer(fit$conf$engine), 1L)) {
    stop("Predicted tag positions are only available for the Kalman filter ",
         "engine (conf$engine = 1). The CTMC engine carries a distribution ",
         "over grid cells instead: see add_tag_dist() and plot_tag_dist().",
         call. = FALSE)
  }

  tags <- .split_tags(fit$dat$tags)
  ids <- .resolve_tag_ids(i, names(tags))

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
  out$d2 <- out$z_x^2 + out$z_y^2
  out$inside_50 <- out$d2 <= stats::qchisq(0.5, df = 2)
  out$inside_95 <- out$d2 <= stats::qchisq(0.95, df = 2)
  out$logdens <- stats::dnorm(out$res_x, 0, out$sd_x, log = TRUE) +
    stats::dnorm(out$res_y, 0, out$sd_y, log = TRUE)
  out$type <- type

  date <- tryCatch(time_2_date(out$t, tref = tr), error = function(e) NULL)
  out$date <- if (is.null(date)) as.POSIXct(rep(NA_real_, nrow(out))) else date

  out[, c("id", "tag_type", "obs", "t", "date", "horizon",
          "x", "y", "x_from", "y_from", "pred_x", "pred_y", "sd_x", "sd_y",
          "res_x", "res_y", "dist", "z_x", "z_y", "d2",
          "inside_50", "inside_95", "logdens", "type")]
}
