
## Main functions ---------------------------------------------------------------------

##' Set up input data for an admove model
##'
##' @description
##' Combine the main input components for an `admove` analysis into a single
##' object of class `admove_data`. The function checks and harmonises spatial and
##' temporal references, prepares covariates and tags, defines spline knots,
##' and constructs default prediction grids and time points used during fitting
##' and plotting.
##'
##' @param grid Optional grid object, typically of class `admove_grid`, as
##'   returned by [create_grid()]. A grid is required for some model engines and
##'   for spatial prediction and plotting.
##' @param cov Optional covariate object or list of covariates. Covariates are
##'   typically prepared with [prep_cov()]. If a single covariate is supplied,
##'   it is coerced internally to a list.
##' @param tags Optional tag data, typically as returned by one or more of
##'   [prep_ctags()], [prep_dtags()], or [prep_stags()].
##' @param trange Optional numeric vector of length two giving the model time
##'   range. If `NULL`, the time range is inferred from available tags and
##'   covariates.
##' @param knots_tax Optional matrix of spline knots for the taxis preference
##'   functions. If `NULL`, default knots are chosen from covariate quantiles.
##' @param knots_dif Optional matrix of spline knots for the diffusion
##'   preference functions. If `NULL`, default knots are chosen from covariate
##'   quantiles.
##' @param sref Optional spatial reference to use as the target spatial
##'   reference for all inputs. If supplied, it should be coercible to an
##'   `admove_sref` object.
##' @param tref Optional time reference to use as the target time reference for
##'   all inputs. If supplied, it should be coercible to an `admove_tref`
##'   object.
##' @param transform_sref Logical; if `TRUE`, spatial components are
##'   transformed or rescaled to a common spatial reference where possible. If
##'   `FALSE` (default), all spatial references must already match.
##' @param shift_tref Logical; if `TRUE`, temporal components are shifted or
##'   harmonised to a common time reference where possible. If `FALSE`
##'   (default), all time references must already match.
##' @param verbose Logical; if `TRUE`, informative messages are printed during
##'   processing. Default is `TRUE`.
##'
##' @details
##' The function first determines a common spatial reference (`sref`) and time
##' reference (`tref`) either from the user-supplied targets or from the input
##' objects. If `transform_sref = FALSE` or `shift_tref = FALSE`, all inputs
##' must already be compatible. Otherwise, inputs are harmonised to the chosen
##' target references.
##'
##' If `trange` is not supplied, it is inferred from the union of tag and
##' covariate time ranges. If no valid time range can be determined, a default
##' range of `c(0, 1)` is used.
##'
##' If spline knots are not supplied, default knots are generated from the
##' marginal covariate distributions. Taxis knots use three quantiles per
##' covariate, while diffusion knots use one quantile per covariate.
##'
##' The returned object also contains default prediction components in
##' `dat$pred`, including:
##' \describe{
##'   \item{`pred$time`}{A sequence of 10 prediction time points over `trange`.}
##'   \item{`pred$cov`}{Prediction ranges for each covariate.}
##'   \item{`pred$grid`}{The prediction grid, if a grid was supplied.}
##' }
##'
##' Additional defaults used later during fitting are also stored in the output,
##' including `eps`, `var_init_kf`, `log2steps`, `min_dt`, `dt`, and `p_init`.
##'
##' @return
##' An object of class `admove_data`, containing the processed grid, covariates,
##' tags, spline knots, prediction settings, and associated spatial and temporal
##' reference information.
##'
##' @examples
##' ctags <- prep_tags(
##'   skjepo$ctags,
##'   tag_type = "c",
##'   names = c(
##'     t0 = "date_time", t1 = "date_caught",
##'     x0 = "rel_lon",   x1 = "recap_lon",
##'     y0 = "rel_lat",   y1 = "recap_lat"
##'   ),
##'   date_origin = "1899-12-30")
##'
##' ## prepare data-storage tags
##' dtags <- prep_tags(
##'   skjepo$dtags,
##'   tag_type = "d",
##'   names = c(t = "time", x = "mptlon", y = "mptlat"),
##'   date_origin = "1899-12-30")
##'
##' grid <- create_grid(x = skjepo$grid, cellsize = c(10, 10))
##'
##' cov <- prep_cov(skjepo$cov)
##'
##' dat <- setup_data(
##'   grid = grid,
##'   cov = cov,
##'   tags = c(dtags, ctags),
##'   transform_sref = TRUE,
##'   shift_tref = TRUE
##' )
##'
##' @export
setup_data <- function(grid = NULL,
                       cov = NULL,
                       tags = NULL,
                       trange = NULL,
                       knots_tax = NULL,
                       knots_dif = NULL,
                       sref = NULL,
                       tref = NULL,
                       transform_sref = FALSE,
                       shift_tref = FALSE,
                       verbose = TRUE) {

  res <- list()

  if (!is.null(cov)) cov <- .make_cov_list(cov)

  ## choose master sref
  if (!is.null(sref)) {
    sref_target <- create_sref(sref$crs, sref$units, sref$crs_scale)
  } else sref_target <- NULL
  if (!is.null(sref_target)) {
    if (!inherits(sref_target, "admove_sref")) {
      stop("'sref_target' must be an admove_sref.")
    }
    master_sref <- sref_target
  } else if (!is.null(grid)) {
    master_sref <- sref(grid)
  } else if (!is.null(cov)) {
    master_sref <- sref(cov[[1]])
  } else if (!is.null(tags)) {
    master_sref <- sref(tags)
  } else master_sref <- NULL

  if (transform_sref) {

    ## harmonise everything to master (units only; CRS must match)
    if (!is.null(grid)) grid <- add_sref(grid,  master_sref, verbose, transform_sref)
    if (!is.null(tags)) tags <- add_sref(tags, master_sref, verbose, transform_sref)
    if (!is.null(cov)){
      cov <- lapply(cov, add_sref, sref = master_sref, verbose = verbose,
                    transform_crs = transform_sref)
      cov <- .add_class(cov, "admove_cov_list")
    }

  } else {

    ## strict check
    bad <- character(0)

    if (!is.null(grid) &&
          !sref_equal(sref(grid), master_sref))  bad <- c(bad, "grid")
    if (!is.null(tags) &&
          !sref_equal(sref(tags), master_sref)) bad <- c(bad, "tags")

    cov_bad <- NULL
    if (!is.null(cov)) cov_bad <- which(vapply(cov,
                                                   function(z)
                                                     !sref_equal(sref(z),
                                                                  master_sref), logical(1)))
    if (length(cov_bad)) bad <- c(bad, paste0("cov[[", cov_bad, "]]"))

    if (length(bad)) {
      if (is.null(sref_target)) {
        stop("Spatial reference mismatch for: ", paste(bad, collapse = ", "),
             ". Provide 'sref' AND set transform_sref=TRUE, or fix inputs.")
      } else {
        stop("Spatial reference mismatch for: ", paste(bad, collapse = ", "),
             ". Set transform_sref=TRUE, or fix inputs.")
      }
    }
  }

  ## choose master tref
  if (!is.null(tref)) {
    tref_target <- create_tref(tref$origin, tref$units, tref$period)
  } else tref_target <- NULL
  if (!is.null(tref_target)) {
    if (!inherits(tref_target, "admove_tref")) {
      stop("'tref_target' must be an admove_tref.")
    }
    master_tref <- tref_target
  } else if (!is.null(cov)) {
    master_tref <- tref(cov[[1]])
  } else if (!is.null(tags)) {
    master_tref <- tref(tags)
  } else master_tref <- NULL

  if (shift_tref) {

    if (!is.null(tags)) tags <- add_tref(tags, master_tref, verbose, shift_tref)
    if (!is.null(cov)){
      ## Warn when a covariate has no tref origin: add_tref will label it with
      ## the new origin but cannot shift the stored time values, so the time
      ## axis will remain in whatever system the raw array was in (e.g. absolute
      ## decimal years). Use date_decimal = TRUE / date_format in prep_cov() to
      ## give the covariate a proper origin before calling setup_data().
      no_origin <- vapply(cov, function(z) {
        tr <- try(tref(z), silent = TRUE)
        if (inherits(tr, "try-error") || is.null(tr)) return(TRUE)
        .is_na_scalar(tr$origin)
      }, logical(1))
      if (any(no_origin)) {
        warning(
          "shift_tref = TRUE but the following covariate(s) have no tref ",
          "origin (tref$origin is NA): cov[[",
          paste(which(no_origin), collapse = ", "),
          "]]. The time values in these covariates cannot be shifted and will ",
          "remain in their original time system. This will likely cause a ",
          "time-axis mismatch with the tags, making the covariate(s) ",
          "inaccessible during fitting (silent zero-gradient). Fix: call ",
          "prep_cov(..., date_decimal = TRUE) or supply date_format / ",
          "date_origin so the covariate gets a proper tref before setup_data().",
          call. = FALSE
        )
      }
      cov <- lapply(cov, add_tref, tref = master_tref, verbose = verbose,
                    shift_origin = shift_tref)
      cov <- .add_class(cov, "admove_cov_list")
    }


  } else {

    ## strict check
    bad <- character(0)

    if (!is.null(tags) &&
          !tref_equal(tref(tags), master_tref)) bad <- c(bad, "tags")

    cov_bad <- NULL
    if (!is.null(cov)) cov_bad <- which(vapply(cov,
                                               function(z)
                                                 !tref_equal(tref(z),
                                                              master_tref), logical(1)))
    if (length(cov_bad)) bad <- c(bad, paste0("cov[[", cov_bad, "]]"))


    if (length(bad)) {
      if (is.null(tref_target)) {
        stop("Time reference mismatch for: ", paste(bad, collapse = ", "),
             ". Provide 'tref' AND set shift_tref=TRUE, or fix inputs.")
      } else {
        stop("Time reference mismatch for: ", paste(bad, collapse = ", "),
             ". Set shift_tref=TRUE, or fix inputs.")
      }
    }
  }


  ## Time range -----------------------------------------
  if (!is.null(tags)) {
    dim_tags <- get_dim_tags(tags)
    trange <- c(trange, range(dim_tags$trange))
  }
  if (!is.null(cov)) {
    trange <- c(trange, range(.get_cov_trange(cov)))
  }
  trange <- range(trange)

  if (is.null(trange) || any(is.infinite(trange))) {
    trange <- c(0,1)
  } else if (trange[2] == trange[1]) {
    trange[2] <- trange[1] + 1
  }

  res$trange <- trange

  ## Space -------------------------------------------
  res$grid <- check_grid(grid)


  ## Covariates --------------------------------------
  res$cov <- check_cov(cov, verbose)

  if (!is.null(res$cov)) {

    ## space
    xyranges_cov <- .get_cov_xyrange(res$cov)
    res$xrange_cov <- xyranges_cov$xr
    res$yrange_cov <- xyranges_cov$yr

    ## Check if any cov NA where needed for interpol given provided grid
    if (!is.null(grid)) {
      xr <- xyranges_cov$xr
      yr <- xyranges_cov$yr
      ncov <- length(res$cov)
      err <- NULL
      for (i in seq_len(ncov)) {
        covi <- res$cov[[i]]
        for (j in seq_len(dim(covi)[3])) {
          liv <- RTMB::interpol2Dfun(covi[,,j],
                                     xlim = round(xr[i,],5),
                                     ylim = round(yr[i,],5),
                                     R = 1)
          tmp <- liv(round(grid$xygrid[,1],5),
                     round(grid$xygrid[,2],5))
          ind <- which(is.na(tmp))
          if (length(ind) > 0) {
            err <- c(err, ind)
          }
        }
      }
      err <- sort(unique(err))

      if (length(err) > 0) {

        message("For the following grid cell(s) the covariate(s) can not be calculated (lead to NA): ", paste(err, collapse = ", "), ". Removing grid cell(s) in order to avoid problems during fitting later. \n")

        ind <- match(err, res$grid$celltable)
        res$grid$celltable[ind] <- NA
        res$grid$celltable[!is.na(res$grid$celltable)] <-
          1:sum(!is.na(res$grid$celltable))
        res$grid$xygrid <- res$grid$xygrid[-err,]
        res$grid$igrid <- res$grid$igrid[-err,]

      }
    }

    ## times
    res$time_cov <- lapply(res$cov, function(x) as.numeric(dimnames(x)[[3]]))

  } else {
    res$xrange_cov <- NULL
    res$yrange_cov <- NULL
    res$time_cov <- NULL
  }

  ## Tags --------------------------------------------
  if (!is.null(tags)) res$tags <- check_tags(tags, res$grid)

  ## Remove tag positions on NA covariate cells -------
  ## check_tags() only guards against NA *grid* cells. A tag can still sit on
  ## (or next to) a covariate cell that is NA, where the bilinear interpolation
  ## used in the likelihood returns NaN and causes NaN objective/gradient
  ## evaluations during minimisation. Evaluate each covariate at the tag
  ## positions, matching how the likelihood accesses them (only time slices with
  ## t2index() > 0 are used), and drop the offending entries.
  if (!is.null(res$cov) && !is.null(res$tags) && nrow(res$tags) > 0) {
    xr <- res$xrange_cov
    yr <- res$yrange_cov
    bad <- rep(FALSE, nrow(res$tags))
    for (i in seq_along(res$cov)) {
      covi <- res$cov[[i]]
      it <- as.integer(t2index(res$tags$t, res$time_cov[[i]]))
      for (j in sort(unique(it[it > 0]))) {
        rows <- which(it == j)
        liv <- RTMB::interpol2Dfun(covi[,,j],
                                   xlim = round(xr[i,], 5),
                                   ylim = round(yr[i,], 5),
                                   R = 1)
        v <- liv(round(res$tags$x[rows], 5), round(res$tags$y[rows], 5))
        bad[rows] <- bad[rows] | is.na(v)
      }
    }
    if (any(bad)) {
      bad_ids <- unique(res$tags$id[bad])
      if (verbose) {
        message(sum(bad), " entr", if (sum(bad) == 1) "y" else "ies",
                " removed because the tag position falls where a covariate is NA",
                " (tag id", if (length(bad_ids) == 1) "" else "s", ": ",
                paste(bad_ids, collapse = ", "), ").")
      }
      res$tags <- res$tags[!bad, , drop = FALSE]
      ## Drop tags left with a single observation (need release + recovery),
      ## matching the recovered-tags filter in check_tags().
      keep_id <- names(which(table(res$tags$id) > 1))
      res$tags <- res$tags[res$tags$id %in% keep_id, , drop = FALSE]
      if (nrow(res$tags) == 0) {
        stop("No tags remain after removing positions on NA covariate cells.")
      }
    }
  }

  ## Check covariate–tag time overlap -----------------
  ## t2index() uses findInterval(..., rightmost.closed = TRUE, left.open = TRUE),
  ## which returns 0 only when t < min(time_cov). Tags above the covariate range
  ## always get the last valid index, and a single-slice covariate (min == max)
  ## is accessible for any t >= that value. So we only need to warn when ALL tag
  ## times fall strictly below the covariate's minimum time.
  if (!is.null(res$time_cov) && !is.null(res$tags) && nrow(res$tags) > 0) {
    tag_max <- max(res$tags$t, na.rm = TRUE)
    for (i in seq_along(res$time_cov)) {
      cov_min <- min(res$time_cov[[i]], na.rm = TRUE)
      if (tag_max < cov_min) {
        tag_min <- min(res$tags$t, na.rm = TRUE)
        cov_max <- max(res$time_cov[[i]], na.rm = TRUE)
        warning(
          "All tag times are below the minimum time of covariate cov[[", i,
          "]]: tags span [", signif(tag_min, 5), ", ", signif(tag_max, 5),
          "], covariate spans [", signif(cov_min, 5), ", ",
          signif(cov_max, 5), "]. ",
          "The covariate will be inaccessible during fitting (t2index returns ",
          "0 for all observations), resulting in a zero gradient and no ",
          "parameter movement. Fix: ensure both use the same time system, ",
          "e.g. call prep_cov(..., date_decimal = TRUE) and ",
          "setup_data(..., shift_tref = TRUE).",
          call. = FALSE
        )
      }
    }
  }


  ## Splines ------------------------------------------
  if (!is.null(res$cov)) {
    res$time_spline <- lapply(seq_along(res$cov), function(x) 0)
  } else {
    res$time_spline <- list()
  }


  res$knots_tax <- knots_tax
  res$knots_dif <- knots_dif

  cov_obs <- cov
  if (is.null(res$knots_tax) && !is.null(cov)) {
    tmp <- sapply(cov,
                  function(x)
                    quantile(as.numeric(x),
                             get_pretty_probs(3),
                             na.rm = TRUE))
    res$knots_tax <- matrix(as.numeric(tmp), nrow = 3, ncol = length(cov))
  }

  if (length(res$knots_tax) != 0 && any(apply(res$knots_tax, 2, duplicated))) {
    warning("Some knots are the same! This will likely give an error!")
  }

  if (is.null(res$knots_dif) && !is.null(cov)) {
    tmp <- sapply(cov,
                  function(x)
                    quantile(as.numeric(x),
                             get_pretty_probs(1),
                             na.rm = TRUE))
    res$knots_dif <- matrix(as.numeric(tmp), nrow = 1, ncol = length(cov))
  }


  ## Prediction ---------------------------------------
  pred <- list()
  pred$time <- seq(trange[1], trange[2], length.out = 10)
  if (!is.null(res$cov)) pred$cov <- sapply(res$cov,
                                        function(x) seq(min(x, na.rm = TRUE),
                                          max(x, na.rm = TRUE),
                                          length.out = 100))
  if (!is.null(res$grid)) pred$grid <- res$grid
  res$pred <- pred

  ## Other variables -----------------------------------
  res$eps <- 0.000001
  res$var_init_kf <- 1e-6
  res$log2steps <- 0
  res$min_dt <- 0.1
  res$dt <- NULL ## if specified than equal ts for KF
  res$p_init <- c(0,0)

  ## Return
  res <- .add_class(res, "admove_data")
  res <- add_sref(res, master_sref)
  res <- add_tref(res, master_tref)

  return(res)
}



##' Summarise admove data
##'
##' @description Summarise data of any `admove` object
##'
##' @param object an object of class `admove_data` (created by `setup_data`) or an
##'   object containing such an object (`admove_data`, `admove_sim`, or
##'   `admove`).
##' @param ... Additional arguments
##'
##' @return Nothing.
##'
##' @examples
##'
##' summarise_data(skjepo$sim$dat)
##'
##' @name summarise_data
##' @export
summarise_data <- function(object, ...) {
  x <- object

  if(inherits(x, "admove_sim")) {
    dat <- x$dat
  } else if(inherits(x, "admove")) {
    dat <- x$dat
  } else if(inherits(x, "admove_data")) {
    dat <- x
  } else stop("Please provide an object of class 'admove_data' or an object containing an such an object (e.g. admove_sim, admove).")

  summarise_grid(dat$grid)
  cat("\n")
  summarise_cov(dat$cov)
  cat("\n")
  summarise_tags(dat$tags)

  invisible(dat)
}




##' Plot components of an `admove_data` object
##'
##' @description
##' Create summary plots for the main components of an object of class
##' `admove_data`. Depending on which components are present, the function plots
##' the spatial grid, covariate fields, and tag data in a multi-panel layout.
##'
##' @param x An object of class `admove_data`, as returned by [setup_data()].
##' @param auto_layout Logical; if `TRUE`, the plotting layout and graphical
##'   parameters are set automatically. Default is `TRUE`.
##' @param ... Additional arguments passed to the underlying plotting functions,
##'   including [plot_grid()], [plot_cov()], and [plot_tags()].
##'
##' @details
##' The function inspects `x` and plots all available data components. If
##' present, the grid is plotted first, followed by each covariate field, and
##' then tag data split by tag type:
##' \describe{
##'   \item{`"d"`}{Archival tags.}
##'   \item{`"s"`}{Mark-resight tags.}
##'   \item{`"c"`}{Conventional tags.}
##' }
##'
##' If `auto_layout = TRUE`, panels are arranged automatically using
##' [n2mfrow()]. Panel labels are added with [add_lab()].
##'
##' @return
##' Invisibly returns `NULL`. Called for its side effect of producing plots.
##'
##' @name plot_data
##' @export
plot_data <- function(x,
                      auto_layout = TRUE,
                      ...) {

  .check_class(x, "admove_data")

  if(auto_layout){
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar))
    n <- sum(!sapply(x[c("grid","cov","tags")], is.null))
    if (!is.null(x$cov)) {
      n <- n + length(x$cov) - 1
    }
    par(mfrow = n2mfrow(n, asp = 2), mar = c(4,4,1,1), oma = c(1,1,1,1))
  }

  i = 1
  if(!is.null(x$grid)){
    plot_grid(x$grid, auto_layout = FALSE, main = "", ...)
    add_lab(LETTERS[i])
    i = i + 1
  }
  if(!is.null(x$cov)){
    for (j in 1:length(x$cov)) {
      plot_cov(x$cov, i = j, select = 1, auto_layout = FALSE,
               main = "", ...)
      add_lab(LETTERS[i])
      i = i + 1
    }
  }
  if(!is.null(x$tags)){
    if (any(x$tags$tag_type == "d")) {
      dtags <- x$tags[x$tags$tag_type == "d",]
      plot_tags(dtags, auto_layout = FALSE, main = "", ...)
      add_lab(LETTERS[i])
      i = i + 1
    }
  }
  if(!is.null(x$tags)){
    if (any(x$tags$tag_type == "s")) {
      stags <- x$tags[x$tags$tag_type == "s",]
      plot_tags(stags, auto_layout = FALSE, main = "", ...)
      add_lab(LETTERS[i])
      i = i + 1
    }
  }
  if (!is.null(x$tags)) {
    if (any(x$tags$tag_type == "c")) {
      ctags <- x$tags[x$tags$tag_type == "c",]
      plot_tags(ctags, auto_layout = FALSE, main = "", ...)
      add_lab(LETTERS[i])
      i = i + 1
    }
  }

}






## s3 methods ------------------------------------------------------------------------

##' @rdname plot_data
##' @export
plot.admove_data <- function(x, ...) {
  plot_data(x, ...)
}


##' @method summary admove_data
##' @rdname summarise_data
##' @export
summary.admove_data <- function(object, ...) {
  summarise_data(object, ...)
}


##' @rdname print-admove
##' @method print admove_data
##' @export
print.admove_data <- function(x, ...) {
  tmp <- x
  attributes(tmp) <- NULL
  NextMethod("print", tmp, ...)
}
