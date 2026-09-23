
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
##'   [prep_ctags()], [prep_dtags()], or [prep_stags()]. Several tag objects can
##'   be supplied combined with `c(dtags, ctags)` or as a list,
##'   `list(dtags, ctags)`; both are merged with [combine_tags()].
##' @param trange Optional numeric vector of length two giving the model time
##'   range. If `NULL`, the time range is inferred from available tags and
##'   covariates.
##' @param knots_tax Optional matrix of spline knots for the taxis preference
##'   functions, with knots in rows and one column per covariate. If `NULL`,
##'   `n_knots_tax` knots are placed at covariate quantiles (see `knots_from`).
##' @param knots_dif Optional matrix of spline knots for the diffusion
##'   preference functions, with knots in rows and one column per covariate. If
##'   `NULL`, `n_knots_dif` knots are placed at covariate quantiles (see
##'   `knots_from`).
##' @param n_knots_tax Number of default knots per covariate for the taxis
##'   preference functions. Default is `3`. Ignored if `knots_tax` is supplied.
##' @param n_knots_dif Number of default knots per covariate for the diffusion
##'   preference functions. Default is `1`, i.e. constant diffusion. Ignored if
##'   `knots_dif` is supplied.
##' @param knots_from Where the default knots are taken from: `"tags"`
##'   (default) uses the covariate values at the tag observations, `"cov"` the
##'   whole covariate field. Without tags, the field is used either way. Ignored
##'   for a knot matrix that is supplied.
##' @param fill_na Number of rings of missing covariate cells to fill next to
##'   non-missing cells with [fill_cov()] before the grid and the tags are
##'   checked against the covariates. `1` fills the cells touching data, e.g.
##'   the coastline; `Inf` fills every reachable cell. `"grid"` fills as many
##'   rings as it takes to keep every grid cell, i.e. until the covariates can
##'   be interpolated at all grid cell centres (one ring without a grid). This
##'   matters when grid cells are larger than covariate cells: the centre of a
##'   coastal grid cell can then lie several covariate cells inland. Default
##'   `0`: no filling.
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
##' If spline knots are not supplied, default knots are placed at quantiles of
##' the covariate values the tags experienced (`knots_from = "tags"`): each
##' covariate is interpolated at the tag observations as in the likelihood, and
##' each tag counts equally, however many observations it has, so a few
##' archival tags do not set the knots for all tags. Knots from the whole field
##' (`knots_from = "cov"`, or when there are no tags) often lie in parts of the
##' covariate range no tag visits, where the preference function is not
##' informed by the data. Use [cov_at_tags()] to inspect the values.
##'
##' The number of knots is `n_knots_tax` per covariate for taxis (default
##' three: the 5%, 50% and 95% quantiles) and `n_knots_dif` for diffusion
##' (default one: the median, i.e. constant diffusion). Two knots
##' are placed at the 25% and 75% quantiles, four at the 5%, 30%, 70% and 95%
##' quantiles, and five or more evenly between the 5% and 95% quantiles. One
##' knot gives a constant function and two knots a linear one (a natural cubic
##' spline through two points is a straight line). All covariates share the
##' same number of knots. The spline coefficients created by [default_par()]
##' are sized from these knots, so set the number of knots here rather than
##' editing `dat$knots_tax` or `dat$knots_dif` afterwards.
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
##' ## five taxis knots and a linear diffusion function
##' dat5 <- setup_data(
##'   grid = grid,
##'   cov = cov,
##'   tags = c(dtags, ctags),
##'   n_knots_tax = 5,
##'   n_knots_dif = 2,
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
                       n_knots_tax = 3,
                       n_knots_dif = 1,
                       knots_from = c("tags", "cov"),
                       fill_na = 0,
                       sref = NULL,
                       tref = NULL,
                       transform_sref = FALSE,
                       shift_tref = FALSE,
                       verbose = TRUE) {

  n_knots_tax <- .check_n_knots(n_knots_tax, "n_knots_tax")
  n_knots_dif <- .check_n_knots(n_knots_dif, "n_knots_dif")
  knots_from <- match.arg(knots_from)
  if (!identical(fill_na, "grid")) fill_na <- .check_n_rings(fill_na, "fill_na")

  res <- list()

  if (!is.null(cov)) cov <- .make_cov_list(cov)

  ## accept a list of tag objects as well as c(dtags, ctags)
  if (!is.null(tags) && !inherits(tags, "admove_tags")) {
    if (is.list(tags) && !is.data.frame(tags)) {
      tags <- combine_tags(tags)
    } else {
      stop("'tags' must be an 'admove_tags' object (see prep_tags()), ",
           "or several combined with c() or list().")
    }
  }

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
  cells_cov_na <- NULL
  if (!is.null(res$cov) && identical(fill_na, "grid")) {
    res$cov <- .fill_cov_for_grid(res$cov, res$grid, verbose)
  } else if (!is.null(res$cov) && fill_na > 0) {
    res$cov <- fill_cov(res$cov, n_rings = fill_na, verbose = verbose)
  }

  if (!is.null(res$cov)) {

    ## space
    xyranges_cov <- .get_cov_xyrange(res$cov)
    res$xrange_cov <- xyranges_cov$xr
    res$yrange_cov <- xyranges_cov$yr

    ## Check if any cov NA where needed for interpol given provided grid
    if (!is.null(grid)) {
      err <- .cov_na_cells(res$cov, res$grid, xyranges_cov$xr,
                           xyranges_cov$yr)

      if (length(err) > 0) {
        ## kept to explain the tags that check_tags() drops from these cells
        grid_before_cov <- res$grid
        cells_cov_na <- err

        message(length(err), " grid cell(s) removed because the covariate(s) can not be calculated there (lead to NA): ", .format_ids(err), ". Removing them in order to avoid problems during fitting later. \n")

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
  if (!is.null(tags)) {
    if (verbose && length(cells_cov_na) > 0) {
      .message_tags_in_cells(tags, grid_before_cov, cells_cov_na, fill_na)
    }
    res$tags <- check_tags(tags, res$grid)
  }

  ## Remove tag positions on NA covariate cells -------
  ## check_tags() only guards against NA *grid* cells. A tag can still sit on
  ## (or next to) a covariate cell that is NA, where the bilinear interpolation
  ## used in the likelihood returns NaN and causes NaN objective/gradient
  ## evaluations during minimisation. Evaluate each covariate at the tag
  ## positions, matching how the likelihood accesses them (only time slices with
  ## t2index() > 0 are used), and drop the offending entries.
  if (!is.null(res$cov) && !is.null(res$tags) && nrow(res$tags) > 0) {
    bad <- rep(FALSE, nrow(res$tags))
    for (v in .cov_at_tags(res)) {
      bad <- bad | (is.na(v) & attr(v, "slice") > 0)
    }
    if (any(bad)) {
      bad_ids <- unique(res$tags$id[bad])
      if (verbose) {
        message(sum(bad), " entr", if (sum(bad) == 1) "y" else "ies",
                " removed because the tag position falls where a covariate is NA",
                " (tag id", if (length(bad_ids) == 1) "" else "s", ": ",
                .format_ids(bad_ids), ").")
      }
      res$tags <- res$tags[!bad, , drop = FALSE]
      ## Drop tags left with a single observation (need release + recovery),
      ## matching the recovered-tags filter in check_tags().
      keep_id <- names(which(table(res$tags$id) > 1))
      res$tags <- res$tags[res$tags$id %in% keep_id, , drop = FALSE]
      ## dropping rows here can also remove candidate positions of an ambiguous
      ## recapture, so restore the probability normalisation (see check_tags())
      res$tags <- .renormalise_events(res$tags, verbose)
      if (nrow(res$tags) == 0) {
        stop("No tags remain after removing positions on NA covariate cells.")
      }
    }
  }

  ## Check covariate–tag time overlap -----------------
  ## t2index() returns 0 only when t < min(time_cov) and otherwise clamps: a tag
  ## above the covariate range silently gets the LAST slice. Both directions are
  ## therefore checked. Above the range some overshoot is normal -- slices are
  ## labelled by their start, so a tag in the last month sits up to one slice
  ## spacing beyond the last label -- so only a substantial overshoot is
  ## flagged, and a single-slice covariate (a climatology) is skipped since it
  ## is legitimately used for every time.
  if (!is.null(res$time_cov) && !is.null(res$tags) && nrow(res$tags) > 0) {
    tag_min <- min(res$tags$t, na.rm = TRUE)
    tag_max <- max(res$tags$t, na.rm = TRUE)
    for (i in seq_along(res$time_cov)) {
      tc <- res$time_cov[[i]]
      cov_min <- min(tc, na.rm = TRUE)
      cov_max <- max(tc, na.rm = TRUE)

      if (tag_max < cov_min) {
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
        next
      }

      if (length(tc) < 2) next

      tol <- stats::median(diff(sort(tc)), na.rm = TRUE)
      if (!is.finite(tol) || tol <= 0) tol <- 0

      above <- which(res$tags$t > cov_max + tol)
      if (length(above) > 0) {
        warning(
          length(above), " of ", nrow(res$tags), " tag observation",
          if (length(above) == 1) "" else "s",
          " lie beyond the last time slice of covariate cov[[", i,
          "]]: those tag times span [", signif(min(res$tags$t[above]), 8),
          ", ", signif(max(res$tags$t[above]), 8), "], covariate spans [",
          signif(cov_min, 5), ", ", signif(cov_max, 5), "]. ",
          "t2index() clamps them, so they are all evaluated against the LAST ",
          "covariate slice regardless of their date, silently and without ",
          "error. Fix: ensure the tags and the covariate use the same time ",
          "system, e.g. give the tags a real time reference in prep_tags() ",
          "via 'date_origin' / 'date_format' / 'date_decimal', or extend the ",
          "covariate in time.",
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
  res$knots_from <- c(tax = "user", dif = "user")

  ## Default knots from the covariate values at the tag observations, one
  ## weight per tag; see dev/code_notes.org, "Default knot placement".
  ## Computed after the pruning above, so dropped positions do not count.
  if ((is.null(knots_tax) || is.null(knots_dif)) && !is.null(cov)) {
    knot_vals <- lapply(cov, function(x) as.numeric(unclass(x)))
    knot_w <- NULL
    from <- "cov"
    if (knots_from == "tags" && !is.null(res$tags) && nrow(res$tags) > 0) {
      at_tags <- .cov_at_tags(res)
      ok <- vapply(at_tags, function(v) any(is.finite(v)), logical(1))
      if (verbose && any(!ok)) {
        message("No tag observation lies within ",
                .knot_cov_labels(cov, which(!ok)),
                ", so its default knots are taken from the whole field.")
      }
      if (any(ok)) {
        w <- .tag_weights(res$tags)
        knot_w <- rep(list(NULL), length(cov))
        knot_vals[ok] <- lapply(at_tags[ok], as.numeric)
        knot_w[ok] <- list(w)
        from <- "tags"
      }
    }
    names(knot_vals) <- names(cov)
    if (is.null(knots_tax)) {
      res$knots_tax <- .default_knots(knot_vals, n_knots_tax, "n_knots_tax",
                                      knot_w)
      res$knots_from[["tax"]] <- from
    }
    if (is.null(knots_dif)) {
      res$knots_dif <- .default_knots(knot_vals, n_knots_dif, "n_knots_dif",
                                      knot_w)
      res$knots_from[["dif"]] <- from
    }
  }
  if (!is.null(knots_tax)) .warn_duplicated_knots(res$knots_tax, "knots_tax")
  if (!is.null(knots_dif)) .warn_duplicated_knots(res$knots_dif, "knots_dif")


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
  .summarise_knots(dat)

  invisible(dat)
}


## Number of knots and where the defaults came from (older data objects do not
## record the source).
.summarise_knots <- function(dat) {

  if (length(dat$knots_tax) == 0 && length(dat$knots_dif) == 0) {
    return(invisible(NULL))
  }
  src <- c(tags = "covariate at tags", cov = "covariate field",
           user = "supplied")
  line <- function(knots, from) {
    paste0(NROW(knots), " per covariate",
           if (!is.null(from)) paste0(", ", src[[from]]))
  }

  cat("<knots>\n")
  cat("  taxis:      ", line(dat$knots_tax, dat$knots_from[["tax"]]), "\n",
      sep = "")
  cat("  diffusion:  ", line(dat$knots_dif, dat$knots_from[["dif"]]), "\n",
      sep = "")

  invisible(NULL)
}



##' Covariate values at the tag observations
##'
##' @description
##' Interpolate each covariate at the positions and times of the tag
##' observations, in the same way as the likelihood does. Use it to see which
##' part of the covariate range the tags actually experienced, e.g. to choose
##' spline knots; by default, [setup_data()] places the knots at quantiles of
##' these values.
##'
##' @param object An object of class `admove_data` (created by [setup_data()])
##'   or an object containing one (`admove_sim` or `admove`).
##'
##' @return A named list with one numeric vector per covariate and one value per
##'   row of the tag data (`dat$tags`). Values are `NA` where the observation
##'   lies before the first time slice of the covariate (such observations are
##'   not evaluated in the likelihood).
##'
##' @examples
##' vals <- cov_at_tags(skjepo$sim)
##' lapply(vals, quantile, probs = c(0.05, 0.5, 0.95), na.rm = TRUE)
##'
##' @export
cov_at_tags <- function(object) {

  if (inherits(object, "admove_sim") || inherits(object, "admove")) {
    dat <- object$dat
  } else if (inherits(object, "admove_data")) {
    dat <- object
  } else {
    stop("Please provide an object of class 'admove_data' or an object ",
         "containing such an object (e.g. admove_sim, admove).", call. = FALSE)
  }
  if (length(dat$cov) == 0 || is.null(dat$tags) || nrow(dat$tags) == 0) {
    stop("The data contain no covariates or no tags.", call. = FALSE)
  }

  lapply(.cov_at_tags(dat), function(v) {
    attr(v, "slice") <- NULL
    v
  })
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
    ## one panel for the grid, one per covariate, and one per tag type
    n <- as.integer(!is.null(x$grid))
    if (!is.null(x$cov)) {
      n <- n + length(.make_cov_list(x$cov))
    }
    if (!is.null(x$tags)) {
      n <- n + sum(c("d", "s", "c") %in% x$tags$tag_type)
    }
    n <- max(n, 1L)
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


## Helpers ----------------------------------------------------------------------

## Validate a number of spline knots: a single whole number of at least one.
.check_n_knots <- function(n, name) {

  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 1 ||
        n != round(n)) {
    stop("'", name, "' must be a single whole number of at least 1.",
         call. = FALSE)
  }

  as.integer(n)
}


## Default spline knots: `n` quantiles of each covariate, as a matrix with knots
## in rows and one column per covariate. `vals` is a list with one numeric
## vector per covariate (the field or the values at the tags) and `w` an
## optional list of matching weights (NULL entries: unweighted). A covariate
## with too few distinct values gives repeated knots, which the natural spline
## cannot pass through, so say which covariate it is instead of failing later
## inside the likelihood.
.default_knots <- function(vals, n, name, w = NULL) {

  knots <- vapply(seq_along(vals),
                  function(i) .wquantile(vals[[i]], w[[i]], get_pretty_probs(n)),
                  numeric(n))
  knots <- matrix(knots, nrow = n, ncol = length(vals))

  dup <- which(apply(knots, 2, function(k) any(duplicated(k))))
  if (length(dup) > 0) {
    warning("With ", name, " = ", n, ", some default knots coincide for ",
            .knot_cov_labels(vals, dup), " (too few distinct covariate values). ",
            "This will likely give an error! Use fewer knots or supply the ",
            "knot matrix directly.", call. = FALSE)
  }

  knots
}


## Indices (rows of grid$xygrid) of the grid cells at whose centre a covariate
## cannot be interpolated in some time slice. setup_data() removes these cells;
## the same interpolation as in the likelihood decides.
.cov_na_cells <- function(cov, grid, xr, yr) {

  err <- NULL
  for (i in seq_along(cov)) {
    covi <- cov[[i]]
    for (j in seq_len(dim(covi)[3])) {
      liv <- RTMB::interpol2Dfun(covi[,,j],
                                 xlim = round(xr[i,], 5),
                                 ylim = round(yr[i,], 5),
                                 R = 1)
      tmp <- liv(round(grid$xygrid[,1], 5), round(grid$xygrid[,2], 5))
      err <- c(err, which(is.na(tmp)))
    }
  }

  sort(unique(err))
}


## setup_data(fill_na = "grid"): fill one ring at a time until no grid cell is
## lost to a missing covariate, or nothing is left to fill.
.fill_cov_for_grid <- function(cov, grid, verbose) {

  n_na <- function(cov) sum(vapply(cov, function(x) sum(is.na(x)), numeric(1)))

  if (is.null(grid)) {
    if (verbose) message("fill_na = \"grid\" without a grid: filling one ring.")
    return(fill_cov(cov, n_rings = 1, verbose = verbose))
  }

  xy <- .get_cov_xyrange(cov)
  n0 <- n_na(cov)
  k <- 0
  repeat {
    err <- .cov_na_cells(cov, grid, xy$xr, xy$yr)
    if (length(err) == 0) break
    before <- n_na(cov)
    cov <- fill_cov(cov, n_rings = 1, verbose = FALSE)
    if (n_na(cov) == before) break
    k <- k + 1
  }

  if (verbose) {
    message("fill_na = \"grid\": filled ", k, " ring", if (k == 1) "" else "s",
            " (", n0 - n_na(cov), " of ", n0, " NA covariate cells)",
            if (length(err) == 0) " so that no grid cell is lost" else
              paste0("; ", length(err), " grid cell(s) have no covariate data ",
                     "within reach"), ".")
  }

  cov
}


## Explain tag entries that check_tags() is about to drop because their grid
## cell was just removed for a missing covariate (its message only says "NA
## grid cell", which reads as if the input grid were masked there).
.message_tags_in_cells <- function(tags, grid, cells, fill_na) {

  ic <- grid$celltable[cbind(as.integer(cut(tags$x, grid$xgr)),
                             as.integer(cut(tags$y, grid$ygr)))]
  hit <- which(ic %in% cells)
  if (length(hit) == 0) return(invisible(NULL))

  ids <- unique(tags$id[hit])
  hint <- if (identical(fill_na, "grid")) "" else
    paste0(" Use setup_data(fill_na = \"grid\") to fill the covariate",
           " until these cells are kept.")
  message(length(hit), " tag entr", if (length(hit) == 1) "y lies" else "ies lie",
          " in the grid cell(s) removed because a covariate is NA at the cell",
          " centre (tag id", if (length(ids) == 1) "" else "s", ": ",
          .format_ids(ids), "); they are dropped below.", hint)

  invisible(NULL)
}


## Weighted quantiles that equal stats::quantile(type = 7) for equal weights:
## the sorted values sit at cumulative-weight positions running from 0 (first)
## to 1 (last), and the quantiles are interpolated linearly between them.
.wquantile <- function(x, w, probs) {

  if (is.null(w)) w <- rep(1, length(x))
  keep <- is.finite(x) & is.finite(w) & w > 0
  x <- x[keep]
  w <- w[keep]
  if (length(x) == 0) return(rep(NA_real_, length(probs)))
  if (length(x) == 1) return(rep(x, length(probs)))

  o <- order(x)
  x <- x[o]
  w <- w[o]
  n <- length(x)
  pos <- cumsum(w) - w / 2 - w[1] / 2
  pos <- pos / (sum(w) - w[1] / 2 - w[n] / 2)

  stats::approx(pos, x, xout = probs, rule = 2, ties = "ordered")$y
}


## One weight per tag row so that every tag counts equally in the default
## knots, however many observations it has. The candidate positions of an
## ambiguous recapture (see check_tags()) share one observation's weight
## according to their probabilities.
.tag_weights <- function(tags) {

  id <- as.character(tags$id)
  if (all(c("event", "prob") %in% colnames(tags))) {
    ev <- paste(id, tags$event)
    n_obs <- tapply(ev, id, function(e) length(unique(e)))
    tags$prob / as.numeric(n_obs[id])
  } else {
    1 / as.numeric(table(id)[id])
  }
}


## Covariate values at the tag observations, interpolated exactly as the
## likelihood does (RTMB::interpol2Dfun, R = 1, on the slice t2index() picks).
## One vector per covariate, one value per row of `dat$tags`, with the slice
## index as attribute "slice": NA where the covariate is NA (the NA pruning in
## setup_data() relies on this) and where slice == 0 (tag time before the first
## slice, never evaluated in the likelihood).
.cov_at_tags <- function(dat) {

  tags <- dat$tags
  xr <- dat$xrange_cov
  yr <- dat$yrange_cov

  res <- lapply(seq_along(dat$cov), function(i) {
    covi <- dat$cov[[i]]
    it <- as.integer(t2index(tags$t, dat$time_cov[[i]]))
    v <- rep(NA_real_, nrow(tags))
    for (j in sort(unique(it[it > 0]))) {
      rows <- which(it == j)
      liv <- RTMB::interpol2Dfun(covi[,,j],
                                 xlim = round(xr[i,], 5),
                                 ylim = round(yr[i,], 5),
                                 R = 1)
      v[rows] <- liv(round(tags$x[rows], 5), round(tags$y[rows], 5))
    }
    attr(v, "slice") <- it
    v
  })
  names(res) <- names(dat$cov)

  res
}


## Warn about repeated knots in a user-supplied knot matrix.
.warn_duplicated_knots <- function(knots, name) {

  if (length(knots) == 0 || !is.matrix(knots)) return(invisible(NULL))

  dup <- which(apply(knots, 2, function(k) any(duplicated(k))))
  if (length(dup) > 0) {
    warning("Some knots in '", name, "' are the same (column(s) ",
            paste(dup, collapse = ", "), ")! This will likely give an error!",
            call. = FALSE)
  }

  invisible(NULL)
}


.knot_cov_labels <- function(cov, i) {
  nms <- names(cov)
  lab <- if (is.null(nms)) rep("", length(cov)) else nms
  lab <- ifelse(is.na(lab) | lab == "", paste0("covariate ", seq_along(cov)),
                paste0("covariate '", lab, "'"))
  paste(lab[i], collapse = ", ")
}
