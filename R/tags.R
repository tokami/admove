
## Main functions ---------------------------------------------------------------------

##' Prepare tagging data for admove
##'
##' @description
##' `prep_tags()` converts raw tagging data into a standardised object of class
##' `admove_tags` used throughout \emph{admove}. It supports data-storage
##' (archival) tags, mark-resight tags, and mark-recapture tags, and accepts
##' input either in long format, wide format, or as a list split by tag.
##'
##' The function harmonises variable names, optionally converts date variables
##' to numeric model time, adds missing tag identifiers when needed, and attaches
##' spatial and temporal reference information.
##'
##' @param x Tagging data. Can be a data frame, a list of data frames
##'   (typically one per tag), or an existing object of class `admove_tags`.
##' @param tag_type Character string specifying the tag type:
##'   `"d"` for data-storage (archival) tags,
##'   `"s"` for mark-resight tags, and
##'   `"c"` for mark-recapture tags.
##' @param names Named character vector giving the column names in `x`.
##'   For long format, provide at least `c(t = "...", x = "...", y = "...",
##'   id = "...")`. For wide format, provide
##'   `c(t0 = "...", t1 = "...", x0 = "...", y0 = "...", x1 = "...", y1 = "...")`.
##'   The order does not matter, but the vector must be named.
##' @param date_decimal Logical; if `TRUE`, interpret the time variable as a
##'   decimal year and convert it to model time. Default: `FALSE`.
##' @param date_format Optional character string passed to [as.Date()] to parse
##'   character dates. Default: `NULL`.
##' @param date_origin Optional origin passed to [as.Date()] when times are
##'   stored numerically. Default: `NULL`.
##' @param keep_only_recaptured Logical; if `TRUE`, only keep tags with at least
##'   two observations. Currently mainly relevant for mark-recapture-style data.
##'   Default: `TRUE`.
##' @param tz Time zone used when converting dates. Default: `"UTC"`.
##' @param sref Optional spatial reference information to attach to the returned
##'   object.
##' @param tref Optional time reference information to attach to the returned
##'   object.
##' @param transform_sref Logical; if `TRUE`, transform coordinates to `sref`
##'   when possible. Default: `FALSE`.
##' @param shift_tref Logical; if `TRUE`, shift times to `tref` when possible.
##'   Default: `FALSE`.
##' @param verbose Logical; if `TRUE`, print informative messages. Default:
##'   `TRUE`.
##'
##' @return
##' A data frame of class `admove_tags` with standardised columns such as
##' `t`, `x`, `y`, `id`, `tag_type`, and `use`, plus any additional columns
##' provided in the input.
##'
##' @details
##' `prep_tags()` accepts three common input structures:
##'
##' \itemize{
##'   \item a long-format data frame with one row per observation,
##'   \item a wide-format data frame with release and recapture columns
##'     (`t0`, `t1`, `x0`, `y0`, `x1`, `y1`), or
##'   \item a list of data frames, usually one element per tag.
##' }
##'
##' If no tag identifier is supplied for list input, the list order is used to
##' create an `id` column automatically.
##'
##' @examples
##' ## prepare data-storage tags
##' dtags <- prep_tags(
##'   skjepo$dtags,
##'   tag_type = "d",
##'   names = c(t = "time", x = "mptlon", y = "mptlat"),
##'   date_origin = "1899-12-30"
##' )
##'
##' ## prepare mark-recapture tags
##' ctags <- prep_tags(
##'   skjepo$ctags,
##'   tag_type = "c",
##'   names = c(
##'     t0 = "date_time", t1 = "date_caught",
##'     x0 = "rel_lon",   x1 = "recap_lon",
##'     y0 = "rel_lat",   y1 = "recap_lat"
##'   ),
##'   date_origin = "1899-12-30"
##' )
##'
##' @name prep_tags
##' @export
prep_tags <- function(x,
                      tag_type = NULL,
                      names = NULL,
                      date_decimal = FALSE,
                      date_format = NULL,
                      date_origin = NULL,
                      keep_only_recaptured = TRUE,
                      tz = "UTC",
                      sref = NULL,
                      tref = NULL,
                      transform_sref = FALSE,
                      shift_tref = FALSE,
                      verbose = TRUE) {

  ## general stuff
  req_wide <- c("t0","t1","x0","y0","x1","y1")
  req_long <- c("t","x","y","id")
  colis <- colnames(x)

  if (inherits(x, "admove_tags")) {
    cols <- colnames(x)
    if (!all(req_long %in% colis)) stop("Not all required variables in the tagging data set. See details and examples in ?prep_tags for more information.")
    stopifnot(any(colis == "tag_type"))
    if (!any(colis == "use")) {
      x$use <- 1
    }
    return(x)
  }

  if (!tag_type %in% c("d","s","c")) stop("Tag type ('tag_type') has to one of the following letters: 'd' = data-storage tags, 's' = mark-resight tags, and 'c' = mark-recapture tags.")

  if (is.null(names)) stop("Please provide the names of the columns with the required minimum information: t = time, x = x position, y = y position, and id = identifier for multiple tags. See details and examples in ?prep_tags for more information.")

  if (is.null(names(names))) stop("The 'names' vector has to be a named vector, e.g. names = c(t = '...', x = '...', y = '...', id = '...'). See details and examples in ?prep_tags for more information.")

  flag_id <- ifelse(any(names(names) == "id"), TRUE, FALSE)
  flag_wide <- ifelse(any(names(names) %in% req_wide), TRUE, FALSE)

  if (flag_wide) {
    req <- req_wide
  } else {
    req <- req_long
  }

  if (inherits(x, "list")) {
    if (!flag_id) {
      if (verbose) message("ID not specified, using order of list elements.")
      x <- lapply(1:length(x), function(i) cbind(x[[i]], id = i))
      names <- c(names, id = "id")
    }
    x <- do.call(rbind, x)
  } else if (!flag_wide && !flag_id) stop("Not clear which entries belong together, please provide either an 'id' column, the input data as a list by tag, or the long format, in which each row corresponds to a tag. See details and examples in ?prep_tags for more information.")

  if (!all(req %in% names(names))) stop("Not all required variables provided. Please provided names = c(t = '...', x = '...', y = '...', id = '...') for the long format or names = c(t0 = '...', t1 = '...', x0 = '...', y0 = '...', x1 = '...', y1 = '...') for the wide format (mostly relevant for mark-recapture tags). (The order doesn't matter). See details and examples in ?prep_tags for more information.")

  ## any column names missing
  col_names <- colnames(x)
  idx <- which(is.na(col_names))
  if (length(idx) > 0) {
    if (verbose) message("Some columns names missing! Setting names for now, but lease check your data.")
    col_names[idx] <- paste0("unkown",seq_along(idx))
  }
  colnames(x) <- col_names

  x_in <- x
  x <- x[,names[req]]
  col_names <- req

  idx <- which(!(colnames(x_in) %in% names[req]))
  if (length(idx) > 0) {
    x <- cbind(x, x_in[,idx])
    col_names <- c(col_names, colnames(x_in)[idx])
  }
  colnames(x) <- col_names

  if (flag_id) colnames(x)[colnames(x) == names["id"]] <- "id"

  if (flag_wide) {
    x <- ctags_wide_2_long(x)
  }

  ## Convert dates
  if (!is.null(date_origin) ||
        !is.null(date_format) ||
         isTRUE(date_decimal)) {

    if (is.null(date_format) && !is.null(date_origin)) {
      dati <- as.Date(x$t,
                      origin = date_origin,
                      tz = tz)

    } else if (is.null(date_origin) && !is.null(date_format)) {
      dati <- as.Date(x$t,
                      format = date_format,
                      tz = tz)

    } else if (isFALSE(date_decimal)) {
      dati <- as.Date(x$t,
                      origin = date_origin,
                      format = date_format,
                      tz = tz)
    } else {
      dati <- .decimal_year_2_date(as.numeric(x$t), tz = tz)
    }

    x$t <- date_2_time(dati, tref)
    tref <- create_tref(attr(x$t, "tref")[["origin"]],
                        attr(x$t, "tref")[["units"]])
    if (isTRUE(attr(x$t, "tref")[["inferred"]])) {
      if (verbose) message("tref (time origin and units) was inferred from dates. Please check and adjust if needed.")
    }
  }

  if (all(!is.numeric(x$t))) stop("The time column is not numeric. Please provide the time as a decimal date or use the 'date_format' and/or 'date_origin' argument to convert the date into a decimal date (see ?as.Date).")

  x_list <- split(x, x$id)

  ## combine
  res <- do.call(rbind, x_list)
  rownames(res) <- NULL

  ## tag type
  res$tag_type <- .get_tag_type(tag_type)

  ## use column to deactivate
  if (!any(colnames(res) == "use")) res$use <- 1

  ## class and sref & tref
  res <- .add_class(res, "admove_tags")
  res <- add_sref(res, sref, transform_crs = transform_sref)
  res <- add_tref(res, tref, shift_origin = shift_tref)

  return(res)
}




##' Check and standardise a tag object
##'
##' @description
##' `check_tags()` validates tagging data, removes invalid entries, optionally
##' checks whether observations fall within the spatial and temporal domain of a
##' grid or data object, and returns a cleaned object of class `admove_tags`.
##'
##' The function can be used on raw tagging data, on an `admove_tags` object,
##' or on higher-level \emph{admove} objects that contain tags.
##'
##' @param x Tagging data to check. Can be an object of class `admove_tags`, a
##'   data frame, a list of tag-specific data frames, or an object that contains
##'   tags such as `admove_data`, `admove_sim`, or `admove`.
##' @param grid Optional spatial grid used to remove observations outside the
##'   spatial domain and to assign grid cells.
##' @param dat Optional `admove_data` object used to check whether observations
##'   fall within the time domain.
##' @param conf Optional configuration list. If provided, only tag types enabled
##'   in `conf` are retained.
##' @param remove_non_recovered_tags Logical; if `TRUE`, remove tags with fewer
##'   than two observations. Default: `TRUE`.
##' @param verbose Logical; if `TRUE`, print informative messages about removed
##'   entries. Default: `TRUE`.
##'
##' @return
##' A cleaned object of class `admove_tags`.
##'
##' @details
##' The function performs several checks, including:
##' \itemize{
##'   \item removal of rows with missing required values in `t`, `x`, or `y`,
##'   \item removal of rows with `use = FALSE`,
##'   \item removal of observations outside the spatial domain of `grid`,
##'   \item removal of observations outside the temporal domain of `dat`, and
##'   \item optional removal of tags with fewer than two observations.
##' }
##'
##' Spatial and temporal reference information are attached to the returned
##' object from `grid` and `dat` when available.
##'
##' @export
check_tags <- function(x, grid = NULL, dat = NULL, conf = NULL,
                       remove_non_recovered_tags = TRUE,
                       verbose = TRUE) {


  if(!inherits(x, "admove_tags")) {
    tags <- x
  } else if(inherits(x, "admove_data")) {
    tags <- x$tags
    if (is.null(grid)) grid <- x$grid
    if (is.null(dat)) dat <- x
  } else if(inherits(x, "admove_sim")) {
    tags <- x$tags
    if (is.null(grid)) grid <- x$dat$grid
    if (is.null(dat)) dat <- x$dat
    if (is.null(conf)) conf <- x$conf
  } else if(inherits(x, "admove")) {
    tags <- dat$tags
    if (is.null(grid)) grid <- dat$grid
    if (is.null(dat)) dat <- x$dat
    if (is.null(conf)) conf <- x$conf
  } else {
    tags <- x
  }

  if (is.null(tags)) stop("No tags found!")

  flag_grid <- ifelse(is.null(grid), FALSE, TRUE)
  flag_dat <- ifelse(is.null(dat), FALSE, TRUE)


  ## as data.frame
  if (inherits(tags, "list")) {
    tags <- lapply(1:length(tags), function(i) {
      if (!any(colnames(tags[[i]]) == "id")) {
        tags[[i]]$id <- paste0(.get_random_id(1), "-", i)
      }
      tags[[i]]
    })
    tags <- do.call(rbind, tags)
  }

  if (!any(colnames(tags) == "id")) {
    tags$id <- .get_random_id(1)
    if (verbose) message("No id found. Assuming all entries belong to same tag.")
  }

  ## scale coords to grid
  if (flag_grid) {
    ## sref(tags)
    ## head(tags)
    ## grid
    ## scale_sref(tags, crs_scale(grid))
  }

  ## redo ic
  if (flag_grid) {
    tags$ic <- grid$celltable[cbind(as.integer(cut(tags$x, grid$xgr)),
                                    as.integer(cut(tags$y, grid$ygr)))]
  }

  ## remove NaN
  ind <- which(apply(tags[,c("t","x","y")], 1, function(x) any(is.na(x))))
  if (length(ind) > 0) {
    tags <- tags[-ind,]
    if (verbose) message(length(ind), " entries removed because required info (t, x, y) is NaN.")
  }

  ## tag identifier missing
  if (!any(colnames(tags) == "tag_type")) stop("No tag_type found. Please add a column called tag_type to the tags and specify if the tags are data-storage (tag_type = 1), mark-resight tags (tag_type = 2), or mark-recapture tags (tag_type = 3). ")

  tags$tag_type <- .get_tag_type(tags$tag_type)


  ## use missing
  if (!any(colnames(tags) == "use")) tags$use <- 1


  ## remove use = 0
  ind <- which(isFALSE(as.logical(tags$use)))
  if (length(ind) > 0) {
    tags <- tags[-ind,]
    if (verbose) message(length(ind), " entries removed because use = FALSE.")
  }


  ## update in ctmc missing
  if (!any(colnames(tags) == "update")) {
    tags$update <- TRUE
    tags$update[tags$tag_type %in% "c"] <- FALSE
  }

  if (any(tags$tag_type %in% "c") &&
        any(isTRUE(tags$update[tags$tag_type %in% "c"]))) {
    warning("You are combining mark-recapture tags (tag_type = c) and updating in CTMC. This is not recommended and can lead to unexpected results in CTMC. Consider setting update = FALSE for mark-recapture tags!")
  }


  ## outside of spatial domain
  if (flag_grid) {
    ind <- which(tags$x < grid$xrange[1] |
                   tags$x > grid$xrange[2] |
                     tags$y < grid$yrange[1] |
                     tags$y > grid$yrange[2])
    if (length(ind) > 0) {
      if (verbose) message(length(ind), " entries removed because outside of spatial dimensions.")
      tags <- tags[-ind,]
    }

    ## ind <- which(is.na(tags$ic) | is.na(grid$celltable[tags$ic]))
    ind <- which(is.na(tags$ic))
    if (length(ind) > 0) {
      tags <- tags[-ind,]
      message(length(ind), " entries removed because ic = NA or grid[ic] = NA.")
    }

  }


  ## outside of time domain
  if (flag_dat) {
    ind <- which(tags$t < dat$trange[1] |
                   tags$t > dat$trange[2])
    if(length(ind) > 0){
      if (verbose) message(length(ind), " entries removed because outside of time dimensions.")
      tags <- tags[-ind,]
    }
  }

  ## keep only recovered tags
  if (remove_non_recovered_tags) {
    tags_list <- split(tags, tags$id)
    tags_list <- tags_list[sapply(tags_list, nrow) > 1]
  }

  tags_out <- do.call(rbind, tags_list)
  if (is.null(tags_out) || nrow(tags_out) == 0) stop("No tags passed the checks!")

  ## strict if conf provided
  if (!is.null(conf)) {

    tags <- tags_out

    ## Combine all tags into a list
    if (!conf$use_dtags || is.null(tags) || !any(tags$tag_type == "d")) {
      dtags <- NULL
    } else {
      dtags0 <- tags[tags$tag_type == "d",]
      dtags <- split(dtags0, dtags0$id)
    }
    if (!conf$use_stags || is.null(tags) || !any(tags$tag_type == "s")) {
      stags <- NULL
    } else {
      stags0 <- tags[tags$tag_type == "s",]
      stags <- split(stags0, stags0$id)
    }

    if(conf$use_ctags && !is.null(tags) && any(tags$tag_type == "c")) {
      ctags0 <- tags[tags$tag_type == "c",]
      ctags <- split(ctags0, ctags0$id)
    } else {
      ctags <- NULL
    }

    tags_out <- c(dtags, stags, ctags)

    if (length(tags_out) < 1) stop("No tags provided or selected (conf$use_*tags). Please provide either conventional, data-logging, or mark-resight tags in admove_data().")

    tags_out <- do.call(rbind, tags_out)
    if (nrow(tags_out) == 0) stop("No tags passed the checks!")
  }


  ## return
  tags_out <- .add_class(tags_out, "admove_tags")
  if (flag_grid) {
    tags_out <- add_sref(tags_out, sref(grid))
  } else {
    sref(tags_out) <- create_sref()
  }
  if (flag_dat) {
    tags_out <- add_tref(tags_out, tref(dat))
} else {
    tref(tags_out) <- create_tref()
  }

  return(tags_out)
}


##' Combine `admove_tags` objects
##'
##' @description
##' `combine_tags()` combines multiple objects of class `admove_tags` into a
##' single `admove_tags` object. It is designed to support workflows such as
##' combining data-storage, mark-resight, and mark-recapture tags into one
##' unified data set.
##'
##' All inputs must inherit from `admove_tags` and must have compatible spatial
##' and temporal reference information.
##'
##' @param ... Objects of class `admove_tags`, or a single list containing such
##'   objects.
##' @param recursive Ignored. Included for compatibility with the generic
##'   [base::c()] interface.
##'
##' @return
##' A single object of class `admove_tags` containing all rows from the input
##' objects.
##'
##' @details
##' Before combining inputs, the function:
##' \itemize{
##'   \item removes `NULL` inputs,
##'   \item checks that all remaining inputs inherit from `admove_tags`,
##'   \item aligns columns across objects by filling missing columns with `NA`,
##'   \item checks that all objects have the same `sref`, and
##'   \item checks that all objects have the same `tref`.
##' }
##'
##' If spatial or temporal reference information differs across inputs, the
##' function stops with an error.
##'
##' @examples
##' ## tags <- combine_tags(dtags, stags, ctags)
##' ## tags <- c(dtags, ctags)
##'
##' @export
combine_tags <- function(..., recursive = FALSE) {
  dots <- list(...)
  if (length(dots) == 0) return(structure(data.frame(), class = c("admove_tags","data.frame")))

  ## Flatten if someone does c(list(dtags, ctags))
  if (length(dots) == 1 && is.list(dots[[1]]) && !inherits(dots[[1]], "data.frame")) {
    dots <- dots[[1]]
  }

  ## Keep only non-NULL
  dots <- Filter(Negate(is.null), dots)
  if (length(dots) == 0) return(structure(data.frame(), class = c("admove_tags","data.frame")))

  ## Validate inputs
  bad <- vapply(dots, function(x) !inherits(x, "admove_tags"), logical(1))
  if (any(bad)) {
    stop("All inputs to c.admove_tags must inherit from 'admove_tags'.",
         call. = FALSE)
  }

  ## Coerce to data.frame and ensure type column exists
  dfs <- lapply(dots, function(x) {
    df <- as.data.frame(x)
    if (!("tag_type" %in% names(df))) stop("tag_type column missing. Cannot savely combine tags.")
    df
  })

  ## Align columns across inputs
  all_cols <- unique(unlist(lapply(dfs, names), use.names = FALSE))
  dfs <- lapply(dfs, function(df) {
    miss <- setdiff(all_cols, names(df))
    if (length(miss)) df[miss] <- NA
    df[all_cols]
  })

  ## Row-bind
  res <- do.call(rbind, dfs)

  ## Clean up row names
  rownames(res) <- NULL
  res <- .add_class(res, "admove_tags")

  ## strict check for sref
  bad <- character(0)
  srefs <- lapply(dots, sref)
  checks <- sapply(srefs, function(x) sref_equal(x, srefs[[1]]))
  if (all(checks)) {
    res <- add_sref(res, srefs[[1]])
  } else {
    stop("The spatial reference infos of the tags are not the same. Check sref(...) of your tags and make sure they align.")
  }

  ## strict check for tref
  bad <- character(0)
  trefs <- lapply(dots, tref)
  checks <- sapply(trefs, function(x) tref_equal(x, trefs[[1]]))
  if (all(checks)) {
    res <- add_tref(res, trefs[[1]])
  } else {
    stop("The time reference infos of the tags are not the same. Check tref(...) of your tags and make sure they align.")
  }


  return(res)
}



##' Get spatial and temporal ranges from tagging data
##'
##' @description
##' `get_dim_tags()` extracts the observed temporal range and the spatial ranges
##' in the `x` and `y` directions from tagging data.
##'
##' @param tags Tagging data as a data frame, an object of class `admove_tags`,
##'   or a list of tag-specific data frames.
##'
##' @return
##' A list with components:
##' \describe{
##'   \item{trange}{Range of observed times.}
##'   \item{xrange}{Range of observed `x` coordinates.}
##'   \item{yrange}{Range of observed `y` coordinates.}
##' }
##'
##' @export
get_dim_tags <- function(tags = NULL) {

  if (inherits(tags, "list")) {
    tags <- do.call(rbind, tags)
  }

  trange <- xrange <- yrange <- NULL
  if (!is.null(tags)) {
    trange <- range(tags[,"t"], na.rm = TRUE)
    xrange <- range(tags[,"x"], na.rm = TRUE)
    yrange <- range(tags[,"y"], na.rm = TRUE)
  }
  res <- list(trange = trange,
              xrange = xrange,
              yrange = yrange)
  return(res)
}


##' Convert wide conventional-tag data to long format
##'
##' @description
##' Internal helper to convert wide-format tag data with release and recapture
##' columns (`t0`, `t1`, `x0`, `x1`, `y0`, `y1`) into long format with columns
##' `t`, `x`, `y`, and `id`.
##'
##' @param x A data frame in wide format.
##'
##' @return
##' A data frame in long format with two rows per tag.
##'
##' @noRd
ctags_wide_2_long <- function(x) {

  ## return NULL if NULL
  if (is.null(x)) return(NULL)

  ## random id if missing
  if (!any(colnames(x) == "id")) {
    x$id <- paste0(.get_random_id(3), "-", 1:nrow(x))
  }

  other_cols <- which(!colnames(x) %in% c("t0","t1","x0","x1","y0","y1","id"))

  ## wide 2 long list
  ntags <- nrow(x)
  res_list <- vector("list", ntags)
  for (i in seq_len(ntags)) {
    tag <- data.frame(
      t = c(x[i,"t0"],
            x[i,"t1"]),
      x = c(x[i,"x0"],
            x[i,"x1"]),
      y = c(x[i,"y0"],
            x[i,"y1"]),
      id = rep(x[i, "id"],2)
    )
    if (length(other_cols) > 0) {
      tmp <- as.data.frame(x[rep(i,2), other_cols])
      colnames(tmp) <- colnames(x)[other_cols]
      tag <- cbind(tag, tmp)
    }
    res_list[[i]] <- tag
  }

  ## long dataframe
  res <- do.call(rbind, res_list)
  rownames(res) <- NULL
  return(res)
}



##' Convert a list of conventional tags to wide format
##'
##' @description
##' Internal helper to convert a list of conventional-tag records into a
##' wide-format data frame with one row per tag.
##'
##' @param x A list of tag-specific data frames.
##'
##' @return
##' A data frame with one row per tag and columns such as `t0`, `t1`, `x0`,
##' `x1`, `y0`, `y1`, `id`, and `tag_type`.
##'
##' @noRd
ctags_list_2_wide <- function(x) {

  if (is.null(x)) return(NULL)
  if (!is.list(x)) stop("x must be a list.")

  sapply(x, nrow)

  ntags <- length(x)
  if (ntags == 0L) return(NULL)

  ## preallocate
  res <- data.frame(
    id = rep(NA_character_, ntags),
    t0 = rep(NA_real_, ntags),
    t1 = rep(NA_real_, ntags),
    x0 = rep(NA_real_, ntags),
    x1 = rep(NA_real_, ntags),
    y0 = rep(NA_real_, ntags),
    y1 = rep(NA_real_, ntags),
    ## itrel = rep(NA_integer_, ntags),
    ## itrec = rep(NA_integer_, ntags),
    ## icrel = rep(NA_integer_, ntags),
    ## icrec = rep(NA_integer_, ntags),
    ## use = rep(NA, ntags),
    tag_type = rep(.get_tag_type("c"), ntags),
    stringsAsFactors = FALSE
  )

  ## helpers
  get_id <- function(tag, i) {
    if ("id" %in% names(tag)) {
      as.character(tag$id[1])
    } else if (!is.null(names(x)) && nzchar(names(x)[i])) {
      as.character(names(x)[i])
    } else {
      NA_character_
    }
  }

  has_all <- function(x, nms) all(nms %in% names(x))

  for (i in seq_len(ntags)) {
    tag <- x[[i]]
    if (!is.data.frame(tag)) stop(sprintf("Element %d is not a data.frame.", i))

    ## accept either:
    ## (A) long: 2 rows with t/x/y
    ## (B) already-wide: 1 row with t0/t1/x0/x1/y0/y1
    if (nrow(tag) == 1L && has_all(tag, c("t0","t1","x0","x1","y0","y1"))) {
      res$id[i] <- get_id(tag, i)
      res$t0[i] <- tag$t0[1]; res$t1[i] <- tag$t1[1]
      res$x0[i] <- tag$x0[1]; res$x1[i] <- tag$x1[1]
      res$y0[i] <- tag$y0[1]; res$y1[i] <- tag$y1[1]

      ## if (has_all(tag, c("itrel","itrec"))) {
      ##   res$itrel[i] <- tag$itrel[1]
      ##   res$itrec[i] <- tag$itrec[1] }
      ## if (has_all(tag, c("icrel","icrec"))) {
      ##   res$icrel[i] <- tag$icrel[1]
      ##   res$icrec[i] <- tag$icrec[1] }
      if ("use" %in% names(tag)) res$use[i] <- tag$use[1]

      next
    }

    if (nrow(tag) != 2L) {
      warning(sprintf("Element %d must have either 2 rows (long) or 1 row (already wide). Removing this tag", i))
      next
    }

    req <- c("t","x","y")
    miss <- setdiff(req, names(tag))
    if (length(miss) > 0L) {
      stop(sprintf("Element %d is missing columns: %s", i, paste(miss, collapse = ", ")))
    }

    res$id[i] <- get_id(tag, i)

    ## core coordinates
    res$t0[i] <- as.numeric(tag$t[1])
    res$t1[i] <- as.numeric(tag$t[2])
    res$x0[i] <- as.numeric(tag$x[1])
    res$x1[i] <- as.numeric(tag$x[2])
    res$y0[i] <- as.numeric(tag$y[1])
    res$y1[i] <- as.numeric(tag$y[2])

    ## strata: prefer per-row it/ic; otherwise per-tag itrel/itrec columns if
    ## present
    ## if ("it" %in% names(tag)) {
    ##   res$itrel[i] <- tag$it[1]
    ##   res$itrec[i] <- tag$it[2]
    ## } else if (has_all(tag, c("itrel","itrec"))) {
    ##   res$itrel[i] <- tag$itrel[1]
    ##   res$itrec[i] <- tag$itrec[1]
    ## }

    ## if ("ic" %in% names(tag)) {
    ##   res$icrel[i] <- tag$ic[1]
    ##   res$icrec[i] <- tag$ic[2]
    ## } else if (has_all(tag, c("icrel","icrec"))) {
    ##   res$icrel[i] <- tag$icrel[1]
    ##   res$icrec[i] <- tag$icrec[1]
    ## }

    ## use: if per-row, take first (warn if inconsistent)
    ## if ("use" %in% names(tag)) {
    ##   if (!identical(tag$use[1], tag$use[2])) {
    ##     warning(sprintf("Element %d has different 'use' values in the two rows; using the first.", i))
    ##   }
    ##   res$use[i] <- tag$use[1]
    ## }

  }

  res
}



##' Summarise tagging data
##'
##' @description
##' `summarise_tags()` prints a compact summary of tagging data stored in an
##' object of class `admove_tags` or in a higher-level \emph{admove} object that
##' contains tagging data.
##'
##' The summary includes the total number of tags, the number of tags by tag
##' type, and simple averages such as the number of observations per tag, tag
##' duration, and average step sizes in time and space.
##'
##' @param object An object of class `admove_tags`, or an object containing tagging
##'   data such as `admove_data`, `admove_sim`, or `admove`.
##' @param ... Additional arguments
##'
##' @return
##' Invisibly returns the corresponding `admove_tags` object.
##'
##' @examples
##' summarise_tags(skjepo$sim$dat$tags)
##'
##' @name summarise_tags
##' @export
summarise_tags <- function(object, ...) {
  x <- object

  if(inherits(x, "admove_sim")) {
    tags <- x$tags
  } else if(inherits(x, "admove_data")) {
    tags <- x$tags
  } else if(inherits(x, "admove")) {
    tags <- x$dat$tags
  } else if(inherits(x, "admove_tags")) {
    tags <- x
  } else stop("Please provide an object of class 'admove_tags' or an object containing an such an object (e.g. admove_data, admove_sim, admove).")

  ## rans <- sprintf("%.2f", range(tags, na.rm = TRUE))
  if (inherits(tags, "list")) {
    tags <- do.call(rbind, tags)
  }
  dims <- dim(tags)

  tags_split <- split(tags, tags$tag_type)
  tags_split2 <- lapply(tags_split, function(x) split(x, x$id))
  n_tags <- lapply(tags_split2, function(x) ifelse(length(x) == 0, NA, length(x)))
  av_nt <- lapply(tags_split2, function(x) ifelse(length(x) == 0, NA, mean(sapply(x, nrow), na.rm = TRUE)))
  av_t <- lapply(tags_split2, function(x) ifelse(length(x) == 0, NA, mean(sapply(x, function(x) x[nrow(x),"t"] - x[1,"t"]), na.rm = TRUE)))
  av_tdiff <- lapply(tags_split2, function(x) ifelse(length(x) == 0, NA, mean(unlist(sapply(x, function(x) diff(x[,"t"]))), na.rm = TRUE)))
  av_xdiff <- lapply(tags_split2, function(x) ifelse(length(x) == 0, NA, mean(unlist(sapply(x, function(x) diff(x[,"x"]))), na.rm = TRUE)))
  av_ydiff <- lapply(tags_split2, function(x) ifelse(length(x) == 0, NA, mean(unlist(sapply(x, function(x) diff(x[,"y"]))), na.rm = TRUE)))

  spinfo <- sref(tags)
  tinfo <- tref(tags)

  labw <- 15

  cat("<admove_tags>\n")
  cat(sprintf(paste0("  %-", labw, "s %s\n"), "tags total:",
              length(unique(paste0(tags$tag_type,":",tags$id)))))

  for (i in 1:3) {
    if (is.na(n_tags[[i]]) || n_tags[[i]] == 0) next

    cat("  ---------------------------------\n")
    cat(" ", paste0(c("data-storage","mark-resight","mark-recapture")[i]," tags\n"))
    cat(sprintf(paste0("  %-", labw, "s %s\n"), "n:", n_tags[[i]]))
    cat("  average over ids:\n")
    cat(sprintf(paste0("  %-", labw, "s %s\n"), "n obs:",
                sprintf("%.2f", av_nt[i])))
    cat(sprintf(paste0("  %-", labw, "s %s\n"), "duration:",
                sprintf("%.2f", av_t[i])))
    cat(sprintf(paste0("  %-", labw, "s %s\n"), "time step:",
                sprintf("%.2f", av_tdiff[i])))
    cat(sprintf(paste0("  %-", labw, "s %s\n"), "x step:",
                sprintf("%.2f", av_xdiff[i])))
    cat(sprintf(paste0("  %-", labw, "s %s\n"), "y step:",
                sprintf("%.2f", av_ydiff[i])))
  }
  cat("  ---------------------------------\n")


  ## cat(sprintf(paste0("  %-", labw, "s %s\n"), "NAs:",
  ##             n_na))

  cat(.format_sref_short(spinfo), sep = "\n")
  cat(.format_tref_short(tinfo), sep = "\n")
  cat("\n")


  invisible(tags)
}


##' Extract mark-recapture tags
##'
##' @description
##' `get_ctags()` extracts conventional mark-recapture tags from an object that
##' directly contains tagging data or from an object with a `tags` element.
##'
##' @param x An `admove_tags` object, a data frame with a `tag_type` column, or
##'   an object containing a `tags` element.
##'
##' @return
##' A subset of the input tagging data containing only mark-recapture tags
##' (`tag_type == "c"`).
##'
##' @export
get_ctags <- function(x) {
  if (any(colnames(x) == "tag_type")){
    ind <- which(x$tag_type == "c")

    x[ind,]
  } else if (any(names(x) == "tags")){
    ind <- which(x$tags$tag_type == "c")

    x$tags[ind,]
  } else stop("No tags found!")
}


##' Extract data-storage tags
##'
##' @description
##' `get_dtags()` extracts data-storage (archival) tags from an object that
##' directly contains tagging data or from an object with a `tags` element.
##'
##' @param x An `admove_tags` object, a data frame with a `tag_type` column, or
##'   an object containing a `tags` element.
##'
##' @return
##' A subset of the input tagging data containing only data-storage tags
##' (`tag_type == "d"`).
##'
##' @export
get_dtags <- function(x) {
  if (any(colnames(x) == "tag_type")){
    ind <- which(x$tag_type == "d")

    x[ind,]
  } else if (any(names(x) == "tags")){
    ind <- which(x$tags$tag_type == "d")

    x$tags[ind,]
  } else stop("No tags found!")
}


##' Extract mark-resight tags
##'
##' @description
##' `get_stags()` extracts mark-resight tags from an object that directly
##' contains tagging data or from an object with a `tags` element.
##'
##' @param x An `admove_tags` object, a data frame with a `tag_type` column, or
##'   an object containing a `tags` element.
##'
##' @return
##' A subset of the input tagging data containing only mark-resight tags
##' (`tag_type == "s"`).
##'
##' @export
get_stags <- function(x) {
  if (any(colnames(x) == "tag_type")){
    ind <- which(x$tag_type == "s")

    x[ind,]
  } else if (any(names(x) == "tags")){
    ind <- which(x$tags$tag_type == "s")

    x$tags[ind,]
  } else stop("No tags found!")
}




##' Group tags by release events
##'
##' @description
##' `use_release_events()` groups tags that were released close in space and time
##' into common release events. This is mainly intended for conventional
##' mark-recapture tags, where multiple tags may share the same approximate
##' release location and release time.
##'
##' The grouping is defined by the supplied spatial grid and time vector.
##'
##' @param x Tagging data as an object of class `admove_tags` or an
##'   `admove_data` object containing tags.
##' @param grid Spatial grid used to aggregate release locations into common
##'   release events.
##' @param time_cont Numeric time vector used to aggregate release times into
##'   common release events.
##' @param tag_types Character vector giving the tag types to group by release
##'   event. Default: `"c"`.
##'
##' @return
##' An object like the input tagging data, but with selected tags grouped by
##' common release events.
##'
##' @details
##' For each selected tag type, the first observation of each tag is treated as
##' the release event. Releases falling into the same spatial grid cell and the
##' same time interval are grouped together and assigned a common release event.
##' The remaining observations are then combined under the new grouped tag id.
##'
##' @examples
##' ## use_release_events(tags, grid = dat$grid, time_cont = dat$time_cont)
##'
##' @export
use_release_events <- function(x, grid, time_cont,
                               tag_types = "c") {
  xin <- x

  if (any(tag_types %in% c("d","s"))) warning("release events are only tested for mark-recapture tags (tag_types = 'c'). I hope you know what you do.")

  if (inherits(x, "admove_tags")) {
    tags <- x
  } else if (inherits(x, "admove_data")) {
    tags <- dat$tags
  } else stop("x must contain tags. Did you run check_tags()?")

  tags_by_type <- split(tags, tags$tag_type)

  if (is.null(grid)) stop("grid cannot be NULL")
  celltable <- grid$celltable
  xygrid <- grid$xygrid
  xgr <- grid$xgr
  ygr <- grid$ygr
  if (is.null(time_cont)) stop("time_cont cannot be NULL")

  tags_new <- vector("list", length(tags_by_type))
  for (i in 1:length(tags_by_type)) {

      tagi <- tags_by_type[[i]]

    if (names(tags_by_type)[i] %in% tag_types) {
      tagi_split <- split(tagi, tagi$id)

      tagi_rel <- do.call(rbind, lapply(tagi_split, function(x) x[1,c("t","x","y")]))

      idx_space <- celltable[cbind(cut(tagi_rel$x, xgr),
                                   cut(tagi_rel$y, ygr))]

      idx_time <- as.integer(cut(tagi_rel$t, time_cont, include.lowest = TRUE))

      rel_all <- data.frame(time_cont[idx_time],
                            xygrid[idx_space,])
      colnames(rel_all) <- c("t","x","y")

      ## Unique release events
      release_events <- rel_all[!duplicated(rel_all),]
      rownames(release_events) <- NULL

      ## Index to match tag to release event
      idx_all <- apply(rel_all, 1, paste, collapse = ":")
      idx_uni <- apply(release_events, 1, paste, collapse = ":")
      idx <- match(idx_all, idx_uni)

      n_rel <- nrow(release_events)
      tagi_new <- vector("list", n_rel)
      for (j in 1:n_rel) {
        idx_j <- which(idx == j)
        tmp <- do.call(rbind, lapply(tagi_split[idx_j], function(x) x[-1,]))
        tmp_rel <- tmp[1,]
        tmp_rel[,c("t","x","y")] <- release_events[j,c("t","x","y")]
        tmp_comb <- rbind(tmp_rel, tmp)
        rownames(tmp_comb) <- NULL
        tmp_comb$id <- j
        tagi_new[[j]] <- tmp_comb
      }

      tagi_new <- do.call(rbind, tagi_new)

    } else {
      tagi_new <- tagi
    }

    tags_new[[i]] <- tagi_new
  }

  tags_out <- do.call(rbind, tags_new)

  return(tags_out)
}


##' Plot tagging data
##'
##' @description
##' `plot_tags()` plots tagging data in space, showing release locations,
##' intermediate observations where available, and final recovery or resight
##' locations.
##'
##' Depending on the tag type, tags are drawn either as trajectories connecting
##' successive observations or as straight lines between release and recovery.
##'
##' @param x An object of class `admove_tags`, or an object containing tagging
##'   data such as `admove_data`, `admove_sim`, or `admove`.
##' @param main Main title of the plot. Default: `"Tags"`.
##' @param plot_land Logical; if `TRUE`, add land to the plot. Default: `FALSE`.
##' @param auto_layout Logical; if `TRUE`, plotting parameters are set
##'   automatically and restored afterwards. Default: `TRUE`.
##' @param xlim Optional x-axis limits.
##' @param ylim Optional y-axis limits.
##' @param add Logical; if `TRUE`, add to an existing plot. Default: `FALSE`.
##' @param xlab Label for the x-axis. Default: `"x"`.
##' @param ylab Label for the y-axis. Default: `"y"`.
##' @param leg_pos Position of the legend. Default: `"topright"`.
##' @param labels Logical; if `TRUE`, label observations by time instead of
##'   plotting intermediate points. Default: `FALSE`.
##' @param bg Optional background colour for the plot. Default: `NULL`.
##' @param by_tag_type Logical; if `TRUE`, create separate panels by tag type.
##'   Default: `TRUE`.
##' @param by_tag Logical; if `TRUE`, create separate panels for individual tags.
##'   Default: `FALSE`.
##' @param col Character vector of length 1 to 3 giving colours for tag paths,
##'   release positions, and recovery or final observation positions.
##' @param ... Additional graphical arguments passed to [plot()].
##'
##' @return
##' No return value. Called for its side effect of producing a plot.
##'
##' @details
##' Data-storage and mark-resight tags with intermediate observations are plotted
##' as trajectories through space. Conventional mark-recapture tags are plotted
##' as straight-line segments between release and recovery positions.
##'
##' When `by_tag_type = TRUE`, separate panels are created for each tag type
##' present in the data. When `by_tag = TRUE`, each tag is shown in a separate
##' panel.
##'
##' @examples
##' plot_tags(skjepo$sim$tags)
##'
##' @name plot_tags
##' @export
plot_tags <- function(x,
                      main = "Tags",
                      plot_land = FALSE,
                      auto_layout = TRUE,
                      xlim = NULL,
                      ylim = NULL,
                      add = FALSE,
                      xlab = "x",
                      ylab = "y",
                      leg_pos = "topright",
                      labels = FALSE,
                      bg = NULL,
                      by_tag_type = TRUE,
                      by_tag = FALSE,
                      col = c(adjustcolor("grey60",0.3), .admove_cols(2)),
                      pch = c(1,0,16),
                      cex = 0.8,
                      ...) {

  pchin <- pch
  pch0 <- c(1,0,16)
  if (length(pchin) != 3) {
    pch <- pch0
    pch[1:length(pchin)] <- pchin
  }

  if (inherits(x, "admove_data")) {
    tags <- x$tags
    sref <- sref(x)
    tref <- tref(x)
  } else   if (inherits(x, "admove_sim")) {
    tags <- x$tags
    sref <- sref(x$dat)
    tref <- tref(x$dat)
  } else if(inherits(x, "admove")) {
    tags <- x$dat$tags
    sref <- sref(x$dat)
    tref <- tref(x$dat)
  } else {
    tags <- x
    .check_class(tags, "admove_tags")
    sref <- sref(tags)
    tref <- tref(tags)
  }

  if (inherits(tags, "data.frame")) {
    tags <- split(tags, tags$id)
  }

  if (is.null(xlim)) {
    xlims <- range(sapply(tags, function(x) range(x[,2], na.rm = TRUE)))
  } else {
    xlims <- xlim
  }
  if (is.null(ylim)) {
    ylims <- range(sapply(tags, function(x) range(x[,3], na.rm = TRUE)))
  } else {
    ylims <- ylim
  }

  if (auto_layout) {
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar))
  }

  cols0 <- c(adjustcolor("grey60",0.3), .admove_cols(2))

  cols_use <- cols0
  if (!missing(col) && !is.null(col)) {
    if (!is.character(col)) stop("'col' must be a character vector (e.g., 'red' or c('red','blue')).")
    if (length(col) > 3) stop("'col' must have length 1, 2, or 3.")

    idx <- seq_len(length(col))
    cols_use[idx] <- col
  }
  cols <- cols_use

  plot_one <- function(tags) {

    if(!add){
      if(!is.null(bg)){
        par(bg = bg)
      }
      plot(0,0, ty = "n", main = main,
           xlim = xlims,
           ylim = ylims,
           xlab = xlab, ylab = ylab,
           xaxt = xaxt,
           yaxt = yaxt,
           asp = 1,
           ...)
    }
    if (plot_land) {
      plot_land(sref)
    }
    start_pos <- t(sapply(tags, function(x) x[1, c(2,3)]))
    end_pos <- t(sapply(tags, function(x) x[nrow(x), c(2,3)]))
    points(start_pos[,1], start_pos[,2],
           col = cols[2], pch = pch[1])
    points(end_pos[,1], end_pos[,2],
           col = cols[3], pch = pch[2])
    tmp <- sapply(tags, function(x) .get_tag_type_integer(x$tag_type)[1])
    idx1 <- which(tmp %in% c(1,2,4))
    idx2 <- which(tmp %in% c(3))
    ## tmp <- sapply(tags, nrow) > 2
    ## idx1 <- as.integer(which(tmp))
    ## idx2 <- as.integer(which(!tmp))
    for (i in idx1) {
      lines(tags[[i]][,2], tags[[i]][,3],
            col = cols[1], ty = "b", pch = NA)
      if (labels) {
        text(tags[[i]][,2], tags[[i]][,3],
             labels = sprintf("%.2f", tags[[i]][,1]),
             col = cols[1], pch = pch[3], cex = cex)
      } else {
        points(tags[[i]][-c(1,nrow(tags[[i]])),2],
               tags[[i]][-c(1,nrow(tags[[i]])),3],
               col = cols[1], pch = pch[3], cex = cex)
      }
    }

    if (length(idx2)) {
      segments(as.numeric(start_pos[,1]), as.numeric(start_pos[idx2,2]),
               as.numeric(end_pos[idx2,1]), as.numeric(end_pos[idx2,2]),
               col = cols[1])
    }

  }

  if (by_tag) {

    tag_types <- do.call(rbind,tags)$tag_type
    n <- length(tags)

    mfrow <- n2mfrow(n, asp = 2)
    if(auto_layout){
      par(mfrow = mfrow,
          mar = c(0.1,0.1,0.1,0.1),
          oma = c(4,4,1,1))
    }
    main <- ""
  } else if (by_tag_type) {
    tag_types <- unique(do.call(rbind,tags)$tag_type)
    n <- length(tag_types)
    mfrow <- n2mfrow(n, asp = 2)
  }else {
    n <- 1

    mfrow <- c(1,1)
    if(auto_layout){
      par(mfrow = mfrow)
    }
    main <- main
  }

  for (i in 1:n) {

    xaxt <- ifelse(i %in% (prod(mfrow) - mfrow[2]+1):prod(mfrow), "s", "n")
    yaxt <- ifelse(i %in% seq(1, prod(mfrow), mfrow[2]), "s", "n")

    if (by_tag) {
      tags2 <- tags[i]
    } else if (by_tag_type) {
      tags2 <- tags[which(sapply(tags, function(x) all(x$tag_type == tag_types[i])))]
    } else {
      tags2 <- tags
    }

    plot_one(tags2)

    if (i == 1 || by_tag_type) {
      if (by_tag_type || by_tag) {
      labo <- list(c("Deployment", "Intermediate obs.", "Recovery"),
                  c("Release", "Resights", "Final resight"),
                  c("Release", "Recovery"),
                  c("Release", "Detection"))[[as.integer(tag_types[i])]]
      pcho <- list(c(1,16,0),
                   c(1,16,0),
                   c(1,0),
                   c(1,16))[[as.integer(tag_types[i])]]
      colo <- list(cols[c(2,1,3)],
                   cols[c(2,1,3)],
                   cols[c(2,3)],
                   cols[c(2,1)])[[as.integer(tag_types[i])]]
      } else if (any(sapply(tags2, nrow) > 2)){
        labo <- c("Release", "Intermediate obs.", "Final obs.")
        pcho <- c(1,16,0)
        colo <- cols[c(2,1,3)]
      } else {
        labo <- c("Release", "Recovery")
        pcho <- c(1,0)
        colo <- cols[c(2,3)]
      }
      legend(leg_pos,
             legend = labo,
             pch = pcho,
             col = colo,
             bg = "white")
    }
    box(lwd = 1.5)
    if (by_tag) {
      legend("topleft",
             legend = round(tags2[[1]][1,1],3),
             cex = 0.8,
             pch = NA,
             bg = "white")
    }
  }

  if (by_tag) {
    mtext(xlab, 1, 2, outer = TRUE)
    mtext(ylab, 2, 2, outer = TRUE)
  }
}


##' @rdname prep_tags
##' @export
prep_dtags <- function(x,
                       names = NULL,
                       date_decimal = FALSE,
                       date_format = NULL,
                       date_origin = NULL,
                       tz = "UTC",
                       sref = NULL,
                       tref = NULL,
                       transform_sref = FALSE,
                       shift_tref = FALSE,
                       verbose = TRUE) {
  prep_tags(
    x = x,
    tag_type = "d",
    names = names,
    date_decimal = date_decimal,
    date_format = date_format,
    date_origin = date_origin,
    tz = tz,
    sref = sref,
    tref = tref,
    transform_sref = transform_sref,
    shift_tref = shift_tref,
    verbose = verbose
  )
}

##' @rdname prep_tags
##' @export
prep_stags <- function(x,
                       names = NULL,
                       date_decimal = FALSE,
                       date_format = NULL,
                       date_origin = NULL,
                       tz = "UTC",
                       sref = NULL,
                       tref = NULL,
                       transform_sref = FALSE,
                       shift_tref = FALSE,
                       verbose = TRUE) {
  prep_tags(
    x = x,
    tag_type = "s",
    names = names,
    date_decimal = date_decimal,
    date_format = date_format,
    date_origin = date_origin,
    tz = tz,
    sref = sref,
    tref = tref,
    transform_sref = transform_sref,
    shift_tref = shift_tref,
    verbose = verbose
  )
}

##' @rdname prep_tags
##' @export
prep_ctags <- function(x,
                       names = NULL,
                       date_decimal = FALSE,
                       date_format = NULL,
                       date_origin = NULL,
                       tz = "UTC",
                       sref = NULL,
                       tref = NULL,
                       transform_sref = FALSE,
                       shift_tref = FALSE,
                       verbose = TRUE) {
  prep_tags(
    x = x,
    tag_type = "c",
    names = names,
    date_decimal = date_decimal,
    date_format = date_format,
    date_origin = date_origin,
    tz = tz,
    sref = sref,
    tref = tref,
    transform_sref = transform_sref,
    shift_tref = shift_tref,
    verbose = verbose
  )
}




## Internal functions -----------------------------------------------------------------

.get_tag_type <- function(type) {

  if (is.numeric(type)) {
    type <- c("d","s","c","a")[type]
  }

  ## normalize
  x <- tolower(trimws(as.character(type)))

  ## map common aliases -> canonical levels
  map <- list(
    d = c("d", "dtag", "dtags",
          "data-storage", "datastorage",
          "data-logging", "datalogging",
          "archival", "archival-tag", "archivaltag"),
    s = c("s", "stag", "stags",
          "mark-resight", "mark-resighting", "mark-resighting-tag",
          "mark-resight-tag", "markresight", "markresighting"),
    c = c("c", "ctag", "ctags",
          "mark-recapture", "mark-recapturing",
          "conventional", "conventional-tag", "conventionaltag"),
    a = c("a", "atag", "atags",
          "acoustic")
  )

  canonical <- rep(NA_character_, length(x))
  for (k in names(map)) {
    canonical[x %in% map[[k]]] <- k
  }

  ## allow already-canonical single letters
  canonical[is.na(canonical) & x %in% c("d","s","c")] <- x[is.na(canonical) & x %in% c("d","s","c")]

  ## error on unknown
  if (anyNA(canonical)) {
    bad <- unique(x[is.na(canonical)])
    stop("Unknown tag type: ", paste(shQuote(bad), collapse = ", "),
         ". Allowed types include: d/s/c, dtags/stags/ctags, archival/data-storage/data-logging, ",
         "mark-resight/mark-resighting, mark-recapture/conventional.",
         call. = FALSE)
  }

  factor(canonical, levels = c("d","s","c","a"))
}


.get_random_id <- function(n = 6) {
  paste(sample(c(letters,0:9), n, replace = TRUE), collapse="")
}



## s3 methods -------------------------------------------------------------------------

##' @rdname plot_tags
##' @export
plot.admove_tags <- function(x, ...) {
  plot_tags(x, ...)
}

##' @rdname combine_tags
##' @export
c.admove_tags <- function(..., recursive = FALSE) {
  combine_tags(..., recursive = recursive)
}

##' @method summary admove_tags
##' @rdname summarise_tags
##' @export
summary.admove_tags <- function(object, ...) {
  summarise_tags(object, ...)
}

##' @rdname print-admove
##' @method print admove_tags
##' @export
print.admove_tags <- function(x, ...) {
  tmp <- x
  attributes(tmp) <- NULL
  NextMethod("print", tmp, ...)
}
