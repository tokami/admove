
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
##' @param names Named character vector (or list) giving the column names in
##'   `x`. For long format, provide at least `c(t = "...", x = "...", y = "...",
##'   id = "...")`, optionally with `event` and `prob` (see
##'   \emph{Ambiguous positions} below). For wide format, provide
##'   `c(t0 = "...", t1 = "...", x0 = "...", y0 = "...", x1 = "...", y1 = "...")`,
##'   optionally with `p1`. The order does not matter, but the vector must be
##'   named. Given as a list, the wide entries `t1`, `x1`, `y1` and `p1` may
##'   each hold several column names -- one per candidate recapture location.
##' @param candidates Optional integer or character vector of candidate indices.
##'   When supplied, the wide entries `t1`, `x1`, `y1` and `p1` of `names` are
##'   read as column \emph{stems} and these indices are appended, so
##'   `names = c(t1 = "date", x1 = "lon", y1 = "lat", p1 = "per")` together with
##'   `candidates = 1:9` maps `date1`, `lon1`, `lat1`, `per1`, `date2`, ... That
##'   is the layout recapture-uncertainty tables usually come in. Default:
##'   `NULL`, which leaves the entries as literal column names.
##' @param date_decimal Logical; if `TRUE`, interpret the time variable as a
##'   decimal year and convert it to model time. Default: `FALSE`.
##' @param date_format Optional format string (see [strptime()]) to parse
##'   character dates or date-times, e.g. `"%m/%d/%Y %H:%M"`. The time of day is
##'   kept. Dates that cannot be parsed are reported in a warning. Default:
##'   `NULL`.
##' @param date_origin Optional origin when times are stored numerically as
##'   (possibly fractional) days since `date_origin`, e.g. `"1899-12-30"` for
##'   spreadsheet dates. Default: `NULL`.
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
##' @section Ambiguous positions:
##'
##' The final position of a tag is sometimes known only up to a set of
##' candidates. The usual case is a mark-recapture tag whose recapturing vessel
##' is known but whose individual set is not: any of that vessel's fishing sets
##' could be the recapture location, each with a probability derived from its
##' effort. Two optional columns express this:
##'
##' \itemize{
##'   \item `event` -- rows of one tag that share an `event` are mutually
##'     exclusive \strong{alternatives} for a single observation, not successive
##'     observations.
##'   \item `prob` -- the probability of each alternative. These must already be
##'     probabilities and sum to 1 within an event; convert an effort measure
##'     first, e.g. `prob = hours / sum(hours)`.
##' }
##'
##' Candidates may sit at different times as well as different positions. Only
##' the \strong{final} observation of a tag may be ambiguous, and the release
##' never can: that restriction is what keeps the likelihood an exact mixture
##' \eqn{\log \sum_k p_k \, f(x_k, t_k \mid \mathrm{release})} in both the
##' Kalman-filter and the CTMC engine, with no Gaussian-mixture posterior to
##' approximate. Both rules are enforced by [check_tags()].
##'
##' No state update is performed at an ambiguous observation, so
##' `conf$do_update` has no effect there. Nothing is propagated past the final
##' observation, which is why skipping the update costs nothing and every
##' candidate is evaluated against the same release-conditioned prediction.
##'
##' Omitting the columns means "no ambiguity" and reproduces the behaviour of
##' earlier versions exactly. See [add_candidates()] for attaching candidates
##' that arrive as a separate table.
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
##' ## mark-recapture tags whose recapture position is ambiguous, given as
##' ## repeated columns date1, lon1, lat1, per1, date2, ... in a wide table
##' ## ctags <- prep_ctags(
##' ##   unc,
##' ##   names = c(id = "fish_id", t0 = "release_date",
##' ##             x0 = "release_lon", y0 = "release_lat",
##' ##             t1 = "date", x1 = "lon", y1 = "lat", p1 = "per"),
##' ##   candidates = 1:9, date_origin = "1899-12-30")
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
                      candidates = NULL,
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

  ## 'names' is a named list internally so that the wide-format entries t1, x1,
  ## y1 and p1 can each carry several column names -- one per candidate
  ## recapture location. A named character vector (the classic input) is just
  ## the case where every entry has length one.
  if (!is.null(names) && !is.list(names)) names <- as.list(names)

  ## 'candidates' is sugar for the common stem+index layout of uncertainty
  ## tables (date1, lat1, lon1, per1, date2, ...): the t1/x1/y1/p1 entries are
  ## read as stems and the indices are appended. With candidates = NULL they
  ## stay literal column names, which is the previous behaviour.
  if (!is.null(candidates)) {
    if (is.null(names)) stop("'candidates' given but 'names' is missing. Provide the column stems, e.g. names = c(t1 = 'date', x1 = 'lon', y1 = 'lat', p1 = 'per') together with candidates = 1:9. See ?prep_tags.")
    for (nm in c("t1","x1","y1","p1")) {
      if (!is.null(names[[nm]])) {
        if (length(names[[nm]]) != 1) stop("With 'candidates' the '", nm, "' entry of 'names' must be a single column stem, not ", length(names[[nm]]), " names. Either drop 'candidates' and list the columns explicitly, or give one stem.")
        names[[nm]] <- paste0(names[[nm]], candidates)
      }
    }
  }

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

  ## Optional canonical entries. For the wide format 'p1' carries the
  ## probability of each candidate recapture location; for the long format the
  ## same information is given directly as 'event' (which rows are alternatives
  ## for one another) and 'prob'.
  opt <- if (flag_wide) "p1" else c("event", "prob")
  opt <- opt[opt %in% base::names(names)]

  ## Rename the requested columns to their canonical names. An entry of 'names'
  ## may hold several column names (t1/x1/y1/p1 in the wide format, one per
  ## candidate recapture location); those become t1.1, t1.2, ... so that
  ## ctags_wide_2_long() can pick the candidates up by position.
  x_in <- x
  sel <- character(0)
  col_names <- character(0)
  for (nm in c(req, opt)) {
    cols_nm <- as.character(unlist(names[[nm]]))
    if (length(cols_nm) == 0) next()
    sel <- c(sel, cols_nm)
    col_names <- c(col_names,
                   if (length(cols_nm) == 1) nm else paste0(nm, ".", seq_along(cols_nm)))
  }

  dups <- unique(sel[duplicated(sel)])
  if (length(dups) > 0) stop("The same input column is mapped more than once in 'names': ", paste(dups, collapse = ", "), ".")

  miss <- setdiff(sel, colnames(x_in))
  if (length(miss) > 0) stop("Column(s) named in 'names' but not found in the input data: ", paste(miss, collapse = ", "), ".")

  x <- x_in[, sel, drop = FALSE]

  idx <- which(!(colnames(x_in) %in% sel))
  if (length(idx) > 0) {
    x <- cbind(x, x_in[, idx, drop = FALSE])
    col_names <- c(col_names, colnames(x_in)[idx])
  }
  colnames(x) <- col_names

  if (flag_id) colnames(x)[colnames(x) == as.character(names[["id"]])] <- "id"

  if (flag_wide) {
    x <- ctags_wide_2_long(x)
  }

  ## Convert dates
  if (!is.null(date_origin) ||
        !is.null(date_format) ||
         isTRUE(date_decimal)) {

    dati <- .parse_dates(x$t, date_format = date_format,
                         date_origin = date_origin,
                         date_decimal = date_decimal, tz = tz)

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
##'   \item removal of observations outside the temporal domain of `dat`,
##'   \item chronological reordering of archival and mark-resight tags whose
##'     observations are not sorted in time,
##'   \item removal of observations that repeat the time of the previous
##'     observation of the same tag (zero time step),
##'   \item removal of mark-recapture tags whose recapture time is at or before
##'     the release time, and
##'   \item optional removal of tags with fewer than two observations.
##' }
##'
##' It also canonicalises the optional `event` and `prob` columns that express an
##' ambiguous final position (see [prep_tags()]), creating them when absent
##' (one event per row, probability 1). Candidate positions of one event are
##' allowed to share a time -- they are alternatives, never successive steps --
##' so the time checks above compare \emph{events} rather than rows. An
##' ambiguous release, an ambiguous observation that is not the last one, or
##' probabilities that do not sum to 1 within an event are errors. When the
##' checks remove some candidates of an event, the probabilities of the
##' survivors are rescaled to sum to 1 again.
##'
##' The time checks matter because the likelihood builds its time axis from the
##' sorted observation times while reading the positions in the order they are
##' stored, so unsorted times would pair times with the wrong positions. Tags
##' with a non-positive time step additionally break the internal time grid.
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


  ## ids shared between tag types
  ## Everything downstream groups observations with split(tags, tags$id), so an
  ## id reused by two different tag types silently merges two physical tags into
  ## one: their positions are interleaved into a single track, the merged tag
  ## takes the tag_type of whichever row comes first, and the internal time step
  ## collapses to the finer of the two -- which can inflate the number of
  ## integration steps by orders of magnitude. This happens easily because
  ## prep_dtags()/prep_stags() fall back to numbering tags 1, 2, 3, ... when no
  ## id column is supplied, which readily collides with numeric ids carried by
  ## mark-recapture tags. Make such ids unique per tag type and say so.
  tags <- .disambiguate_ids(tags, verbose)

  ## observation events: canonicalise / validate the optional event + prob
  ## columns that express an ambiguous (multi-candidate) final position
  tags <- .check_events(tags, verbose)

  ## duplicate ids within the mark-recapture tags (the one within-type case that
  ## is detectable: a "c" tag must have exactly two observation events)
  .check_ctag_rows(tags, verbose)


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

    ## Tag positions that fall on an NA (masked / removed) grid cell have no
    ## valid cell index and must be dropped before fitting. Entries outside the
    ## grid range are already handled above, so a remaining NA 'ic' means the
    ## position sits on an NA cell.
    ind <- which(is.na(tags$ic))
    if (length(ind) > 0) {
      na_ids <- unique(tags$id[ind])
      tags <- tags[-ind,]
      if (verbose) {
        message(length(ind), " entr", if (length(ind) == 1) "y" else "ies",
                " removed because the tag position falls on an NA grid cell (tag id",
                if (length(na_ids) == 1) "" else "s", ": ",
                .format_ids(na_ids), ").")
      }
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

  ## time ordering within tags
  ## Row order carries information: build_time() sorts the times internally, but
  ## the rows (positions) keep their original order, so unsorted times silently
  ## mismatch times and positions. Non-positive time steps additionally make
  ## dt_min <= 0, which errors in build_time().
  ## Comparisons are between EVENTS, not rows: the candidate positions of one
  ## ambiguous observation share an event and may tie in time or run in any
  ## order among themselves, which is legitimate rather than unsorted input.
  ## .check_events() keys the events by order of appearance, so with no
  ## ambiguity every event holds exactly one row and the checks below reduce to
  ## the per-row ones they replace.
  row_list <- split(seq_len(nrow(tags)), tags$id)
  ord <- seq_len(nrow(tags))
  bad_ids <- NULL
  resorted_ids <- NULL
  for (nm in names(row_list)) {
    rows <- row_list[[nm]]
    if (length(rows) < 2) next
    ti <- tags$t[rows]
    ei <- tags$event[rows]
    ue <- unique(ei)
    t_start <- as.numeric(tapply(ti, ei, min))[match(ue, sort(unique(ei)))]
    t_end <- as.numeric(tapply(ti, ei, max))[match(ue, sort(unique(ei)))]
    if (tags$tag_type[rows][1] == "c") {
      ## for mark-recapture tags the row order is meaningful (release first,
      ## recapture second), so a recapture time at or before the release time is
      ## a data error rather than unsorted input
      if (length(ue) < 2 || any(t_start[-1] <= t_end[-length(t_end)])) {
        bad_ids <- c(bad_ids, nm)
      }
    } else {
      ## sort by event time, keeping the candidate rows of one event together.
      ## With no ambiguity every event is a single row and this is order(ti).
      oo <- order(t_start[match(ei, ue)], ti)
      if (!identical(oo, seq_along(rows))) {
        ord[rows] <- rows[oo]
        resorted_ids <- c(resorted_ids, nm)
      }
    }
  }

  if (!is.null(resorted_ids)) {
    tags <- tags[ord, , drop = FALSE]
    rownames(tags) <- NULL
    ## rows moved, so re-key the events by their new order of appearance
    tags$event <- stats::ave(tags$event, tags$id,
                             FUN = function(z) match(z, unique(z)))
    if (verbose) message(length(resorted_ids), " tag", if (length(resorted_ids) == 1) "" else "s",
                         " reordered because the observations were not in chronological order (id",
                         if (length(resorted_ids) == 1) "" else "s", ": ",
                         .format_ids(resorted_ids), ").")
  }

  ## duplicated times (dt = 0) within a tag; positions are contiguous per id
  ## Only duplicates ACROSS events are a zero time step. Two candidate positions
  ## of the same ambiguous observation are allowed to share a time -- they are
  ## alternatives, never consecutive steps.
  drop_rows <- integer(0)
  for (nm in names(row_list)) {
    rows <- row_list[[nm]]
    if (length(rows) < 2 || nm %in% bad_ids) next
    ti <- tags$t[rows]
    ei <- tags$event[rows]
    dup <- which(diff(ti) == 0 & diff(ei) != 0)
    if (length(dup) > 0) drop_rows <- c(drop_rows, rows[dup + 1])
  }
  if (length(drop_rows) > 0) {
    if (verbose) message(length(drop_rows), " entr", if (length(drop_rows) == 1) "y" else "ies",
                         " removed because the time is identical to the previous observation of the same tag (zero time step).")
  }
  if (!is.null(bad_ids)) {
    if (verbose) message(length(bad_ids), " mark-recapture tag", if (length(bad_ids) == 1) "" else "s",
                         " removed because the recapture time is at or before the release time (id",
                         if (length(bad_ids) == 1) "" else "s", ": ",
                         .format_ids(bad_ids), ").")
    drop_rows <- c(drop_rows, which(tags$id %in% bad_ids))
  }
  if (length(drop_rows) > 0) {
    tags <- tags[-unique(drop_rows), , drop = FALSE]
    if (nrow(tags) == 0) stop("No tags passed the checks!")
  }

  ## the filters above may have removed candidate positions of an ambiguous
  ## recapture, which breaks the probability normalisation
  tags <- .renormalise_events(tags, verbose)

  ## keep only recovered tags
  tags_list <- split(tags, tags$id)
  if (remove_non_recovered_tags) {
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


## Make tag ids unique across tag types.
##
## Only ids that actually occur under more than one tag_type are rewritten, as
## `<tag_type>-<id>`; every other id is left untouched so that plot labels and
## any user-side matching keep working. Repeats until the new ids are free, in
## the unlikely event that a constructed name is already taken.
##
## Note what is deliberately NOT flagged here: an id repeated *within* one tag
## type is the normal representation, not an error -- that is exactly how the
## successive positions of a data-storage tag are grouped. Reused ids within a
## type are therefore undetectable in general; the one exception is handled by
## .check_ctag_rows() below.
.disambiguate_ids <- function(tags, verbose = TRUE) {

  if (is.null(tags) || nrow(tags) == 0) return(tags)

  ids <- as.character(tags$id)
  types <- as.character(tags$tag_type)

  ntype <- tapply(types, ids, function(z) length(unique(z)))
  shared <- names(ntype)[!is.na(ntype) & ntype > 1]
  if (length(shared) == 0) return(tags)

  rows <- which(ids %in% shared)
  new <- paste0(types[rows], "-", ids[rows])

  ## avoid colliding with an id that already exists
  taken <- setdiff(unique(ids), shared)
  pre <- ""
  while (any(new %in% taken)) {
    pre <- paste0(pre, "_")
    new <- paste0(types[rows], pre, "-", ids[rows])
  }

  tags$id <- ids
  tags$id[rows] <- new

  if (verbose) {
    message(length(shared), " tag id", if (length(shared) == 1) "" else "s",
            " used by more than one tag type (",
            .format_ids(shared),
            "). Observations of different tag types would otherwise be merged ",
            "into a single tag; the affected ids have been made unique per tag ",
            "type (e.g. \"", shared[1], "\" -> \"", sort(unique(new))[1],
            "\"). Supply distinct ids across tag types to avoid this.")
  }

  tags
}


## Canonicalise and validate the observation-event structure.
##
## Rows of one tag that share an 'event' are mutually exclusive ALTERNATIVES for
## a single observation, not successive observations: the tag was recovered at
## exactly one of them, and 'prob' says how likely each is. That is how an
## ambiguous recapture is expressed -- the recapturing vessel is known and the
## tag could have been taken at any of its fishing sets, weighted by effort.
##
## Absent columns mean "no ambiguity": every row becomes its own event with
## probability 1, which reproduces the previous behaviour exactly.
##
## Event ids are re-keyed to 1..n per tag in order of first appearance, so the
## input row order is preserved. That matters because row order already carries
## meaning for mark-recapture tags (release first, recovery second) -- keying by
## time instead would quietly repair a reversed recapture that check_tags() is
## supposed to flag. Archival and mark-resight tags are sorted by time further
## down in check_tags(), which renumbers the events again afterwards.
## Downstream code relies on the result: nll() treats max(event) as the final
## observation, and the time-ordering checks compare events rather than rows.
.check_events <- function(tags, verbose = TRUE) {

  if (is.null(tags) || nrow(tags) == 0) return(tags)

  if (!any(colnames(tags) == "event")) {
    tags$event <- stats::ave(seq_len(nrow(tags)), tags$id, FUN = seq_along)
  }
  if (!any(colnames(tags) == "prob")) tags$prob <- 1

  if (!is.numeric(tags$prob)) {
    pr <- suppressWarnings(as.numeric(as.character(tags$prob)))
    if (any(is.na(pr) & !is.na(tags$prob))) {
      stop("The 'prob' column of the tags is not numeric and could not be converted. It must hold one probability per candidate position.", call. = FALSE)
    }
    tags$prob <- pr
  }

  if (any(is.na(tags$event))) {
    stop("Missing values in the 'event' column of the tags. Every observation must belong to an event; leave the column out entirely if the tags carry no ambiguous positions.", call. = FALSE)
  }

  rows_by_id <- split(seq_len(nrow(tags)), tags$id)

  bad_first <- bad_many <- bad_late <- bad_prob <- bad_neg <- dup_pos <- NULL
  prob_sums <- NULL
  ev_new <- integer(nrow(tags))

  for (nm in names(rows_by_id)) {

    rows <- rows_by_id[[nm]]
    ev <- as.character(tags$event[rows])

    ## key the events by order of first appearance
    ev_int <- match(ev, unique(ev))
    ev_new[rows] <- ev_int

    n_ev <- max(ev_int)
    sizes <- tabulate(ev_int, nbins = n_ev)
    amb <- which(sizes > 1)

    if (length(amb) == 0) next()

    ## an ambiguous RELEASE has no anchor to propagate from
    if (1L %in% amb) bad_first <- c(bad_first, nm)
    if (length(amb) > 1) bad_many <- c(bad_many, nm)
    if (!all(amb %in% n_ev)) bad_late <- c(bad_late, nm)

    pr <- tags$prob[rows][ev_int %in% amb]
    if (any(is.na(pr)) || any(pr < 0)) {
      bad_neg <- c(bad_neg, nm)
    } else {
      for (k in amb) {
        sk <- sum(tags$prob[rows][ev_int == k])
        if (abs(sk - 1) > 1e-6) {
          bad_prob <- c(bad_prob, nm)
          prob_sums <- c(prob_sums, sk)
        }
      }
    }

    ## a repeated position within one candidate set usually means the candidate
    ## list was mis-paired upstream; worth a look but not an error
    for (k in amb) {
      kk <- rows[ev_int == k]
      if (anyDuplicated(paste(tags$x[kk], tags$y[kk]))) dup_pos <- c(dup_pos, nm)
    }
  }

  tags$event <- ev_new

  if (!is.null(bad_first)) {
    stop(length(bad_first), " tag",
         if (length(bad_first) == 1) " gives" else "s give",
         " several alternatives for the FIRST observation (id",
         if (length(bad_first) == 1) "" else "s", ": ", .format_ids(bad_first),
         "). The release is the anchor the model propagates from and must be a ",
         "single position; only the final observation may be ambiguous.",
         call. = FALSE)
  }
  if (!is.null(bad_many)) {
    stop(length(bad_many), " tag",
         if (length(bad_many) == 1) " has" else "s have",
         " more than one ambiguous observation (id",
         if (length(bad_many) == 1) "" else "s", ": ", .format_ids(bad_many),
         "). Only the last observation of a tag may have several candidate ",
         "positions, so that the likelihood stays an exact mixture in both ",
         "engines.", call. = FALSE)
  }
  if (!is.null(bad_late)) {
    stop(length(bad_late), " tag",
         if (length(bad_late) == 1) " has" else "s have",
         " an ambiguous observation that is not the last one (id",
         if (length(bad_late) == 1) "" else "s", ": ", .format_ids(bad_late),
         "). Only the final observation of a tag may have several candidate ",
         "positions.", call. = FALSE)
  }
  if (!is.null(bad_neg)) {
    stop(length(bad_neg), " tag",
         if (length(bad_neg) == 1) " has" else "s have",
         " missing or negative candidate probabilities in 'prob' (id",
         if (length(bad_neg) == 1) "" else "s", ": ", .format_ids(bad_neg), ").",
         call. = FALSE)
  }
  if (!is.null(bad_prob)) {
    stop(length(bad_prob), " tag",
         if (length(bad_prob) == 1) " has" else "s have",
         " candidate probabilities that do not sum to 1 (id",
         if (length(bad_prob) == 1) "" else "s", ": ", .format_ids(bad_prob),
         "; observed sum",
         if (length(unique(signif(prob_sums, 6))) == 1)
           paste0(" ", signif(prob_sums[1], 6))
         else
           paste0("s from ", paste(signif(range(prob_sums), 6), collapse = " to ")),
         "). 'prob' must be a probability per candidate position: convert an ",
         "effort measure first, e.g. prob = hours / sum(hours) within each ",
         "recapture event.", call. = FALSE)
  }
  if (!is.null(dup_pos) && verbose) {
    warning(length(dup_pos), " tag",
            if (length(dup_pos) == 1) " lists" else "s list",
            " the same candidate position twice within one recapture event (id",
            if (length(dup_pos) == 1) "" else "s", ": ", .format_ids(dup_pos),
            "). That usually means positions and probabilities were paired up ",
            "wrongly upstream. Please check the input data.", call. = FALSE)
  }

  tags
}


## Re-key events and restore the probability normalisation after rows have been
## dropped.
##
## check_tags() and setup_data() both remove observations (outside the grid, on
## an NA cell, outside the time range, on an NA covariate). Losing one candidate
## of an ambiguous recapture leaves the remaining probabilities summing to less
## than 1, which would silently down-weight that tag in the likelihood, so the
## survivors are re-scaled. A tag that loses its anchor -- i.e. whose first
## remaining event is the ambiguous one -- cannot be fitted and is dropped.
.renormalise_events <- function(tags, verbose = TRUE) {

  if (is.null(tags) || nrow(tags) == 0) return(tags)
  if (!any(colnames(tags) == "event")) return(tags)

  rows_by_id <- split(seq_len(nrow(tags)), tags$id)

  ev_new <- integer(nrow(tags))
  pr_new <- tags$prob
  rescaled <- NULL
  orphan <- NULL

  for (nm in names(rows_by_id)) {
    rows <- rows_by_id[[nm]]
    ev <- tags$event[rows]
    ev_int <- match(ev, unique(ev))
    ev_new[rows] <- ev_int

    sizes <- tabulate(ev_int, nbins = max(ev_int))
    if (sizes[1] > 1) {
      orphan <- c(orphan, nm)
      next()
    }
    for (k in which(sizes > 1)) {
      kk <- rows[ev_int == k]
      sk <- sum(tags$prob[kk])
      if (sk > 0 && abs(sk - 1) > 1e-9) {
        pr_new[kk] <- tags$prob[kk] / sk
        rescaled <- c(rescaled, nm)
      }
    }
  }

  tags$event <- ev_new
  tags$prob <- pr_new

  if (!is.null(rescaled)) {
    rescaled <- unique(rescaled)
    if (verbose) message(length(rescaled), " tag", if (length(rescaled) == 1) "" else "s",
                         " had candidate recapture positions removed by the checks; the ",
                         "probabilities of the remaining candidates were rescaled to sum to 1 (id",
                         if (length(rescaled) == 1) "" else "s", ": ", .format_ids(rescaled), ").")
  }
  if (!is.null(orphan)) {
    if (verbose) message(length(orphan), " tag", if (length(orphan) == 1) "" else "s",
                         " removed because the release observation was dropped by the checks, ",
                         "leaving only ambiguous candidate positions with nothing to propagate ",
                         "from (id", if (length(orphan) == 1) "" else "s", ": ",
                         .format_ids(orphan), ").")
    tags <- tags[!(tags$id %in% orphan), , drop = FALSE]
  }

  tags
}


## Mark-recapture tags carry exactly one release and one recovery, so a "c" tag
## with more than two observation EVENTS means two physical tags share an id.
## (Several rows within one event are the candidate positions of a single
## ambiguous recovery and are counted once here -- see .check_events().) Unlike
## the cross-type case this cannot be repaired automatically -- there is no way
## to tell which row belongs to which tag -- so warn and leave the data alone
## rather than silently dropping observations.
.check_ctag_rows <- function(tags, verbose = TRUE) {

  if (is.null(tags) || nrow(tags) == 0) return(invisible(NULL))
  if (!any(tags$tag_type %in% "c")) return(invisible(NULL))

  ctg <- tags[tags$tag_type %in% "c", , drop = FALSE]
  if (any(colnames(ctg) == "event")) {
    n <- tapply(ctg$event, as.character(ctg$id), function(z) length(unique(z)))
  } else {
    n <- table(as.character(ctg$id))
  }
  bad <- names(n)[n > 2]

  if (length(bad) > 0 && verbose) {
    warning(length(bad), " mark-recapture tag id",
            if (length(bad) == 1) "" else "s",
            " with more than two observations (", .format_ids(bad),
            "). A mark-recapture tag has exactly one release and one recovery, ",
            "so these ids are most likely shared by several tags. They will be ",
            "fitted as single tags with interleaved positions. Please check the ",
            "ids in the input data.", call. = FALSE)
  }

  invisible(NULL)
}


##' Attach candidate recapture locations to tags
##'
##' @description
##' `add_candidates()` replaces the final observation of each tag with a set of
##' candidate positions, turning a single assumed recapture into an ambiguous
##' one. This is the situation where the recapturing vessel is known but the
##' individual set is not: each of the vessel's fishing sets is a possible
##' recapture location, with a probability derived from its effort.
##'
##' Use this when the candidates arrive as their own long table, one row per
##' candidate, keyed by tag. When they arrive as repeated columns of a wide
##' table instead (`date1`, `lat1`, `lon1`, `per1`, `date2`, ...), map them
##' directly in [prep_tags()] via `candidates` and skip this function.
##'
##' @param x An object of class `admove_tags` (or one containing tags, such as
##'   `admove_data`).
##' @param candidates A data frame of candidate positions, one row per
##'   candidate, with a column identifying the tag.
##' @param by Named character scalar linking tags to candidates, given as
##'   `c(<tags column> = "<candidates column>")`, e.g. `c(id = "REC_ID")`. A
##'   single unnamed string is used for both sides.
##' @param names Named character vector giving the candidate columns:
##'   `c(t = "...", x = "...", y = "...", prob = "...")`. `prob` is optional; if
##'   omitted the candidates of a tag are treated as equally likely.
##' @param verbose Logical; if `TRUE`, report how many tags were matched.
##'   Default: `TRUE`.
##'
##' @return
##' An object of class `admove_tags` carrying `event` and `prob` columns. Tags
##' with no entry in `candidates` are returned unchanged, with their single
##' recapture at probability 1.
##'
##' @details
##' Probabilities must already sum to 1 within each tag; convert an effort
##' measure beforehand, for example `prob = hours / sum(hours)` per tag. The
##' check itself happens in [check_tags()], together with the rule that only the
##' final observation of a tag may be ambiguous.
##'
##' @seealso [prep_tags()], [check_tags()]
##'
##' @examples
##' ## tags <- add_candidates(ctags, sets,
##' ##                        by = c(id = "REC_ID"),
##' ##                        names = c(t = "set_date", x = "set_lon",
##' ##                                  y = "set_lat", prob = "set_prob"))
##'
##' @export
add_candidates <- function(x, candidates, by = c(id = "id"),
                           names = NULL, verbose = TRUE) {

  if (inherits(x, "admove_data")) {
    tags <- x$tags
  } else {
    tags <- x
  }
  if (is.null(tags) || nrow(tags) == 0) stop("No tags found!")
  if (is.null(candidates) || nrow(candidates) == 0) stop("'candidates' is empty.")

  if (is.null(names) || is.null(base::names(names))) {
    stop("Please provide 'names' as a named vector giving the candidate columns, e.g. names = c(t = 'set_date', x = 'set_lon', y = 'set_lat', prob = 'set_prob').")
  }
  req <- c("t","x","y")
  if (!all(req %in% base::names(names))) {
    stop("'names' must contain at least t, x and y. See ?add_candidates.")
  }

  if (length(by) != 1) stop("'by' must link exactly one column, e.g. by = c(id = 'REC_ID').")
  by_tags <- if (is.null(base::names(by))) as.character(by) else base::names(by)
  by_cand <- as.character(by)

  if (!any(colnames(tags) == by_tags)) stop("Column '", by_tags, "' not found in the tags.")
  if (!any(colnames(candidates) == by_cand)) stop("Column '", by_cand, "' not found in 'candidates'.")

  miss <- setdiff(as.character(names[intersect(c(req, "prob"), base::names(names))]),
                  colnames(candidates))
  if (length(miss) > 0) stop("Column(s) named in 'names' but not found in 'candidates': ", paste(miss, collapse = ", "), ".")

  ## make sure the existing event structure is present and canonical
  tags <- .check_events(tags, verbose = FALSE)

  cand <- data.frame(.key = as.character(candidates[[by_cand]]),
                     t = candidates[[names[["t"]]]],
                     x = candidates[[names[["x"]]]],
                     y = candidates[[names[["y"]]]],
                     stringsAsFactors = FALSE)
  cand$prob <- if (!is.null(names[["prob"]])) candidates[[names[["prob"]]]] else NA_real_

  keys <- as.character(tags[[by_tags]])
  matched <- unique(keys[keys %in% cand$.key])

  if (length(matched) == 0) {
    if (verbose) message("None of the tags appear in 'candidates'; nothing changed.")
    return(x)
  }

  rows_by_id <- split(seq_len(nrow(tags)), keys)
  out <- vector("list", length(rows_by_id))
  nms <- base::names(rows_by_id)

  for (j in seq_along(rows_by_id)) {

    rows <- rows_by_id[[j]]
    tg <- tags[rows, , drop = FALSE]

    if (!(nms[j] %in% matched)) {
      out[[j]] <- tg
      next()
    }

    cj <- cand[cand$.key == nms[j], , drop = FALSE]

    ## the last event is the one being replaced; keep everything before it
    last_ev <- max(tg$event)
    head_rows <- tg[tg$event != last_ev, , drop = FALSE]
    template <- tg[tg$event == last_ev, , drop = FALSE][1, , drop = FALSE]

    new <- template[rep(1L, nrow(cj)), , drop = FALSE]
    new$t <- cj$t
    new$x <- cj$x
    new$y <- cj$y
    new$event <- last_ev
    new$prob <- if (all(is.na(cj$prob))) rep(1 / nrow(cj), nrow(cj)) else cj$prob

    out[[j]] <- rbind(head_rows, new)
  }

  tags_out <- do.call(rbind, out)
  rownames(tags_out) <- NULL

  if (verbose) {
    message(length(matched), " tag", if (length(matched) == 1) "" else "s",
            " given candidate recapture positions (",
            length(rows_by_id) - length(matched), " left with a single recapture); ",
            "median ", stats::median(as.numeric(table(cand$.key[cand$.key %in% matched]))),
            " candidates per tag.")
  }

  tags_out <- .add_class(tags_out, "admove_tags")
  tags_out <- add_sref(tags_out, sref(tags))
  tags_out <- add_tref(tags_out, tref(tags))

  if (inherits(x, "admove_data")) {
    x$tags <- tags_out
    return(x)
  }
  tags_out
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

  cn <- colnames(x)

  ## Candidate recapture slots. Either the plain t1/x1/y1 (and optionally p1) of
  ## a single, known recapture, or the indexed t1.k/x1.k/y1.k/p1.k that
  ## prep_tags() writes when several candidate recapture locations were mapped
  ## -- the layout of uncertainty tables such as the IATTC
  ## date1/lat1/lon1/per1, date2, ... blocks.
  slot_cols <- function(stem) {
    if (stem %in% cn) return(stem)
    hit <- grep(paste0("^", stem, "\\.[0-9]+$"), cn, value = TRUE)
    hit[order(as.integer(sub(paste0("^", stem, "\\."), "", hit)))]
  }
  t_cand <- slot_cols("t1")
  x_cand <- slot_cols("x1")
  y_cand <- slot_cols("y1")
  p_cand <- slot_cols("p1")

  nk <- length(t_cand)
  if (length(x_cand) != nk || length(y_cand) != nk) {
    stop("The wide format needs the same number of columns for t1, x1 and y1, but found ",
         nk, ", ", length(x_cand), " and ", length(y_cand),
         ". Check the 'names' argument of prep_tags().")
  }
  if (length(p_cand) > 0 && length(p_cand) != nk) {
    stop("'p1' maps ", length(p_cand), " column(s) but t1/x1/y1 map ", nk,
         ". Give one probability column per candidate recapture location.")
  }

  ## Only carry the candidate machinery when it is actually used, so a classic
  ## single-recapture input returns exactly the two-row frame it always did.
  flag_cand <- nk > 1 || length(p_cand) > 0

  used <- c("t0","x0","y0","id", t_cand, x_cand, y_cand, p_cand)
  other_cols <- which(!cn %in% used)

  ## Wide to long, assembled column-wise rather than row-wise. Building whole
  ## columns keeps the class of t0/t1 intact (they are often Date, and
  ## prep_tags() parses them further down), which a per-row c() or unlist()
  ## would silently strip. It is also far faster on the ~10^4-10^5 row tag
  ## tables these files come in.
  ntags <- nrow(x)
  blocks <- vector("list", nk + 1L)

  mk <- function(tt, xx, yy, pp, slot) {
    out <- data.frame(t = tt, x = xx, y = yy, id = x[["id"]])
    if (flag_cand) {
      out$event <- if (slot == 0L) 1L else 2L
      out$prob <- pp
    }
    if (length(other_cols) > 0) {
      tmp <- as.data.frame(x[, other_cols, drop = FALSE])
      colnames(tmp) <- colnames(x)[other_cols]
      out <- cbind(out, tmp)
    }
    out$.tag <- seq_len(ntags)
    out$.slot <- slot
    out
  }

  blocks[[1]] <- mk(x[["t0"]], x[["x0"]], x[["y0"]], rep(1, ntags), 0L)
  for (k in seq_len(nk)) {
    blocks[[k + 1L]] <- mk(x[[t_cand[k]]], x[[x_cand[k]]], x[[y_cand[k]]],
                           if (length(p_cand) > 0) x[[p_cand[k]]] else rep(1, ntags),
                           k)
  }

  res <- do.call(rbind, blocks)

  ## An empty candidate slot means "this tag has fewer candidates than the table
  ## is wide", which is how such tables pad to a fixed number of columns. Drop
  ## those, but keep one (NA) recapture row for a tag that has none at all, so a
  ## non-recovered tag is still removed downstream by check_tags() exactly as
  ## before.
  empty <- res$.slot > 0L & is.na(res$t) & is.na(res$x) & is.na(res$y)
  if (any(empty)) {
    n_left <- tapply(!empty & res$.slot > 0L, res$.tag, sum)
    orphan <- as.integer(names(n_left)[n_left == 0])
    keep <- !empty | (res$.slot == 1L & res$.tag %in% orphan)
    res <- res[keep, , drop = FALSE]
  }

  res <- res[order(res$.tag, res$.slot), , drop = FALSE]
  res$.tag <- NULL
  res$.slot <- NULL
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

  ## tags whose final observation has several candidate positions
  n_amb <- lapply(tags_split2, function(z) {
    if (length(z) == 0) return(NA)
    sum(sapply(z, function(w) {
      ev <- .tag_events(w)
      sum(ev == max(ev)) > 1L
    }))
  })
  k_amb <- lapply(tags_split2, function(z) {
    if (length(z) == 0) return(NA)
    k <- sapply(z, function(w) {
      ev <- .tag_events(w)
      sum(ev == max(ev))
    })
    if (all(k < 2)) NA else mean(k[k > 1])
  })

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
    if (!is.na(n_amb[[i]]) && n_amb[[i]] > 0) {
      cat(sprintf(paste0("  %-", labw, "s %s\n"), "ambiguous:",
                  paste0(n_amb[[i]], " tag", if (n_amb[[i]] == 1) "" else "s",
                         " with ", sprintf("%.1f", k_amb[[i]]),
                         " candidate positions on average")))
    }
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
    tags <- x$tags
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
        ## The merged tag holds the recoveries of SEVERAL physical tags that
        ## happen to share a release, so they are successive observations of
        ## the pooled tag -- not candidate positions for one of them. Renumber
        ## the events accordingly; without this every recovery would inherit
        ## event 2 from its original tag and the likelihood would read them as
        ## mutually exclusive alternatives (see .check_events()).
        if (!is.null(tmp_comb[["event"]])) {
          tmp_comb$event <- seq_len(nrow(tmp_comb))
          tmp_comb$prob <- 1
        }
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
##' @param xlab Label for the x-axis. If `NULL` (default), `"x"` with the
##'   spatial units in brackets, e.g. `"x [km]"` (`"lon [°]"` for degrees).
##' @param ylab Label for the y-axis. If `NULL` (default), `"y"` with the
##'   spatial units in brackets, e.g. `"y [km]"` (`"lat [°]"` for degrees).
##' @param leg_pos Position of the legend. Default: `"topright"`.
##' @param labels Logical; if `TRUE`, label observations by time instead of
##'   plotting intermediate points. Default: `FALSE`.
##' @param bg Optional background colour for the plot. Default: `NULL`.
##' @param by_tag_type Logical; if `TRUE`, create separate panels by tag type.
##'   Default: `TRUE`.
##' @param by_tag Logical; if `TRUE`, create separate panels for individual tags.
##'   Default: `FALSE`.
##' @param show Character vector naming the elements to draw: any of
##'   `"release"` (release positions), `"recovery"` (recovery or final
##'   observation positions) and `"path"` (trajectories, intermediate
##'   observations and release-recovery segments). Legend entries of elements
##'   not drawn are dropped. Default: all three.
##' @param col Character vector of length 1 to 3 giving colours for tag paths,
##'   release positions, and recovery or final observation positions.
##' @param pch Integer vector of length 1 to 3 giving plotting symbols for
##'   intermediate observations, release positions, and recovery or final
##'   observation positions. Default: `c(1, 0, 16)`.
##' @param cex Numeric character expansion factor for plotted points.
##'   Default: `0.8`.
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
##' ## release and recovery positions in separate panels
##' op <- par(mfrow = c(1, 2))
##' plot_tags(skjepo$sim$tags, show = "release", main = "Release",
##'           by_tag_type = FALSE, auto_layout = FALSE)
##' plot_tags(skjepo$sim$tags, show = "recovery", main = "Recovery",
##'           by_tag_type = FALSE, auto_layout = FALSE)
##' par(op)
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
                      xlab = NULL,
                      ylab = NULL,
                      leg_pos = "topright",
                      labels = FALSE,
                      bg = NULL,
                      by_tag_type = TRUE,
                      by_tag = FALSE,
                      show = c("release", "recovery", "path"),
                      col = c(adjustcolor("grey60",0.3), .admove_cols(2)),
                      pch = c(1,0,16),
                      cex = 0.8,
                      ...) {

  map_labs <- .map_labs(x)
  if (is.null(xlab)) xlab <- map_labs[1]
  if (is.null(ylab)) ylab <- map_labs[2]

  show <- match.arg(show, several.ok = TRUE)
  show_release <- "release" %in% show
  show_recovery <- "recovery" %in% show
  show_path <- "path" %in% show

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

  if (!inherits(tags, "data.frame")) {
    tags <- do.call(rbind, tags)
  }

  ## Work on the long table with per-tag bookkeeping done vectorised: splitting
  ## the data frame into one data frame per tag costs minutes for tens of
  ## thousands of tags. Tags are drawn in the order split() would give.
  idf <- factor(tags$id)
  o <- order(as.integer(idf))
  ti <- as.integer(idf)[o]  ## tag index of each row, 1..ntag, non-decreasing
  ntag <- nlevels(idf)
  tt <- .subset2(tags, 1L)[o]
  tx <- .subset2(tags, 2L)[o]
  ty <- .subset2(tags, 3L)[o]
  first_row <- which(!duplicated(ti))
  last_row <- which(!duplicated(ti, fromLast = TRUE))
  nrow_tag <- tabulate(ti, ntag)
  tag_type_int <- .get_tag_type_integer(tags$tag_type[o][first_row])

  ## Rows of the final observation event, as .tag_events(): events are numbered
  ## by first appearance within a tag, so the last event is the one that
  ## appears last.
  if (is.null(tags[["event"]])) {
    is_last <- seq_along(ti) %in% last_row
  } else {
    key <- paste(ti, tags[["event"]][o], sep = "\r")
    ev_new <- !duplicated(key)
    ev_ord <- cumsum(ev_new)[match(key, key)]
    is_last <- ev_ord == cumsum(ev_new)[last_row][ti]
  }
  ## The final observation may be ambiguous: several candidate positions, one
  ## of which is the true recovery. Represent the tag by its most likely
  ## candidate (first one on ties) and draw the alternatives as a fan.
  cand <- which(is_last)
  cand_prob <- .na_zero(tags[["prob"]][o][cand])
  if (length(cand_prob) != length(cand)) cand_prob <- rep(0, length(cand))
  oc <- order(ti[cand], -cand_prob, cand)
  best <- oc[!duplicated(ti[cand][oc])]
  end_row <- cand[best]
  n_cand <- tabulate(ti[cand], ntag)
  fan <- n_cand[ti[cand]] >= 2L
  fan_row <- cand[fan]
  fan_pr <- cand_prob[fan]
  fan_max <- cand_prob[best][ti[fan_row]]
  fan_pr <- ifelse(fan_max <= 0, 1, fan_pr / fan_max)

  if (is.null(xlim)) {
    xlims <- range(tx, na.rm = TRUE)
  } else {
    xlims <- xlim
  }
  if (is.null(ylim)) {
    ylims <- range(ty, na.rm = TRUE)
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

  tag_type_labs <- c("Data-storage tags",
                     "Mark-resight tags",
                     "Mark-recapture tags",
                     "Acoustic tags")

  ## sel: logical per tag, the tags drawn in this panel
  plot_one <- function(sel) {

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
    st <- first_row[sel]
    en <- end_row[sel]
    if (show_release) {
      points(tx[st], ty[st],
             col = cols[2], pch = pch[1])
    }
    if (show_recovery) {
      points(tx[en], ty[en],
             col = cols[3], pch = pch[2])
    }

    if (show_path) {
      ## trajectories: few tags, so a per-tag loop over plain vectors is fine
      traj <- sel & tag_type_int %in% c(1,2,4)
      rows <- which(traj[ti])
      if (length(rows)) {
        g <- ti[rows]
        xs <- split(tx[rows], g)
        ys <- split(ty[rows], g)
        ts <- split(tt[rows], g)
        for (i in seq_along(xs)) {
          lines(xs[[i]], ys[[i]],
                col = cols[1], ty = "b", pch = NA)
          if (labels) {
            text(xs[[i]], ys[[i]],
                 labels = sprintf("%.2f", ts[[i]]),
                 col = cols[1], pch = pch[3], cex = cex)
          } else {
            nr <- length(xs[[i]])
            points(xs[[i]][-c(1,nr)], ys[[i]][-c(1,nr)],
                   col = cols[1], pch = pch[3], cex = cex)
          }
        }
      }

      segm <- sel & tag_type_int == 3
      if (any(segm)) {
        segments(tx[first_row[segm]], ty[first_row[segm]],
                 tx[end_row[segm]], ty[end_row[segm]],
                 col = cols[1])
      }
    }

    ## fan out the alternative recapture positions, shaded by probability
    k <- which(sel[ti[fan_row]])
    if (length(k)) {
      fr <- fan_row[k]
      pr <- fan_pr[k]
      if (show_path) {
        ## adjustcolor() takes a single alpha.f, so scale the alpha directly
        cc <- grDevices::col2rgb(cols[3], alpha = TRUE) / 255
        segments(tx[first_row[ti[fr]]], ty[first_row[ti[fr]]],
                 tx[fr], ty[fr],
                 col = grDevices::rgb(cc[1], cc[2], cc[3], cc[4] * pmax(0.15, pr)),
                 lty = 3, lwd = 0.5 + 1.5 * pr)
      }
      if (show_recovery) {
        points(tx[fr], ty[fr],
               col = cols[3], pch = pch[2], cex = cex * (0.5 + 0.8 * pr))
      }
    }

  }

  use_layout <- FALSE

  if (by_tag) {

    tag_types <- tag_type_int
    n <- ntag

    mfrow <- n2mfrow(n, asp = 2)
    if(auto_layout && !add){
      par(mfrow = mfrow,
          mar = c(0.1,0.1,0.1,0.1),
          oma = c(4,4,1,1))
      use_layout <- TRUE
    }
    main <- ""
  } else if (by_tag_type) {
    tag_types <- unique(tag_type_int)
    n <- length(tag_types)
    mfrow <- n2mfrow(n, asp = 2)
    if(auto_layout && !add && n > 1){
      par(mfrow = mfrow,
          mar = c(0.1,0.1,0.1,0.1),
          oma = c(4,4,1,1))
      use_layout <- TRUE
      main <- ""
    }
  }else {
    n <- 1

    mfrow <- c(1,1)
    if(auto_layout){
      par(mfrow = mfrow)
    }
    main <- main
  }

  for (i in 1:n) {

    if (use_layout) {
      xaxt <- ifelse(i %in% (prod(mfrow) - mfrow[2]+1):prod(mfrow), "s", "n")
      yaxt <- ifelse(i %in% seq(1, prod(mfrow), mfrow[2]), "s", "n")
    } else {
      xaxt <- "s"
      yaxt <- "s"
    }

    if (by_tag) {
      sel <- seq_len(ntag) == i
    } else if (by_tag_type) {
      sel <- tag_type_int == tag_types[i]
    } else {
      sel <- rep(TRUE, ntag)
    }

    plot_one(sel)

    if (i == 1 || by_tag_type) {
      if (by_tag_type || by_tag) {
      labo <- list(c("Deployment", "Intermediate obs.", "Recovery"),
                  c("Release", "Resights", "Final resight"),
                  c("Release", "Recovery"),
                  c("Release", "Detection"))[[as.integer(tag_types[i])]]
      roleo <- list(c("release", "path", "recovery"),
                    c("release", "path", "recovery"),
                    c("release", "recovery"),
                    c("release", "path"))[[as.integer(tag_types[i])]]
      } else if (any(nrow_tag[sel] > 2)){
        labo <- c("Release", "Intermediate obs.", "Final obs.")
        roleo <- c("release", "path", "recovery")
      } else {
        labo <- c("Release", "Recovery")
        roleo <- c("release", "recovery")
      }
      ## symbols and colours as drawn in plot_one()
      pcho <- c(release = pch[1], path = pch[3], recovery = pch[2])[roleo]
      colo <- c(release = cols[2], path = cols[1], recovery = cols[3])[roleo]
      keep <- roleo %in% show
      if (any(keep)) {
        legend(leg_pos,
               legend = labo[keep],
               pch = pcho[keep],
               col = colo[keep],
               bg = "white")
      }
    }
    box(lwd = 1.5)
    if (by_tag) {
      legend("topleft",
             legend = round(tt[first_row[i]],3),
             cex = 0.8,
             pch = NA,
             bg = "white")
    } else if (by_tag_type && n > 1) {
      legend("topleft",
             legend = tag_type_labs[as.integer(tag_types[i])],
             cex = 0.8,
             pch = NA,
             bg = "white")
    }
  }

  if (use_layout) {
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
                       candidates = NULL,
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
    candidates = candidates,
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
