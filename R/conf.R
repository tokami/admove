##' Default model configuration for admove
##'
##' @description
##' Creates a configuration list with default settings for the \code{admove}
##' model based on the supplied input data.
##'
##' @param dat A data list containing model input data, as produced by
##'   [setup_data()].
##' @param n_seasons Number of seasons per covariate, as a single number applied
##'   to all covariates or one value per covariate. \code{1} (the default) means
##'   no seasonality. Values above 1 split the seasonal cycle
##'   (\code{period(dat)}) into that many equal seasons. See [set_seasons()].
##' @param verbose Logical; if \code{TRUE}, informative messages may be printed.
##'
##' @details
##' The default configuration is determined from the available data. In
##' particular, the function detects which tag types are present and sets
##' corresponding model flags.
##'
##' The returned list contains logical flags controlling which data sources and
##' movement components are used, settings for observation-variance estimation,
##' the estimation engine, the CTMC approximation method, the drift
##' discretisation scheme, and default seasonal settings for covariates and
##' spline effects.
##'
##' `drift_scheme` selects how the drift term (taxis *and* advection) is
##' discretised on the grid when assembling the generator. `"upwind"` (the
##' default) is first-order upstream: off-diagonal rates are guaranteed
##' non-negative, so the generator is always a valid CTMC generator, at the cost
##' of some numerical diffusion. `"central"` is second-order central difference:
##' it removes that numerical diffusion but can produce negative off-diagonal
##' rates when drift dominates diffusion (grid-Peclet > 2), yielding an invalid
##' generator and possibly negative transition probabilities. Use `"central"`
##' only when the grid is fine relative to the drift.
##'
##' Seasonality is specified here rather than in the data. `n_seasons` gives the
##' number of seasons per covariate (`1`, the default, means no seasonality) and
##' is the only setting needed to make a model seasonal: the breakpoints of the
##' seasonal spline basis are derived from it and from the cycle length
##' `period(dat)` when the model is set up, so they never have to be written into
##' `dat` by hand. `seasonal_spline` is filled in from `n_seasons` and is
##' retained for models that set breakpoints manually via `dat$time_spline`.
##' `seasonal_cov` is a separate setting controlling whether the covariate
##' *fields* themselves repeat over the cycle, and stays `FALSE` by default.
##'
##' `seasonal_dif` is a single logical controlling whether
##' diffusion is estimated season by season, and is **`FALSE` by default** even
##' when a seasonal spline basis is in use: diffusion is a second-moment
##' quantity, so estimating it per season divides the information available per
##' season and is easily traded off against a seasonal taxis acting on the same
##' covariate. Taxis and advection remain season-specific. See [default_map()]
##' for how these settings translate into fixed and estimated coefficients.
##'
##' Observation uncertainty is **off by default** (`obs_var_type = c(0L, 0L, 0L)`),
##' meaning tag locations are treated as exact. To estimate observation variance
##' for a tag type, set the corresponding element to `1L` (all but last
##' observation) or `2L` (all observations). For example, to estimate observation
##' variance for data-storage tags: `conf$obs_var_type[1] <- 1L`.
##'
##' The smooth used for the habitat preference functions is selected by
##' `conf$smooth_method`:
##' * `"natural"` (default) — a natural cubic spline through the knot values:
##'   piecewise cubic, twice continuously differentiable, with local support and
##'   **linear extrapolation** beyond the outer knots (robust in the covariate
##'   tails).
##' * `"poly"` — the legacy single global interpolating polynomial of degree
##'   `nknots - 1`; retained for reproducibility of older fits, but prone to
##'   oscillation and unbounded extrapolation.
##'
##' @return
##' A named list of default model configuration settings.
##'
##' @examples
##' conf <- default_conf(skjepo$sim$dat)
##'
##' ## four seasons per year for every covariate
##' conf <- default_conf(skjepo$sim$dat, n_seasons = 4)
##'
##' @export
default_conf <- function(dat, n_seasons = 1, verbose = TRUE) {

  flag_dtags <- !is.null(dat$tags) && any(dat$tags$tag_type == "d")
  flag_stags <- !is.null(dat$tags) && any(dat$tags$tag_type == "s")
  flag_ctags <- !is.null(dat$tags) && any(dat$tags$tag_type == "c")

  conf <- list()

  ## Flags
  conf$use_ctags <- flag_ctags
  conf$use_dtags <- flag_dtags
  conf$use_stags <- flag_stags
  conf$use_taxis <- TRUE
  conf$use_advection <- FALSE

  ## Observation uncertainty (per tag type: dtags, stags, ctags)
  ## 0 = not estimating observation uncertainty (default; treat locations as exact)
  ## 1 = estimating observation uncertainty for all but last observation
  ## 2 = estimating observation uncertainty for all observations
  ## 3 = fixing observation uncertainty to imported values (tag$sdx, tag$sdy)
  ## Enable estimation via conf$obs_var_type[1] <- 1L (dtags) etc.
  conf$obs_var_type <- c(0L, 0L, 0L)

  ## CTMC/KF updating
  ## 0 = no updating (default for ctags)
  ## 1 = updating (default for d and stags)
  conf$do_update <- c(TRUE, TRUE, FALSE)

  ## Estimation engine
  conf$engine <- 1

  ## CTMC method
  ## 0 = Matrix::expm
  ## 1 = expAv (uni = FALSE)
  ## 2 = expAv (uni = TRUE)
  conf$ctmc_method <- 0

  ## Discretisation of the drift term (taxis + advection) in the generator
  ## "upwind"  = first-order upstream; off-diagonal rates are always >= 0, so
  ##             the CTMC generator is always valid (default; robust).
  ## "central" = central difference; second order, no upwind numerical
  ##             diffusion, but off-diagonal rates can turn negative when
  ##             drift dominates diffusion (grid-Peclet > 2), which may yield
  ##             invalid generators / negative probabilities.
  conf$drift_scheme <- "upwind"

  ## Seasonality (no seasonality by default)
  if (!is.null(dat$cov)) {
    ncov <- length(dat$cov)
  } else {
    ncov <- 1
  }

  ## Number of seasons per covariate. 1 = no seasonality. The breakpoints of the
  ## seasonal spline basis are derived from this and from period(dat), so the
  ## seasonal structure of the model lives entirely in conf.
  conf$n_seasons <- .expand_n_seasons(n_seasons, ncov)

  conf$seasonal_cov <- rep(FALSE, ncov)
  conf$seasonal_spline <- conf$n_seasons > 1L

  if (any(conf$n_seasons > 1L)) {
    ## fail here rather than deep inside the likelihood
    .require_period(dat, conf$n_seasons)
  }

  ## Seasonal diffusion
  ## FALSE = diffusion spline coefficients are shared across seasons (default),
  ##         so only the taxis and advection fields vary seasonally.
  ## TRUE  = diffusion spline coefficients are estimated separately per season.
  ## Diffusion is a second-moment quantity, so splitting it by season divides
  ## the information available per season and is easily traded off against a
  ## seasonal taxis acting on the same covariate. It is therefore opt-in, even
  ## when conf$seasonal_spline enables a seasonal basis.
  conf$seasonal_dif <- FALSE

  ## Smooth used for the habitat preference functions
  ## "natural" = natural cubic spline (default): piecewise cubic, C2, local support,
  ##             linear extrapolation beyond the outer knots.
  ## "poly"    = legacy single global interpolating polynomial (retained for
  ##             reproducibility of older fits).
  conf$smooth_method <- "natural"

  conf
}


##' Check and complete a model configuration list
##'
##' @description
##' Checks whether a configuration list contains all required settings for
##' \code{admove}. Missing settings are filled in from [default_conf()].
##'
##' @param conf An optional configuration list. If \code{NULL}, the default
##'   configuration generated by [default_conf()] is returned.
##' @param dat A data list containing model input data, as produced by
##'   [setup_data()].
##' @param verbose Logical; if \code{TRUE}, informative messages are printed
##'   when missing settings are added or when default settings are used.
##'
##' @details
##' The function compares the user-supplied configuration list with the default
##' configuration produced by [default_conf()]. Any missing named entries are
##' added, while existing entries in \code{conf} are preserved unchanged.
##'
##' This function does not validate whether supplied configuration values are
##' internally consistent; it only ensures that all required settings are
##' present.
##'
##' @return
##' A named configuration list containing all required settings.
##'
##' @examples
##' conf <- list(use_taxis = FALSE)
##' conf <- check_conf(conf, skjepo$sim$dat)
##'
##' @export
check_conf <- function(conf = NULL, dat, verbose = TRUE) {

  conf_default <- default_conf(dat = dat, verbose = FALSE)

  ## If no configuration was supplied, return the defaults
  if (is.null(conf)) {
    if (verbose) {
      message("No configuration list provided, using default settings.")
    }
    return(conf_default)
  }

  ## Check input
  if (!is.list(conf)) {
    stop("'conf' must be a list or NULL.", call. = FALSE)
  }

  ## Fill in missing settings
  missing_names <- setdiff(names(conf_default), names(conf))

  if (length(missing_names) > 0) {
    for (nm in missing_names) {
      conf[[nm]] <- conf_default[[nm]]
    }

    if (verbose) {
      message(
        "Added missing configuration setting(s): ",
        paste(missing_names, collapse = ", ")
      )
    }
  }

  conf <- .check_seasonal_lengths(conf, dat)

  if (!is.null(conf$n_seasons)) {
    ncov <- if (!is.null(dat$cov)) length(dat$cov) else 1L
    conf$n_seasons <- .expand_n_seasons(conf$n_seasons, ncov)
    if (any(conf$n_seasons > 1L)) .require_period(dat, conf$n_seasons)
  }

  if (!is.null(conf$seasonal_dif) &&
        (length(conf$seasonal_dif) != 1L ||
           !is.logical(conf$seasonal_dif) ||
           is.na(conf$seasonal_dif))) {
    stop("'conf$seasonal_dif' must be a single TRUE or FALSE.", call. = FALSE)
  }

  if (!is.null(conf$drift_scheme) &&
        !identical(conf$drift_scheme, "upwind") &&
        !identical(conf$drift_scheme, "central")) {
    stop("'conf$drift_scheme' must be either \"upwind\" or \"central\", not ",
         deparse(conf$drift_scheme), ".", call. = FALSE)
  }

  ## Validate the preference smooth method
  if (!is.character(conf$smooth_method) || length(conf$smooth_method) != 1L ||
        !conf$smooth_method %in% c("natural", "poly")) {
    stop("'conf$smooth_method' must be one of \"natural\" or \"poly\".",
         call. = FALSE)
  }

  conf
}

##' Set the number of seasons of a model configuration
##'
##' @description
##' Makes a model seasonal by setting the number of seasons per covariate in a
##' configuration list. This is the only step needed to turn a fitted model
##' seasonal: the breakpoints of the seasonal spline basis are derived from the
##' number of seasons and the cycle length \code{period(dat)} when the model is
##' set up, so the data list never has to be edited.
##'
##' @param conf A configuration list, as produced by [default_conf()].
##' @param dat A data list containing model input data, as produced by
##'   [setup_data()]. Used to look up the covariates and the seasonal period.
##' @param n Number of seasons, as a single number or one value per selected
##'   covariate. \code{1} switches seasonality off again.
##' @param cov Covariates to apply \code{n} to, given as names or indices. If
##'   \code{NULL} (the default), \code{n} is applied to all covariates.
##' @param verbose Logical; if \code{TRUE}, a summary of the resulting seasonal
##'   structure is printed.
##'
##' @details
##' Seasons are equal divisions of the cycle, anchored at the time origin of
##' \code{dat}: with \code{n = 4} and a period of 12 months, the seasons run
##' \code{[0, 3)}, \code{[3, 6)}, \code{[6, 9)} and \code{[9, 12)} months after
##' the origin. The cycle length must be set on the data beforehand, e.g. through
##' \code{create_tref(period = 12)} or \code{period(dat) <- 12}.
##'
##' Seasonality of the covariate *fields* is a separate question, controlled by
##' \code{conf$seasonal_cov}: set that when a covariate is a climatology that
##' should repeat every cycle, rather than a time series covering the whole
##' study period.
##'
##' Diffusion stays constant across seasons unless \code{conf$seasonal_dif} is
##' set to \code{TRUE}; see [default_map()].
##'
##' @return
##' The configuration list with \code{n_seasons} and \code{seasonal_spline}
##' updated.
##'
##' @examples
##' dat <- skjepo$sim$dat
##' period(dat) <- 12
##' conf <- set_seasons(default_conf(dat), dat, n = 4)
##'
##' @export
set_seasons <- function(conf, dat, n, cov = NULL, verbose = TRUE) {

  if (!is.list(conf)) stop("'conf' must be a configuration list.", call. = FALSE)

  ncov <- if (!is.null(dat$cov)) length(dat$cov) else 1L
  cov_names <- names(dat$cov)

  if (is.null(cov)) {
    idx <- seq_len(ncov)
  } else if (is.character(cov)) {
    idx <- match(cov, cov_names)
    if (anyNA(idx)) {
      stop("Unknown covariate(s): ", paste(cov[is.na(idx)], collapse = ", "),
           if (!is.null(cov_names)) paste0(". Available: ", paste(cov_names, collapse = ", ")),
           call. = FALSE)
    }
  } else {
    idx <- as.integer(cov)
    if (anyNA(idx) || any(idx < 1L) || any(idx > ncov)) {
      stop("'cov' must index covariates 1:", ncov, ".", call. = FALSE)
    }
  }

  n <- .check_n_seasons_values(n, "n")
  if (length(n) == 1L) n <- rep(n, length(idx))
  if (length(n) != length(idx)) {
    stop("'n' has length ", length(n), " but ", length(idx),
         " covariate(s) were selected. Supply a single value or one per covariate.",
         call. = FALSE)
  }

  ns <- .expand_n_seasons(conf$n_seasons, ncov)
  ns[idx] <- n
  conf$n_seasons <- ns

  ## n > 1 drives the basis; n == 1 keeps a manually specified one, if any
  ss <- .expand_seasonal_spline(conf$seasonal_spline, ncov)
  conf$seasonal_spline <- ifelse(ns > 1L, TRUE, ss & .has_manual_breaks(dat, ncov))

  if (any(ns > 1L)) .require_period(dat, ns)

  if (isTRUE(verbose)) {
    per <- .get_period(dat)
    lab <- if (!is.null(cov_names)) cov_names else paste0("cov ", seq_len(ncov))
    sel <- ns > 1L
    if (!any(sel)) {
      message("Seasonality switched off for all covariates.")
    } else {
      message("Seasonal cycle of ", signif(per, 6), " ",
              tryCatch(units_time(dat), error = function(e) ""),
              " split into: ",
              paste0(lab[sel], ": ", ns[sel], " seasons", collapse = ", "), ".")
    }
  }

  conf
}


## Internal functions ---------------------------------------------------------------

## Validate a number-of-seasons specification without knowing the covariates yet
.check_n_seasons_values <- function(n, nm = "conf$n_seasons") {
  if (!is.numeric(n) || length(n) < 1L || any(!is.finite(n))) {
    stop("'", nm, "' must be one or more finite numbers.", call. = FALSE)
  }
  if (any(n < 1) || any(abs(n - round(n)) > 1e-8)) {
    stop("'", nm, "' must be whole numbers >= 1 (1 means no seasonality).",
         call. = FALSE)
  }
  as.integer(round(n))
}

## Recycle a number-of-seasons specification to one value per covariate
.expand_n_seasons <- function(n, ncov) {
  if (is.null(n)) return(rep(1L, ncov))
  n <- .check_n_seasons_values(n)
  if (length(n) == 1L) return(rep(n, ncov))
  if (length(n) != ncov) {
    stop("'conf$n_seasons' has length ", length(n), " but there ",
         if (ncov == 1L) "is 1 covariate" else paste0("are ", ncov, " covariates"),
         ". Supply a single value (recycled to all covariates) or one value per covariate.",
         call. = FALSE)
  }
  n
}

.expand_seasonal_spline <- function(v, ncov) {
  if (is.null(v)) return(rep(FALSE, ncov))
  if (length(v) == 1L) return(rep(as.logical(v), ncov))
  as.logical(v)
}

## Covariates whose spline breakpoints were written into dat by hand
.has_manual_breaks <- function(dat, ncov) {
  ts <- dat$time_spline
  if (is.null(ts) || length(ts) != ncov) return(rep(FALSE, ncov))
  vapply(ts, function(x) length(x) > 1L && !isTRUE(attr(x, "seasonal")),
         logical(1L))
}

## Seasonal cycle length of a data list, or NULL when none is defined
.get_period <- function(dat) {
  per <- tryCatch(period(dat), error = function(e) NULL)
  if (is.null(per) || length(per) != 1L) per <- dat$period
  if (is.null(per) || length(per) != 1L || !is.numeric(per) ||
        is.na(per) || !is.finite(per) || per <= 0) {
    return(NULL)
  }
  as.numeric(per)
}

.require_period <- function(dat, n_seasons) {
  per <- .get_period(dat)
  if (!is.null(per)) return(invisible(per))
  stop("A seasonal model was requested (n_seasons = ",
       paste(n_seasons, collapse = ", "),
       ") but no seasonal cycle length is defined for the data. Set one before ",
       "configuring the model, e.g. period(dat) <- 12 for monthly time units ",
       "with an annual cycle, or create_tref(..., period = 12) when preparing ",
       "the inputs.", call. = FALSE)
}

## Breakpoints of n equal seasons, as phases within the cycle starting at 0
.season_breaks <- function(period, n) {
  br <- (seq_len(n) - 1L) * period / n
  attr(br, "seasonal") <- TRUE
  br
}

## Derive the spline breakpoints implied by conf$n_seasons.
##
## Seasonality is specified in conf, in seasons rather than in time units, so
## that it survives changes of time origin and time unit and so that the same
## data list can be fitted with different seasonal structures. This translates
## that specification into the dat$time_spline breakpoints the likelihood
## indexes, and marks them as phases within the cycle rather than absolute
## times (see .shift_time_values()). It is idempotent: re-resolving an already
## resolved data list is a no-op.
.resolve_seasons <- function(dat, conf) {

  ncov <- if (!is.null(dat$cov)) length(dat$cov) else 1L
  ns <- .expand_n_seasons(conf$n_seasons, ncov)
  ss <- .expand_seasonal_spline(conf$seasonal_spline, ncov)

  if (any(ns > 1L)) {
    per <- .require_period(dat, ns)

    if (is.null(dat$time_spline) || length(dat$time_spline) != ncov) {
      dat$time_spline <- rep(list(0), ncov)
    }

    for (i in seq_len(ncov)) {
      if (ns[i] <= 1L) next
      old <- dat$time_spline[[i]]
      if (length(old) > 1L && !isTRUE(attr(old, "seasonal"))) {
        warning("Covariate ", i, " has manual breakpoints in dat$time_spline ",
                "and conf$n_seasons[", i, "] = ", ns[i],
                "; the configuration takes precedence and the manual ",
                "breakpoints are replaced.", call. = FALSE)
      }
      dat$time_spline[[i]] <- .season_breaks(per, ns[i])
      ss[i] <- TRUE
    }
  }

  ## Mark manually specified seasonal breakpoints too, so that they are also
  ## treated as phases rather than absolute times.
  if (!is.null(dat$time_spline) && length(dat$time_spline) == ncov) {
    for (i in seq_len(ncov)) {
      if (isTRUE(ss[i])) attr(dat$time_spline[[i]], "seasonal") <- TRUE
    }
  }

  conf$n_seasons <- ns
  conf$seasonal_spline <- ss

  list(dat = dat, conf = conf)
}

.check_seasonal_lengths <- function(conf, dat) {
  ncov <- if (!is.null(dat$cov)) length(dat$cov) else 1L
  for (nm in c("seasonal_spline", "seasonal_cov")) {
    v <- conf[[nm]]
    if (is.null(v)) next
    if (length(v) == 1L) {
      conf[[nm]] <- rep(v, ncov)
    } else if (length(v) != ncov) {
      stop("'conf$", nm, "' has length ", length(v), " but there ",
           if (ncov == 1L) "is 1 covariate" else paste0("are ", ncov, " covariates"),
           ". Supply a single value (recycled to all covariates) or one logical value per covariate.",
           call. = FALSE)
    }
  }
  conf
}
