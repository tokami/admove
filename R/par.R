##' Default initial parameter values for admove
##'
##' @description
##' Creates a named list of model parameters with default initial values for
##' estimation in \code{admove}.
##'
##' @param dat A data list containing model input data, as produced by
##'   [setup_data()].
##' @param conf An optional configuration list, typically created by
##'   [default_conf()]. If \code{NULL}, a default configuration is generated
##'   from \code{dat}.
##' @param cov_taxis Covariates that carry the habitat preference, given as
##'   names or indices into `dat$cov`. Only these are used to scale `logKappa`.
##'   Defaults to all covariates. Covariates that only enter the model as
##'   advection inputs -- typically the zonal and meridional current components,
##'   whose taxis coefficients are fixed in the `map` -- should be excluded:
##'   their range and gradient say nothing about the taxis scale, but they
##'   contribute to the median otherwise, so switching advection on and off
##'   changes `kappa` for no good reason.
##' @param verbose Logical; if \code{TRUE}, informative messages are printed.
##'
##' @details
##' This function generates a named list of initial parameter values based on
##' the supplied data and model configuration. These values can be used as
##' starting values for model fitting and are intended to provide reasonable
##' defaults for the selected model setup.
##'
##' The taxis scaling parameter `logKappa` is fixed during estimation (see
##' [default_map()]) and therefore its initial value is its final value. This is
##' not a restriction: the likelihood depends on the taxis only through
##' `kappa * grad(h)`, and `h` is linear in `alpha`, so `kappa * alpha` is the
##' only identifiable combination. `kappa` is therefore a pure scale factor, and
##' its only job is to put `alpha` on an \eqn{O(1)} scale so that the objective
##' is well conditioned. It cannot change the fitted movement field.
##'
##' Because `kappa` has units of \eqn{[\text{distance}^2 / \text{time}]}, a
##' value of 1 is only appropriate when coordinates are already on a unit scale.
##' The default is instead anchored on observed movement: a preference function
##' with unit range over the tag-sampled covariate range \eqn{R} has
##' \eqn{|dh/dcov| \approx 1/R}, so the drift accumulated over a time \eqn{T} is
##' \eqn{\kappa (G/R) T}, with \eqn{G} the typical covariate gradient at the tag
##' positions. Equating that to a characteristic displacement \eqn{L} gives
##'
##' \deqn{\kappa = L R / (G T),}
##'
##' where \eqn{L} and \eqn{T} are the median displacement and the median time
##' between successive observations, and \eqn{R} and \eqn{G} are evaluated at the
##' tag positions and reduced across the `cov_taxis` covariates by their median.
##' Only the tag types enabled in `conf` contribute, since `dat$tags` may still
##' carry types that `conf` switches off.
##'
##' \eqn{L} and \eqn{T} are computed **per tag type** and the resulting
##' \eqn{\kappa_{type}} combined by their geometric mean, rather than pooling all
##' steps into one median. Tag types differ in their step scale by construction
##' -- a mark-recapture tag contributes one net displacement over months, a
##' data-storage tag hundreds of daily increments -- so a median over the pooled
##' steps takes \eqn{T} from whichever type has more steps and \eqn{L} from the
##' mixture, and can land outside the range of the per-type values. The
##' geometric mean is the natural average for a multiplicative scale factor and
##' always lies between them.
##'
##' If the quantities above cannot be computed (no covariates, no tags, or a
##' covariate with no spatial variation), the function falls back to the earlier
##' rule `kappa = cellsize^2 / median_dt` and, failing that, to `kappa = 1`.
##'
##' Override via `par$logKappa <- log(<value>)` after calling this function. The
##' fit is insensitive to the exact value -- anything within roughly an order of
##' magnitude gives the same optimum -- so this only needs to be in the right
##' range, not finely tuned.
##'
##' @return
##' A named list of initial parameter values.
##'
##' @examples
##' par <- with(skjepo$sim, default_par(dat, conf))
##'
##' @export
default_par <- function(dat, conf = NULL, cov_taxis = NULL, verbose = TRUE) {

  if (is.null(conf)) {
    if (verbose) message("No configuration list provided, using default_conf(dat). ")
    conf <- default_conf(dat, verbose = verbose)
  }

  par <- list()

  ## Seasonal structure is specified in conf, so derive the spline breakpoints
  ## it implies before the array dimensions are taken from them
  dat <- .resolve_seasons(dat, conf)$dat

  ## Number of spline slices per covariate. The third dimension of alpha, beta
  ## and gamma is shared, so it is sized by the covariate that uses the most
  ## slices; default_map() fixes the slices each covariate never evaluates.
  max_seasonal <- max(.get_nsea(dat))

  ## Taxis -----------------------------------------

  if (is.null(dat$knots_tax)) {
    knots_tax <- matrix(NA, 1, 1)
  } else {
    knots_tax <- dat$knots_tax
  }

  par$alpha <- array(rep(0, length(knots_tax)),
                     dim = c(nrow(knots_tax),
                             ncol(knots_tax),
                             max_seasonal))


  ## Diffusion --------------------------------------

  if (is.null(dat$knots_dif)) {
    knots_dif <- matrix(NA, 1, 1)
  } else {
    knots_dif <- dat$knots_dif
  }

  par$beta <- array(rep(0, length(knots_dif)),
                    dim = c(nrow(knots_dif),
                            ncol(knots_dif),
                            max_seasonal))

  ## Advection ----------------------------------------

  if (is.null(dat$cov)) {
    cov <- 1
  } else {
    cov <- dat$cov
  }

  par$gamma <- array(rep(0, length(cov)),
                     dim = c(2,
                             length(cov),
                             max_seasonal))


  ## Taxis scaling -------------------------------------
  ## kappa is a pure scale factor (only kappa * alpha is identifiable), so its
  ## job is to put alpha on an O(1) scale. Anchor it on observed movement,
  ## kappa = L * R / (G * T), rather than on the integration step -- see the
  ## @details section above. Only the tag types enabled in conf are used, only
  ## the cov_taxis covariates enter the median across covariates, and the tag
  ## types are reduced separately and then combined.
  tags_use <- .get_tags_in_use(dat, conf)
  idx_tax <- .resolve_cov_taxis(cov_taxis, dat)
  steps <- .get_tag_steps(tags_use)

  kappa <- NA_real_
  kappa_type <- .kappa_by_tag_type(dat, tags_use, idx_tax)

  if (length(kappa_type) > 0) kappa <- exp(mean(log(kappa_type)))

  ## Fallback: the earlier grid- and time-step-based rule, then 1.
  if (!is.finite(kappa) || kappa <= 0) {
    cs <- if (!is.null(dat$grid)) dat$grid$cellsize[1] else 1
    med_dt <- if (!is.null(steps)) median(steps$dt[steps$dt > 0], na.rm = TRUE) else NA_real_
    if (!is.finite(med_dt) || med_dt <= 0) med_dt <- 1
    kappa <- cs^2 / med_dt
    if (!is.finite(kappa) || kappa <= 0) kappa <- 1
    if (verbose) {
      message("Could not derive kappa from the tag displacements and covariate ",
              "gradients; falling back to cellsize^2 / median_dt = ",
              signif(kappa, 4), ". Check par$logKappa.")
    }
  } else if (verbose) {
    per_type <- ""
    if (length(kappa_type) > 1) {
      per_type <- paste0(" (geometric mean of ",
                         paste0(.tag_type_label(names(kappa_type)), ": ",
                                signif(kappa_type, 4), collapse = ", "), ")")
    }
    message("kappa set to ", signif(kappa, 4),
            " from the observed tag displacements and covariate gradients",
            per_type, ".")
  }

  par$logKappa <- log(kappa)



  ## Observation uncertainty ---------------------------
  ## 3 tag types
  par$logSdO <- matrix(rep(0,6), 2, 3)

  ## return
  par
}


## Internal helpers for the kappa default ------------------------------------

## Human-readable names for the tag type letters, for messages.
.tag_type_label <- function(ty) {
  lab <- c(d = "data-storage", s = "mark-resight", c = "mark-recapture")
  out <- unname(lab[as.character(ty)])
  out[is.na(out)] <- as.character(ty)[is.na(out)]
  out
}


## Covariates that carry the habitat preference, as indices into dat$cov.
## Everything else (typically the current components used for advection only)
## is excluded from the kappa scale.
.resolve_cov_taxis <- function(cov_taxis, dat) {

  ncov <- if (!is.null(dat$cov)) length(dat$cov) else 0L
  if (ncov == 0L) return(integer(0))
  if (is.null(cov_taxis)) return(seq_len(ncov))

  if (is.character(cov_taxis)) {
    idx <- match(cov_taxis, names(dat$cov))
    if (anyNA(idx)) {
      stop("Unknown covariate(s) in 'cov_taxis': ",
           paste(cov_taxis[is.na(idx)], collapse = ", "),
           if (!is.null(names(dat$cov)))
             paste0(". Available: ", paste(names(dat$cov), collapse = ", ")),
           call. = FALSE)
    }
  } else {
    idx <- suppressWarnings(as.integer(cov_taxis))
    if (anyNA(idx) || any(idx < 1L) || any(idx > ncov)) {
      stop("'cov_taxis' must index covariates 1:", ncov, ".", call. = FALSE)
    }
  }

  sort(unique(idx))
}


## kappa = L * R / (G * T) evaluated separately for each tag type present, named
## by the tag type letter. Pooling the steps of different tag types would take T
## from whichever type contributes the most steps (a handful of data-storage
## tags sampled daily outnumber thousands of mark-recapture displacements) and L
## from the mixture of both, giving a value that need not lie between the
## per-type ones.
.kappa_by_tag_type <- function(dat, tags, idx) {

  out <- numeric(0)

  if (is.null(tags) || nrow(tags) == 0 || is.null(dat$cov) ||
        is.null(dat$time_cov) || length(idx) == 0) {
    return(out)
  }

  for (ty in unique(as.character(tags$tag_type))) {

    tt <- tags[as.character(tags$tag_type) == ty, , drop = FALSE]

    steps <- .get_tag_steps(tt)
    if (is.null(steps)) next

    tl <- median(steps$dt[steps$dt > 0], na.rm = TRUE)
    ll <- median(steps$dl[steps$dl > 0], na.rm = TRUE)
    if (!is.finite(tl) || tl <= 0 || !is.finite(ll) || ll <= 0) next

    cs_sum <- .cov_scales_at_tags(dat, tt, idx)
    ki <- ll * cs_sum$range / (cs_sum$grad * tl)
    ki <- ki[is.finite(ki) & ki > 0]
    if (length(ki) > 0) out[ty] <- median(ki)
  }

  out
}

## Subset the tags to the types enabled in conf. dat$tags may still carry types
## that conf switches off (check_tags() only drops them inside admove()), and
## including them skews any scale derived from the tag timing -- e.g. a handful
## of data-storage tags sampled hourly would otherwise dominate the median time
## step of thousands of mark-recapture tags.
.get_tags_in_use <- function(dat, conf) {

  if (is.null(dat$tags) || nrow(dat$tags) == 0) return(NULL)

  keep <- c("d", "s", "c")[c(isTRUE(conf$use_dtags),
                             isTRUE(conf$use_stags),
                             isTRUE(conf$use_ctags))]

  if (length(keep) == 0) return(dat$tags)

  out <- dat$tags[dat$tags$tag_type %in% keep, , drop = FALSE]
  if (nrow(out) == 0) return(dat$tags)

  out
}


## Time steps and displacements between successive observations of each tag.
## Grouped by tag type as well as id: default_par() runs before check_tags(),
## so an id shared between tag types has not been disambiguated yet and would
## otherwise merge two tags into one spurious track.
.get_tag_steps <- function(tags) {

  if (is.null(tags) || nrow(tags) < 2) return(NULL)

  grp <- paste0(as.character(tags$tag_type), "-", as.character(tags$id))

  out <- lapply(split(tags, grp), function(tg) {
    if (nrow(tg) < 2) return(NULL)
    o <- order(tg$t)
    data.frame(dt = diff(tg$t[o]),
               dl = sqrt(diff(tg$x[o])^2 + diff(tg$y[o])^2))
  })

  out <- do.call(rbind, out)
  if (is.null(out) || nrow(out) == 0) return(NULL)

  out
}


## Central differences of a matrix along its first dimension, one-sided at the
## edges. Vectorised: this is called once per covariate time slice, so a
## per-cell loop (as in .dxfield()) would be prohibitively slow here.
.grad_x <- function(m, d) {

  nr <- nrow(m)
  if (nr < 2 || !is.finite(d) || d == 0) return(array(NA_real_, dim = dim(m)))

  g <- (rbind(m[-1,, drop = FALSE], NA) -
          rbind(NA, m[-nr,, drop = FALSE])) / (2 * d)
  g[1,] <- (m[2,] - m[1,]) / d
  g[nr,] <- (m[nr,] - m[nr - 1,]) / d

  g
}


.grad_y <- function(m, d) t(.grad_x(t(m), d))


## Per covariate: the range of the values the tags actually sampled and the
## median gradient magnitude at the tag positions. Uses nearest-cell lookup
## rather than interpolation -- kappa only has to be right to within an order of
## magnitude, and this avoids building interpolators for every time slice.
## `idx` restricts the work to the covariates that carry the taxis; the returned
## vectors follow `idx`.
.cov_scales_at_tags <- function(dat, tags, idx = seq_along(dat$cov)) {

  rng <- rep(NA_real_, length(idx))
  grd <- rep(NA_real_, length(idx))

  if (is.null(tags) || nrow(tags) == 0) return(list(range = rng, grad = grd))

  for (k in seq_along(idx)) {

    i <- idx[k]
    a <- unclass(dat$cov[[i]])
    dn <- dimnames(a)
    if (is.null(dn)) next

    xc <- suppressWarnings(as.numeric(dn[[1]]))
    yc <- suppressWarnings(as.numeric(dn[[2]]))
    if (length(xc) < 2 || length(yc) < 2) next
    if (any(!is.finite(xc)) || any(!is.finite(yc))) next

    dx <- mean(diff(xc))
    dy <- mean(diff(yc))

    ix <- .nearest_index(tags$x, xc)
    iy <- .nearest_index(tags$y, yc)
    it <- as.integer(t2index(tags$t, dat$time_cov[[i]]))

    vals <- rep(NA_real_, nrow(tags))
    gmag <- rep(NA_real_, nrow(tags))

    for (j in sort(unique(it[it > 0]))) {
      rows <- which(it == j & !is.na(ix) & !is.na(iy))
      if (length(rows) == 0) next
      m <- a[,,j]
      ind <- cbind(ix[rows], iy[rows])
      vals[rows] <- m[ind]
      gmag[rows] <- sqrt(.grad_x(m, dx)[ind]^2 + .grad_y(m, dy)[ind]^2)
    }

    if (any(is.finite(vals))) rng[k] <- diff(range(vals[is.finite(vals)]))
    gok <- gmag[is.finite(gmag) & gmag > 0]
    if (length(gok) > 0) grd[k] <- median(gok)
  }

  list(range = rng, grad = grd)
}


## Index of the nearest cell centre, clamped to the field.
.nearest_index <- function(v, centres) {

  step <- mean(diff(centres))
  if (!is.finite(step) || step == 0) return(rep(NA_integer_, length(v)))

  idx <- as.integer(round((v - centres[1]) / step)) + 1L
  idx[!is.finite(idx)] <- NA_integer_

  pmin(pmax(idx, 1L), length(centres))
}




##' Check parameter dimensions for admove
##'
##' @description
##' Checks whether the parameters supplied in \code{par} have the expected
##' names and dimensions for the provided data and model configuration.
##'
##' @param par A named list of model parameters to be checked.
##' @param dat A data list containing model input data, as produced by
##'   [setup_data()].
##' @param conf An optional configuration list, typically created by
##'   [default_conf()]. If \code{NULL}, a default configuration is generated
##'   from \code{dat}.
##' @param verbose Logical; if \code{TRUE}, informative messages are printed
##'   when \code{conf} is generated internally.
##'
##' @details
##' The function constructs the expected parameter structure using
##' [default_par()] and compares it with the user-supplied \code{par} list.
##' It checks that all required parameters are present, flags unexpected
##' parameters, and verifies that each parameter has the correct dimensions.
##'
##' If any mismatch is found, the function stops with an informative error
##' message describing the problem.
##'
##' @return
##' Invisibly returns \code{TRUE} if all parameter names and dimensions are
##' valid.
##'
##' @examples
##' ## If all checks passed, the function returns an invisible TRUE:
##' with(skjepo$sim, check_par(par, dat, conf))
##'
##' ## If there is a problem, the function returns an error:
##' \dontrun{
##' par <- with(skjepo$sim, default_par(dat, conf))
##' par$alpha <- matrix(0, 3, 1)
##' check_par(par, skjepo$sim$dat, skjepo$sim$conf)
##' }
##'
##' @export
check_par <- function(par, dat, conf = NULL, verbose = TRUE) {

  if (is.null(conf)) {
    if (verbose) {
      message("No configuration list provided, using default_conf(dat).")
    }
    conf <- default_conf(dat, verbose = verbose)
  }

  ## expected parameter structure from defaults
  par0 <- default_par(dat = dat, conf = conf, verbose = FALSE)

  ## collect errors
  errors <- character()

  ## check that par is a list
  if (!is.list(par)) {
    stop("'par' must be a list.", call. = FALSE)
  }

  ## check missing parameters
  missing_par <- setdiff(names(par0), names(par))
  if (length(missing_par) > 0) {
    errors <- c(
      errors,
      paste0(
        "Missing parameter(s) in 'par': ",
        paste(missing_par, collapse = ", ")
      )
    )
  }

  ## check unexpected parameters
  extra_par <- setdiff(names(par), names(par0))
  if (length(extra_par) > 0) {
    errors <- c(
      errors,
      paste0(
        "Unknown parameter(s) in 'par': ",
        paste(extra_par, collapse = ", ")
      )
    )
  }

  ## only check dimensions for parameters that exist in both
  common_par <- intersect(names(par0), names(par))

  for (nm in common_par) {

    x <- par[[nm]]
    x0 <- par0[[nm]]

    ## compare dimensions
    dx <- dim(x)
    dx0 <- dim(x0)

    ## vectors/scalars may have NULL dim
    lx <- length(x)
    lx0 <- length(x0)

    if (is.null(dx0) && is.null(dx)) {
      ## both are plain vectors/scalars: compare length
      if (!identical(lx, lx0)) {
        errors <- c(
          errors,
          paste0(
            "Parameter '", nm, "' has wrong length. Expected ",
            lx0, " but got ", lx, "."
          )
        )
      }
    } else if (!identical(dx, dx0)) {
      errors <- c(
        errors,
        paste0(
          "Parameter '", nm, "' has wrong dimensions. Expected ",
          paste(dx0, collapse = " x "),
          " but got ",
          if (is.null(dx)) {
            paste0("length ", lx)
          } else {
            paste(dx, collapse = " x ")
          },
          "."
        )
      )
    }
  }

  ## stop if any problem was found
  if (length(errors) > 0) {
    stop(paste(errors, collapse = "\n"), call. = FALSE)
  }

  invisible(TRUE)
}
