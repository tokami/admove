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
##'   Defaults to all covariates. Covariates that only inform diffusion, with
##'   their taxis coefficients fixed in the `map`, should be excluded: their
##'   range and gradient say nothing about the taxis scale. Advection fields
##'   (`setup_data(adv = )`) are not covariates and never enter `kappa`.
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
##' The default instead balances the information the data carry about `alpha`
##' against that about the diffusion. A unit change of one `alpha` coefficient
##' changes the preference slope by about \eqn{1/\Delta k}, with \eqn{\Delta k}
##' the spacing of the taxis knots, and so shifts the predicted position by
##' \eqn{\kappa G T / \Delta k} over a time \eqn{T}, with \eqn{G} the covariate
##' gradient. Against the diffusive spread \eqn{2 D T} per axis, that carries
##' information \eqn{(\kappa G / \Delta k)^2 T / (2 D)}, while the log
##' diffusion carries information of order one per step. Equating the two gives
##'
##' \deqn{\kappa = \Delta k \sqrt{2 D_0 / T} / G,}
##'
##' where \eqn{T} is the median time between successive observations, \eqn{D_0}
##' the diffusion starting value (see below), and \eqn{G} the root mean square
##' gradient at the tag positions (the information adds up \eqn{G^2}, so a few
##' steep-gradient positions count for more than the median suggests). Without
##' taxis knots, \eqn{\Delta k} falls back to the tag-sampled covariate range.
##' Per covariate the values are reduced by their median across the `cov_taxis`
##' covariates. Only the tag types enabled in `conf` contribute, since
##' `dat$tags` may still carry types that `conf` switches off.
##'
##' \eqn{D_0} and \eqn{T} are computed **per tag type** and the resulting
##' \eqn{\kappa_{type}} combined by their geometric mean, rather than pooling all
##' steps into one median. Tag types differ in their step scale by construction
##' -- a mark-recapture tag contributes one net displacement over months, a
##' data-storage tag hundreds of daily increments -- so a median over the pooled
##' steps takes \eqn{T} from whichever type has more steps and \eqn{D_0} from
##' the mixture, and can land outside the range of the per-type values. The
##' geometric mean is the natural average for a multiplicative scale factor and
##' always lies between them.
##'
##' The diffusion coefficients `beta` start at a flat diffusion
##' \eqn{D_0 = \mathrm{median}(d^2 / \Delta t) / (4 \log 2)}, from the
##' displacements \eqn{d} over the times \eqn{\Delta t} between successive
##' observations: under the model \eqn{d^2 / \Delta t} is exponential with mean
##' \eqn{4 D}. Drift and observation error bias it upwards. Starting at
##' \eqn{D = 1} in data units instead can leave the optimizer stuck orders of
##' magnitude away from the fitted diffusion, in particular when coordinates are
##' in km or m.
##'
##' If the quantities above cannot be computed (no covariates, no tags, or a
##' covariate with no spatial variation), the function falls back to the earlier
##' rule `kappa = cellsize^2 / median_dt` and, failing that, to `kappa = 1`.
##'
##' Override via `par$logKappa <- log(<value>)` after calling this function.
##' The value does not change the likelihood, but it does change the path of the
##' optimizer: when the objective has several local optima in `alpha`, fits with
##' different `kappa` can end in different ones. Compare their objective values
##' rather than taking either at face value.
##'
##' Advection has two parameters. `gamma` is an array
##' `[direction (x, y), field, season]` of entrainment coefficients, one set per
##' advection field of the data (see [prep_adv()]); it is absent when the data
##' have no advection field. `adv_const` is a `2 x season` matrix with a constant
##' drift (x, y) in space units per time unit. Both start at 0, the model being
##' linear in them. With `conf$adv_gamma = "shared"` [default_map()] ties the x
##' and y rows of `gamma` together; `adv_const` is only estimated when
##' `conf$adv_const` is `TRUE`.
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

  ## Number of spline slices per covariate. The third dimension of alpha and
  ## beta is shared, so it is sized by the covariate that uses the most slices;
  ## default_map() fixes the slices each covariate never evaluates.
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

  ## Start at a flat diffusion D0 matched to the observed displacements. The
  ## smooths add up across covariates and default_map() frees the first knot of
  ## covariate j0 only, so log(D0) goes on every knot of j0 and nowhere else;
  ## a start at 0 (D = 1 in data units) can leave nlminb stuck at a diffusion
  ## that is orders of magnitude off.
  tags_use <- .get_tags_in_use(dat, conf)
  D0_type <- .diffusion_by_tag_type(tags_use)
  if (length(D0_type) > 0) {
    D0 <- exp(mean(log(D0_type)))
    j0 <- which.max(rep_len(.get_nsea(dat), ncol(knots_dif)))
    par$beta[, j0, ] <- log(D0)
    if (verbose) {
      us <- tryCatch(units_space(dat), error = function(e) NA)
      ut <- tryCatch(units_time(dat), error = function(e) NA)
      unit_txt <- if (length(us) == 1L && length(ut) == 1L &&
                        !is.na(us) && !is.na(ut)) {
        paste0(" ", us, "^2/", ut)
      } else ""
      message("Diffusion started at ", signif(D0, 4), unit_txt,
              " from the observed tag displacements.")
    }
  }

  ## Advection ----------------------------------------
  ## One entrainment coefficient per field and direction (x, y) and a constant
  ## drift, per advection season. Both enter linearly, so 0 is a fine start.
  conf <- .adv_conf(conf, dat)
  nsea_adv <- conf$n_seasons_adv
  if (length(dat$adv) > 0L) {
    par$gamma <- array(0, dim = c(2L, length(dat$adv), nsea_adv))
  }
  par$adv_const <- matrix(0, 2L, nsea_adv)


  ## Taxis scaling -------------------------------------
  ## kappa is a pure scale factor (only kappa * alpha is identifiable), so its
  ## job is to put alpha on an O(1) scale. Balance the information on alpha
  ## against that on the diffusion, kappa = dk * sqrt(2 * D0 / T) / G. Only the
  ## tag types enabled in conf are used, only the cov_taxis covariates enter the
  ## median across covariates, and the tag types are reduced separately and then
  ## combined. See dev/code_notes.org, "Starting values: kappa and diffusion".
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


## kappa = dk * sqrt(2 * D0 / T) / G evaluated separately for each tag type
## present, named by the tag type letter. Pooling the steps of different tag
## types would take T from whichever type contributes the most steps (a handful
## of data-storage tags sampled daily outnumber thousands of mark-recapture
## displacements) and D0 from the mixture of both, giving a value that need not
## lie between the per-type ones.
.kappa_by_tag_type <- function(dat, tags, idx) {

  out <- numeric(0)

  if (is.null(tags) || nrow(tags) == 0 || is.null(dat$cov) ||
        is.null(dat$time_cov) || length(idx) == 0) {
    return(out)
  }

  dk <- .knot_spacing(dat$knots_tax, idx)
  D0_type <- .diffusion_by_tag_type(tags)

  for (ty in names(D0_type)) {

    tt <- tags[as.character(tags$tag_type) == ty, , drop = FALSE]

    steps <- .get_tag_steps(tt)
    tl <- median(steps$dt[steps$dt > 0], na.rm = TRUE)
    if (!is.finite(tl) || tl <= 0) next

    cs_sum <- .cov_scales_at_tags(dat, tt, idx)
    ## no knot spacing (knots not set yet): one unit of alpha over the sampled
    ## covariate range instead
    step_cov <- ifelse(is.finite(dk) & dk > 0, dk, cs_sum$range)
    ki <- step_cov * sqrt(2 * D0_type[[ty]] / tl) / cs_sum$grad_rms
    ki <- ki[is.finite(ki) & ki > 0]
    if (length(ki) > 0) out[ty] <- median(ki)
  }

  out
}


## Mean spacing of the taxis knots per covariate, following `idx`. NA where the
## knots are missing or constant (a covariate without a taxis smooth).
.knot_spacing <- function(knots, idx) {

  out <- rep(NA_real_, length(idx))
  if (is.null(knots)) return(out)
  knots <- as.matrix(knots)

  for (k in seq_along(idx)) {
    if (idx[k] > ncol(knots)) next
    kn <- knots[, idx[k]]
    kn <- kn[is.finite(kn)]
    if (length(kn) < 2) next
    d <- diff(range(kn)) / (length(kn) - 1)
    if (d > 0) out[k] <- d
  }

  out
}


## Flat diffusion matching the displacements, per tag type. Under the KF each
## axis gains variance 2 * D * dt, so dl^2 / dt is exponential with mean 4 * D
## and median 4 * D * log(2). The median keeps a few long, directed steps from
## dominating; drift and observation error still bias D0 upwards, which is the
## safe side for a starting value.
.diffusion_by_tag_type <- function(tags) {

  out <- numeric(0)
  if (is.null(tags) || nrow(tags) == 0) return(out)

  for (ty in unique(as.character(tags$tag_type))) {
    steps <- .get_tag_steps(tags[as.character(tags$tag_type) == ty, , drop = FALSE])
    if (is.null(steps)) next
    ok <- is.finite(steps$dt) & steps$dt > 0 & is.finite(steps$dl)
    if (!any(ok)) next
    d0 <- median(steps$dl[ok]^2 / steps$dt[ok]) / (4 * log(2))
    if (is.finite(d0) && d0 > 0) out[ty] <- d0
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

  ## Rows sharing an event are candidate positions for one ambiguous
  ## observation, not successive steps: differencing across them would invent
  ## zero-length time steps and displacements between alternatives, corrupting
  ## the diffusion and kappa starting values. Keep the most likely candidate of
  ## each event and difference those.
  if (!is.null(tags[["event"]])) {
    keep <- unlist(lapply(split(seq_len(nrow(tags)),
                                paste0(grp, "-", tags$event)),
                          function(k) if (length(k) == 1L) k else
                            k[which.max(.na_zero(tags[["prob"]][k]))]),
                   use.names = FALSE)
    keep <- sort(keep)
    tags <- tags[keep, , drop = FALSE]
    grp <- grp[keep]
  }

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
## median and root mean square gradient magnitude at the tag positions. Uses nearest-cell lookup
## rather than interpolation -- kappa only has to be right to within an order of
## magnitude, and this avoids building interpolators for every time slice.
## `idx` restricts the work to the covariates that carry the taxis; the returned
## vectors follow `idx`.
.cov_scales_at_tags <- function(dat, tags, idx = seq_along(dat$cov)) {

  rng <- rep(NA_real_, length(idx))
  grd <- rep(NA_real_, length(idx))
  grd_rms <- rep(NA_real_, length(idx))

  if (is.null(tags) || nrow(tags) == 0) {
    return(list(range = rng, grad = grd, grad_rms = grd_rms))
  }

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
    it <- .cov_slice(tags$t, dat, i)

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
    if (length(gok) > 0) {
      grd[k] <- median(gok)
      grd_rms[k] <- sqrt(mean(gok^2))
    }
  }

  list(range = rng, grad = grd, grad_rms = grd_rms)
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
