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
##' @param verbose Logical; if \code{TRUE}, informative messages are printed.
##'
##' @details
##' This function generates a named list of initial parameter values based on
##' the supplied data and model configuration. These values can be used as
##' starting values for model fitting and are intended to provide reasonable
##' defaults for the selected model setup.
##'
##' The taxis scaling parameter `logKappa` is fixed during estimation (see
##' [default_map()]) and therefore its initial value is its final value. Because
##' `kappa` has units of \eqn{[\text{distance}^2 / \text{time}]}, a value of 1
##' is only appropriate when coordinates are already on a unit scale. For
##' projected coordinates such as UTM (metres), `kappa = 1` makes the taxis
##' contribution negligible and renders the taxis spline coefficients
##' unidentifiable. The default is therefore set to
##' `kappa = cellsize^2 / median_dt`, which ensures that a unit covariate
##' gradient over one grid-cell width produces movement of one cell width per
##' median time step. Override via `par$logKappa <- log(<value>)` after calling
##' this function.
##'
##' @return
##' A named list of initial parameter values.
##'
##' @examples
##' par <- with(skjepo$sim, default_par(dat, conf))
##'
##' @export
default_par <- function(dat, conf = NULL, verbose = TRUE) {

  if (is.null(conf)) {
    if (verbose) message("No configuration list provided, using default_conf(dat). ")
    conf <- default_conf(dat, verbose = verbose)
  }

  par <- list()

  ## Taxis -----------------------------------------

  if (is.null(dat$knots_tax)) {
    knots_tax <- matrix(NA, 1, 1)
  } else {
    knots_tax <- dat$knots_tax
  }

  max_seasonal <- 1
  if (any(conf$seasonal_spline)) {
    for (i in seq_along(dat$time_spline)) {
      if (isTRUE(conf$seasonal_spline[i])) {
        max_seasonal <- max(max_seasonal, length(dat$time_spline[[i]]))
      }
    }
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

  max_seasonal <- 1
  if (any(conf$seasonal_spline)) {
    for (i in seq_along(dat$time_spline)) {
      if (isTRUE(conf$seasonal_spline[i])) {
        max_seasonal <- max(max_seasonal, length(dat$time_spline[[i]]))
      }
    }
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

  max_seasonal <- 1
  if (any(conf$seasonal_spline)) {
    for (i in seq_along(dat$time_spline)) {
      if (isTRUE(conf$seasonal_spline[i])) {
        max_seasonal <- max(max_seasonal, length(dat$time_spline[[i]]))
      }
    }
  }

  par$gamma <- array(rep(0, length(cov)),
                     dim = c(2,
                             length(cov),
                             max_seasonal))


  ## Taxis scaling -------------------------------------
  ## kappa has units [distance^2 / time]; scale it so that a unit covariate
  ## gradient over one grid-cell width produces movement of one cell width per
  ## median time step: kappa = cellsize^2 / median_dt.
  ## Falls back to log(1) when grid or tag timing are unavailable.
  cs <- if (!is.null(dat$grid)) dat$grid$cellsize[1] else 1
  all_dts <- if (!is.null(dat$tags)) {
    tags_split <- split(dat$tags, dat$tags$id)
    unlist(lapply(tags_split, function(tg) diff(tg$t)))
  } else {
    NULL
  }
  med_dt <- if (length(all_dts) > 0) median(all_dts, na.rm = TRUE) else 1
  if (is.na(med_dt) || med_dt <= 0) med_dt <- 1
  par$logKappa <- log(cs^2 / med_dt)



  ## Observation uncertainty ---------------------------
  ## 3 tag types
  par$logSdO <- matrix(rep(0,6), 2, 3)

  ## return
  par
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
