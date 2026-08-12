
## Main functions ----------------------------------------------------------------------

##' Fit an admove movement model
##'
##' @description
##' Fits the \code{admove} movement model to tagging data to estimate movement
##' processes and, where applicable, habitat preference and advection effects.
##' The model supports different types of tagging data and can be fitted using
##' either a Kalman filter (continuous space, discrete time) or a continuous-time
##' Markov chain (CTMC; discrete space, continuous time) formulation.
##'
##' @param dat A data list containing model input data, as produced by
##'   [setup_data()].
##' @param conf An optional configuration list, typically created by
##'   [default_conf()]. If \code{NULL}, a default configuration is generated
##'   from \code{dat}.
##' @param par An optional named list of initial parameter values, typically
##'   created by [default_par()]. If \code{NULL}, default initial values are
##'   generated from \code{dat} and \code{conf}.
##' @param map An optional parameter map, typically created by [default_map()].
##'   If \code{NULL}, a default map is generated from \code{dat}, \code{conf},
##'   and \code{par}.
##' @param engine Optional integer to override \code{conf$engine}. Use
##'   \code{1} for the Kalman filter and \code{2} for the CTMC formulation.
##'   If \code{NULL}, the value in \code{conf} is used.
##' @param run Logical; if \code{TRUE} (default), the model is optimized. If
##'   \code{FALSE}, only the RTMB objective object is constructed and returned.
##' @param lower Optional lower bounds for optimization. If \code{NULL}, no
##'   explicit lower bounds are supplied.
##' @param upper Optional upper bounds for optimization. If \code{NULL}, no
##'   explicit upper bounds are supplied.
##' @param rel_tol Relative convergence tolerance passed to [stats::nlminb()].
##'   Default is \code{1e-10}.
##' @param do_predictions Logical; if \code{TRUE} (default), model predictions
##'   are computed after fitting. If \code{FALSE}, prediction-related outputs are
##'   skipped, and some plotting methods may not be available.
##' @param do_sdreport Logical; if \code{TRUE} (default), [RTMB::sdreport()] is
##'   run to obtain uncertainty estimates for model parameters and derived
##'   quantities.
##' @param do_report Logical; if \code{TRUE} (default), \code{obj$report()} is
##'   run to extract reported RTMB quantities.
##' @param do_tag_dist Logical; if \code{TRUE}, predicted location distributions
##'   are precomputed for all tags via [add_tag_dist()] and stored in
##'   \code{fit$tag_dist} (consumed by [plot_tag_dist()]). Default is
##'   \code{FALSE}, as this is per-tag expensive and only needed for tag-location
##'   visualisation; it requires a prediction grid and runs after the steps
##'   above. Compute distributions for selected tags later with
##'   \code{add_tag_dist(fit, i = ...)}.
##' @param save_covariance Logical; if \code{TRUE}, the covariance matrix from
##'   [RTMB::sdreport()] is retained. This may substantially increase memory use.
##' @param dbg Logical; if \code{TRUE}, the function is run in debugging mode.
##'   Default is \code{FALSE}.
##' @param control An optional named list of control settings passed to the
##'   optimizer.
##' @param verbose Logical; if \code{TRUE}, progress messages are printed.
##' @param ... Additional arguments passed to [RTMB::MakeADFun()].
##'
##' @details
##' This is the main model-fitting function in \code{admove}. It combines the
##' supplied data, configuration, parameter values, and parameter mapping,
##' constructs an RTMB objective function, and optionally optimizes it using
##' [stats::nlminb()].
##'
##' If configuration settings, initial parameter values, or parameter maps are
##' not supplied, they are generated automatically using [default_conf()],
##' [default_par()], and [default_map()], respectively.
##'
##' @return
##' A fitted model object of class \code{"admove"}. Depending on the function
##' arguments, the returned object may include the RTMB objective function, the
##' optimization output, reported quantities, predictions, and uncertainty
##' estimates.
##'
##' @examples
##' fit <- admove(skjepo$sim, do_sdreport = FALSE)
##'
##' @importFrom RTMB MakeADFun sdreport
##' @importFrom stats nlminb
##'
##' @export
admove <- function(dat,
                   conf = NULL,
                   par = NULL,
                   map = NULL,
                   engine = NULL,
                   run = TRUE,
                   lower = NULL,
                   upper = NULL,
                   rel_tol = 1e-10,
                   do_predictions = TRUE,
                   do_sdreport = TRUE,
                   do_report = TRUE,
                   do_tag_dist = FALSE,
                   save_covariance = FALSE,
                   dbg = FALSE,
                   control = NULL,
                   verbose = TRUE,
                   ...) {

  ## Flags
  sim_flag <- ifelse(inherits(dat, "admove_sim"), TRUE, FALSE)

  ## RTMB does not allow you to pass the dat to MakeADFun
  cmb <- function(f, d) function(p) f(p, d)

  if (sim_flag) {
    sim <- dat
    dat <- sim$dat
    if(is.null(conf)) conf <- sim$conf
    if(!is.null(engine)) conf$engine <- engine
    if(is.null(par)) par <- sim$par
    if(is.null(map)) map <- sim$map
  }

  ## Assume defaults if not provided
  if(is.null(conf)) conf <- default_conf(dat)
  if(!is.null(engine)) conf$engine <- engine

  ## Derive the seasonal spline breakpoints implied by conf$n_seasons before
  ## anything reads them: the parameter dimensions, the map and the likelihood
  ## all index dat$time_spline.
  res_sea <- .resolve_seasons(dat, conf)
  dat <- res_sea$dat
  conf <- res_sea$conf

  if(is.null(par)) par <- default_par(dat, conf)
  if(is.null(map)) map <- default_map(dat, conf, par)

  conf$engine <- .get_engine_integer(conf$engine)
  conf <- .check_seasonal_lengths(conf, dat)

  ## check and clean tags
  dat$tags <- check_tags(dat$tags, dat$grid, dat, conf, TRUE, verbose)

  ## period for seasonality
  dat$period <- period(dat)

  ## check time_spline starts at 0 for seasonal covariates
  if (any(conf$seasonal_spline) && !is.null(dat$time_spline)) {
    for (i in seq_along(conf$seasonal_spline)) {
      if (isTRUE(conf$seasonal_spline[i])) {
        ts_i <- dat$time_spline[[i]]
        if (!is.null(ts_i) && length(ts_i) >= 1L && ts_i[1L] != 0) {
          stop("dat$time_spline[[", i, "]] starts at ", ts_i[1L],
               " but must start at 0 when conf$seasonal_spline[", i, "] is TRUE. ",
               "The seasonal spline uses time modulo the period, so the first ",
               "breakpoint must be 0 to ensure all observations are covered. ",
               "Either set the number of seasons in the configuration instead, ",
               "with conf <- set_seasons(conf, dat, n = ", max(length(ts_i), 2L),
               "), which derives the breakpoints for you, or fix them by hand: ",
               "dat$time_spline[[", i, "]] <- c(0, ", ts_i[-1L], ")",
               call. = FALSE)
        }
      }
    }
  }

  ## check that kappa is not estimated alongside alpha
  .check_kappa_map(par, map, conf)

  ## check that mapping in line with obs_var_type
  ind_t_use <- c(conf$use_dtags, conf$use_stags, conf$use_ctags)
  obs_var_type_map <- !apply(matrix(map$logSdO, 2, 3)[,ind_t_use, drop = FALSE], 2, function(x) any(is.na(x)))
  obs_var_type_conf <- sapply(conf$obs_var_type[ind_t_use], function(x) ifelse(x %in% c(1,2), TRUE, FALSE))
  if(any(obs_var_type_conf != obs_var_type_map)) stop("conf$obs_var_type and mapped parameters (map) do not agree. Did you manipulate map, but not conf? Please check!")

  ## check that sdx and sdy in tags if obs_var_type == 3
  if (any(conf$obs_var_type == 3) &&
        (!any(colnames(dat$tags) == "sdx") ||
           !any(colnames(dat$tags) == "sdy"))) stop("Option to use imported observation uncertainty specified (conf$obs_var_type = 3), but columns 'sdx' and 'sdy' not provided in dat$tags! Please add these columns with the respective information.")

  ## extra checks for CTMC
  if (conf$engine == 2) {
    if (is.null(dat$grid)) stop("No grid provided! CTMC (engine = 2) requires a grid (dat$grid). See create_grid()!")
    if (!any(colnames(dat$tags) == "ic")) stop("Tags are not matched to the grid cells (column tags$ic is missing). Run check_tags()!")

    if (verbose && identical(conf$drift_scheme, "central")) {
      message("Using central-difference drift scheme (conf$drift_scheme). ",
              "If the optimizer fails to converge or predictions warn about ",
              "negative generator rates, the grid is too coarse for the drift ",
              "(grid-Peclet > 2): refine the grid or switch to conf$drift_scheme = \"upwind\".")
    }
  }

  ## Combine conf and dat
  tmb_all <- c(dat, conf)
  tmb_all$tags <- split(dat$tags, dat$tags$id)
  tmb_all$dbg <- dbg

  if(verbose) message("Building the model, that can take a few minutes.")

  t1 <- Sys.time()
  obj <- RTMB::MakeADFun(func = cmb(nll, tmb_all),
                         parameters = par,
                         map = map,
                         silent = TRUE,
                         ...)
  t2 <- Sys.time()


  ## parameter bounds
  if (is.null(lower)) lower <- .get_lower_bounds(par)
  if (is.null(upper)) upper <- .get_upper_bounds(par)
  lower2 <- .get_lower_bounds(par)
  for(nn in names(lower)) lower2[names(lower2) == nn] <- lower[nn]
  lower2 <- .get_non_na_from_map(lower2, map)
  upper2 <- .get_upper_bounds(par)
  for(nn in names(upper)) upper2[names(upper2) == nn] <- upper[nn]
  upper2 <- .get_non_na_from_map(upper2, map)

  if(!run) return(list(sdrep = NA,
                       pl = obj$par,
                       plsd = NA,
                       dat = dat,
                       conf = conf,
                       par = par,
                       map = map,
                       opt = NA,
                       obj = obj))

  if(verbose) message(paste0("Model built (",
                                signif(as.numeric(difftime(t2, t1,
                                                           units = "mins")),2),
                                "min). Minimizing neg. loglik."))

  ## default list
  ctrl <- list(trace = as.integer(verbose),
               eval.max = 2000,
               iter.max = 1000,
               rel.tol = rel_tol)

  if(is.null(control)){
    ind0 <- names(control) %in% names(ctrl)
    ind <- match(names(control), names(ctrl[ind0]))
    ## overwrite
    if(length(ind) > 0){
      ctrl[ind] <- ctrl[ind0]
    }
    ## add
    if(length(which(!ind0)) > 0){
      ctrl <- c(ctrl,
                control[!(names(control) %in% names(ctrl[!ind0]))])
    }
  }

  opt <- stats::nlminb(obj$par, obj$fn, obj$gr,
                       control = ctrl,
                       lower = lower2,
                       upper = upper2)
  t3 <- Sys.time()

  if(verbose) message(paste0("Minimisation done (",
                                signif(as.numeric(difftime(t3, t2,
                                                           units = "mins")),2),
                                "min). Model ", "not "[opt$convergence],
                                "converged."))

  res <- list(dat = dat,
              conf = conf,
              par = par,
              map = map,
              opt = opt,
              obj = obj,
              low = lower,
              hig = upper)

  res$times <- c(makeadfun = signif(as.numeric(difftime(t2, t1, units = "mins")),2),
                 nlminb = signif(as.numeric(difftime(t3, t2, units = "mins")),2))

  attr(res, "RemoteSha") <- substr(packageDescription("admove")$RemoteSha, 1, 12)
  attr(res, "Version") <- packageDescription("admove")$Version
  res <- .add_class(res, "admove")
  res <- add_sref(res, sref(res$dat))
  res <- add_tref(res, tref(res$dat))


  if (do_predictions) {

    if (is.null(dat$pred$grid$igrid)) {

      if(verbose) message(paste0("No prediction grid provided; skipping predictions."))

    } else {

      if (verbose) message(paste0("Predicting movement rates."))

      res <- add_predictions(res)

      if(verbose) message(paste0("Predictions done (",
                                    res$times[which(names(res$times) == "predictions")],
                                    "min)."))

    }
  }


  if(do_sdreport){

    if(verbose) message(paste0("Estimating uncertainty."))

    res <- add_sdreport(res, save_covariance)

    if(verbose) message(paste0("SDreporting done (",
                                  res$times[which(names(res$times) == "sdreport")],
                                  "min)."))

  } else {
    res$pl <- .get_pl_from_opt(par, map, opt)
  }

  if(do_report){

    if(verbose) message(paste0("Reporting variables."))

    res <- add_report(res)

    if(verbose) message(paste0("Reporting done (",
                                  res$times[which(names(res$times) == "report")],
                                  "min)."))
  }

  if (do_tag_dist) {

    if (is.null(dat$pred$grid$igrid)) {

      if (verbose) message("No prediction grid provided; skipping tag distributions.")

    } else {

      if (verbose) message("Computing tag location distributions.")

      res <- add_tag_dist(res)

    }
  }

  return(res)
}


##' Add an RTMB sdreport to a fitted admove model
##'
##' @description
##' Runs [RTMB::sdreport()] for a fitted \code{admove} model and adds the result
##' to the fitted object. Estimated values and standard deviations are also
##' extracted as named lists and stored in the returned object.
##'
##' @param fit A fitted model object of class \code{"admove"}, as returned by
##'   [admove()].
##' @param save_covariance Logical; if \code{TRUE}, the full covariance matrix
##'   from the sdreport is retained. If \code{FALSE} (default), the covariance
##'   matrix is removed to reduce memory use.
##' @param ad_hessian Logical; if \code{TRUE} (default), the exact AD Hessian
##'   from \code{obj$he()} is supplied to [RTMB::sdreport()] instead of letting
##'   it difference the gradient numerically. See Details. Set to \code{FALSE}
##'   to restore the [RTMB::sdreport()] default behaviour.
##'
##' @details
##' This function adds three components to the fitted object:
##'
##' \itemize{
##'   \item \code{sdrep}: the full sdreport object returned by [RTMB::sdreport()]
##'   \item \code{pl}: a named list of estimates extracted with \code{as.list(sdrep, "Est")}
##'   \item \code{plsd}: a named list of standard deviations extracted with \code{as.list(sdrep, "Std")}
##' }
##'
##' If these components already exist in \code{fit}, they are overwritten.
##'
##' For models without random effects -- which is every \code{admove} model --
##' [RTMB::sdreport()] builds the Hessian with
##' \code{stats::optimHess(par, obj$fn, obj$gr)}, i.e. by differencing the
##' gradient with a fixed step (\code{ndeps}, default \code{1e-3}). That is both
##' inaccurate and fragile here:
##'
##' \itemize{
##'   \item the step is absolute, so it is far too coarse whenever the
##'     parameters sit on a small scale, and the resulting standard errors can be
##'     wrong by factors of several;
##'   \item if the likelihood is not finite \code{ndeps} away from the optimum --
##'     which happens readily when a predicted track leaves the covariate domain
##'     or crosses a masked (\code{NA}) cell -- the whole Hessian becomes
##'     \code{NaN}, and with it every standard error, including those of the
##'     \code{ADREPORT}ed quantities.
##' }
##'
##' \code{obj$he()} instead returns the exact second derivatives from the AD
##' tape, at no extra cost (the runtime of [RTMB::sdreport()] is dominated by the
##' delta-method step for the \code{ADREPORT}ed quantities). It does not mask
##' genuine problems: an indefinite Hessian is still reported via
##' \code{sdrep$pdHess}, and if the objective really is \code{NaN} at the
##' optimum, \code{obj$he()} returns \code{NaN} too. If \code{obj$he()} is
##' unavailable or fails, the function silently falls back to the
##' [RTMB::sdreport()] default.
##'
##' @return
##' An updated object of class \code{"admove"} with sdreport results added.
##'
##' @export
add_sdreport <- function(fit, save_covariance = FALSE, ad_hessian = TRUE) {

  .check_class(fit, "admove")

  res <- fit

  hess <- if (isTRUE(ad_hessian)) .get_ad_hessian(fit$obj) else NULL

  t1 <- Sys.time()
  if (is.null(hess)) {
    sdrep <- RTMB::sdreport(obj = fit$obj)
  } else {
    sdrep <- RTMB::sdreport(obj = fit$obj, hessian.fixed = hess)
  }
  t2 <- Sys.time()

  pl <- as.list(sdrep, "Est")
  plsd <- as.list(sdrep, "Std")

  if (!save_covariance) {
    sdrep$cov <- NULL ## save memory
  }

  ## overwrite
  if (any(names(res) == "sdrep")) {
    res$sdrep <- sdrep
  } else {
    res <- c(res, list(sdrep = sdrep))
  }
  if (any(names(res) == "pl")) {
    res$pl <- pl
  } else {
    res <- c(res, list(pl = pl))
  }
  if (any(names(res) == "plsd")) {
    res$plsd <- plsd
  } else {
    res <- c(res, list(plsd = plsd))
  }

  res$times <- c(res$times,
                 sdreport = signif(as.numeric(difftime(t2, t1, units = "mins")),2))


  res <- .add_class(res, "admove")
  res <- add_sref(res, sref(res$dat))
  res <- add_tref(res, tref(res$dat))

  return(res)
}


##' Add reported RTMB quantities to a fitted admove model
##'
##' @description
##' Runs \code{obj$report()} for a fitted \code{admove} model and adds the
##' reported quantities to the fitted object.
##'
##' @param fit A fitted model object of class \code{"admove"}, as returned by
##'   [admove()].
##'
##' @details
##' The reported quantities are stored in the \code{rep} component of the
##' returned object. If this component already exists, it is overwritten.
##'
##' @return
##' An updated object of class \code{"admove"} with reported quantities added.
##'
##' @export
add_report <- function(fit) {

  .check_class(fit, "admove")

  t1 <- Sys.time()
  rep <- fit$obj$report()
  t2 <- Sys.time()
  fit$rep <- rep

  fit$times <- c(fit$times,
                 report = signif(as.numeric(difftime(t2, t1, units = "mins")),2))

  fit <- .add_class(fit, "admove")
  return(fit)
}



##' Add model predictions to a fitted admove model
##'
##' @description
##' Computes and adds model predictions for a fitted \code{admove} model.
##'
##' @param fit A fitted model object of class \code{"admove"}, as returned by
##'   [admove()].
##' @param grid Optional \code{"admove_grid"} (from [create_grid()]) giving a new
##'   spatial grid to predict on. If \code{NULL} (default), the grid stored in
##'   \code{fit$dat$pred$grid} is used. The grid must be in the same spatial
##'   reference as the fit (no reprojection is done) and, when covariates are
##'   present, should lie within their spatial coverage; a warning is issued
##'   otherwise, as habitat fields extrapolate to \code{NA} beyond it.
##' @param time Optional strictly increasing numeric vector of prediction times
##'   on the model-time scale (see [tref()]). If \code{NULL} (default), the times
##'   in \code{fit$dat$pred$time} are used.
##'
##' @details
##' This function evaluates the model prediction step and stores the resulting
##' predicted quantities in \code{fit$pred}, including the habitat, diffusion,
##' taxis, and advection fields on the prediction grid.
##'
##' Predictions are a pure recomputation from the fitted parameters and
##' covariates evaluated on \code{fit$dat$pred$grid} at \code{fit$dat$pred$time}.
##' Supplying \code{grid} and/or \code{time} re-targets the prediction (updating
##' \code{fit$dat$pred$grid} / \code{fit$dat$pred$time} in the returned object)
##' without refitting, so a single fit can be predicted onto different grids or
##' time sequences. Because habitat is obtained by local interpolation of the
##' fitted covariate fields, predictions are only valid within the covariate
##' spatial and temporal coverage.
##'
##' It also stores \code{fit$pred$mstar}, a list of length \code{nt} holding the
##' continuous-time Markov chain (CTMC) generator matrices, one sparse
##' (\code{"dgCMatrix"}) \code{nc x nc} matrix per prediction time slice (see
##' [calc_mstar()]). Off-diagonal entries are instantaneous cell-to-cell
##' movement rates in units of \strong{1 / time} (the reciprocal of the model
##' time unit); diagonal entries are the negative row sums. The generator does
##' not include a time step, so movement (transition) probabilities over a step
##' \code{dt} are obtained by exponentiating it, e.g.
##' \code{Matrix::expm(fit$pred$mstar[[t]] * dt)}.
##'
##' @return
##' An updated object of class \code{"admove"} with model predictions added in
##' \code{fit$pred}, including the CTMC generator list \code{fit$pred$mstar}
##' (units 1 / time).
##'
##' @export
add_predictions <- function(fit, grid = NULL, time = NULL) {

  .check_class(fit, "admove")

  res <- fit
  dat <- fit$dat
  conf <- fit$conf

  ## optionally re-target the prediction grid / time before computing, and
  ## persist the new target in the returned object
  if (!is.null(grid) || !is.null(time)) {
    dat <- .set_pred_target(dat, grid = grid, time = time)
    res$dat <- dat
  }


  ## dimensions
  ncov <- length(dat$cov)
  nt <- length(dat$time)
  nc <- nrow(dat$grid$igrid)
  ntp <- length(dat$pred$time)
  ncp <- nrow(dat$pred$grid$igrid)


  t1 <- Sys.time()
  par_est <- get_par_est(fit$par, fit$map, fit$opt)
  kappa <- exp(par_est$logKappa)


  ## Preference functions, local interpolation and habi objects --------------
  hb <- .build_habi(dat, conf, par_est, period(fit))
  pref_funcs <- hb$pref_funcs
  habi_tax <- hb$habi$tax
  habi_dif <- hb$habi$dif
  habi_adv_x <- hb$habi$adv_x
  habi_adv_y <- hb$habi$adv_y


  hT_pred <- hTdx_pred <- hTdy_pred <- hD_pred <-
    hAx_pred <- hAy_pred <- matrix(0, ncp, ntp)
  for (t in 1:ntp) {
    hD_pred[,t] <- habi_dif$val(dat$pred$grid$xygrid,
                                dat$pred$time[t])
    hT_pred[,t] <- habi_tax$val(dat$pred$grid$xygrid,
                                dat$pred$time[t])
    tmp <- habi_tax$grad(dat$pred$grid$xygrid,
                         dat$pred$time[t])
    hTdx_pred[,t] <- kappa * tmp[,1]
    hTdy_pred[,t] <- kappa * tmp[,2]
    hAx_pred[,t] <- habi_adv_x$val(dat$pred$grid$xygrid,
                                   dat$pred$time[t])
    hAy_pred[,t] <- habi_adv_y$val(dat$pred$grid$xygrid,
                                   dat$pred$time[t])
  }
  t2 <- Sys.time()

  ## habitat fields are NA where the prediction grid/time falls outside the
  ## fitted covariate coverage (interpolation cannot extrapolate); flag it
  na_cells <- which(!stats::complete.cases(
    cbind(as.vector(hD_pred), as.vector(hTdx_pred), as.vector(hAx_pred))))
  if (length(na_cells) > 0L) {
    warning("Predictions contain ", length(na_cells), " NA value(s): the ",
            "prediction grid/time extends beyond the fitted covariate coverage, ",
            "so habitat fields (and the CTMC generator) are undefined there. ",
            "Restrict 'grid'/'time' to the covariate domain.", call. = FALSE)
  }


  pred <- list()
  pred$pref_funcs <- pref_funcs
  pred$habi <- list(tax = habi_tax,
                    dif = habi_dif,
                    adv_x = habi_adv_x,
                    adv_y = habi_adv_y)
  pred$hTdx <- hTdx_pred
  pred$hTdy <- hTdy_pred
  pred$hD <- hD_pred
  pred$hAx <- hAx_pred
  pred$hAy <- hAy_pred


  ## overwrite
  if (any(names(res) == "pred")) {
    res$pred <- pred
  } else {
    res <- c(res, list(pred = pred))
  }

  res$pred$mstar <- calc_mstar(res)

  res$times <- c(res$times,
                 predictions = signif(as.numeric(difftime(t2, t1, units = "mins")),2))

  res <- .add_class(res, "admove")
  res <- add_sref(res, sref(res$dat))
  res <- add_tref(res, tref(res$dat))

  return(res)
}


## Validate a user-supplied prediction grid / time and write it into dat$pred.
## Hard errors for structurally invalid input; warnings for likely-invalid input
## (sref mismatch, extrapolation beyond covariate coverage).
.set_pred_target <- function(dat, grid = NULL, time = NULL) {

  if (!is.null(grid)) {

    if (!inherits(grid, "admove_grid"))
      stop("'grid' must be an 'admove_grid' object created by create_grid().",
           call. = FALSE)
    if (is.null(grid$igrid) || is.null(grid$xygrid))
      stop("'grid' is missing required fields (igrid/xygrid); rebuild it with ",
           "create_grid().", call. = FALSE)

    ## spatial reference must match the fitted data (admove does not reproject)
    sref_new <- tryCatch(sref(grid), error = function(e) NULL)
    sref_dat <- tryCatch(sref(dat),  error = function(e) NULL)
    if (!is.null(sref_new) && !is.null(sref_dat) &&
          !sref_equal(sref_new, sref_dat)) {
      warning("The spatial reference of 'grid' differs from the fitted data. ",
              "admove does not reproject; rebuild the grid in the fit's sref.",
              call. = FALSE)
    }

    ## drop cells where the fitted covariates are undefined (outside coverage or
    ## masked), exactly as setup_data() prunes the fitting grid, so the resulting
    ## prediction grid is valid everywhere and the CTMC generator is well-defined
    if (!is.null(dat$cov) && !is.null(dat$xrange_cov)) {
      n0 <- nrow(grid$xygrid)
      pr <- .prune_grid_to_cov(grid, dat$cov, dat$xrange_cov, dat$yrange_cov)
      if (length(pr$removed) >= n0)
        stop("None of the cells in 'grid' fall within the fitted covariate ",
             "coverage; nothing to predict on.", call. = FALSE)
      if (length(pr$removed) > 0L)
        message(length(pr$removed), " of ", n0, " prediction grid cell(s) ",
                "removed because the fitted covariates are undefined there ",
                "(outside coverage or masked).")
      grid <- pr$grid
    }

    dat$pred$grid <- grid
  }

  if (!is.null(time)) {

    if (!is.numeric(time) || length(time) < 2L)
      stop("'time' must be a numeric vector of length >= 2 on the model-time ",
           "scale.", call. = FALSE)
    if (any(diff(time) <= 0))
      stop("'time' must be strictly increasing.", call. = FALSE)

    ## Out-of-range times are not pre-flagged: the habitat interpolation clamps
    ## rather than always returning NA, so a heuristic here would false-alarm.
    ## Genuine NA (e.g. undefined habitat) is caught by the post-hoc check in
    ## add_predictions().
    dat$pred$time <- time
  }

  dat
}


## Remove grid cells where any fitted covariate interpolates to NA (outside the
## covariate extent or over masked regions). Mirrors the pruning setup_data()
## applies to the fitting grid, and rebuilds celltable so neighbour lookups and
## cell indexing stay consistent. Returns the pruned grid and removed cell ids.
.prune_grid_to_cov <- function(grid, cov, xrange_cov, yrange_cov) {

  err <- NULL
  for (i in seq_along(cov)) {
    covi <- cov[[i]]
    for (j in seq_len(dim(covi)[3L])) {
      liv <- RTMB::interpol2Dfun(covi[, , j],
                                 xlim = round(xrange_cov[i, ], 5),
                                 ylim = round(yrange_cov[i, ], 5),
                                 R = 1)
      tmp <- liv(round(grid$xygrid[, 1L], 5), round(grid$xygrid[, 2L], 5))
      ind <- which(is.na(tmp))
      if (length(ind) > 0L) err <- c(err, ind)
    }
  }
  err <- sort(unique(err))

  if (length(err) == 0L) return(list(grid = grid, removed = integer(0)))

  ind <- match(err, grid$celltable)
  grid$celltable[ind] <- NA
  grid$celltable[!is.na(grid$celltable)] <- seq_len(sum(!is.na(grid$celltable)))
  grid$xygrid <- grid$xygrid[-err, , drop = FALSE]
  grid$igrid <- grid$igrid[-err, , drop = FALSE]

  list(grid = grid, removed = err)
}


##' Compute predicted location distributions for one or more tags
##'
##' @description
##' Propagates the model's predicted spatial location distribution for one or
##' more archival tags from release to recovery time, using either the CTMC
##' forward-pass (engine 2) or a Kalman-filter forward pass (engine 1). Results
##' are stored in `fit$tag_dist` (a named list keyed by tag index) and consumed
##' by [plot_tag_dist()]. Separating the expensive computation from rendering
##' means plot aesthetics can be changed without re-running the model.
##'
##' Successive calls accumulate into the same list, so distributions can be
##' added in batches without discarding earlier results.
##'
##' Tags recaptured within a single time step (engine 1) or with fewer than two
##' position records are skipped with a warning rather than stopping.
##'
##' @param fit A fitted object of class `admove`, as returned by [admove()].
##' @param i Integer index or vector of indices of the tags to process. `NULL`
##'   (default) processes all tags.
##' @param dt Time step used for the Kalman-filter forward pass (engine 1 only).
##'   Default is `0.5`.
##' @param engine Optional integer overriding the engine stored in
##'   `fit$conf$engine`. `1` = Kalman filter, `2` = CTMC.
##' @param xrel0,yrel0 Optional coordinates overriding the release location
##'   for CTMC-based predictions (applied to every tag in `i`).
##'
##' @return
##' A copy of `fit` with `$tag_dist` set to a named list (one entry per
##' successfully processed tag index) containing the precomputed distributions
##' and metadata required by [plot_tag_dist()].
##'
##' @seealso [plot_tag_dist()]
##'
##' @export
add_tag_dist <- function(fit, i = NULL, dt = 0.5,
                         engine = NULL, xrel0 = NULL, yrel0 = NULL) {

  .check_class(fit, "admove")

  dat <- fit$dat
  conf <- fit$conf

  if (inherits(dat$tags, "data.frame")) {
    tags <- split(dat$tags, dat$tags$id)
  } else {
    tags <- dat$tags
  }

  if (is.null(i)) i <- seq_along(tags)

  if (any(i < 1L | i > length(tags)))
    stop("Tag index 'i' contains values outside [1, ", length(tags), "].")

  if (is.null(engine)) engine <- conf$engine

  ## accumulate into existing list; convert old single-entry format if needed
  tag_dist_list <- fit$tag_dist
  if (is.null(tag_dist_list)) {
    tag_dist_list <- list()
  } else if (!is.null(tag_dist_list$engine)) {
    old_i <- tag_dist_list$i
    tag_dist_list <- setNames(list(tag_dist_list), as.character(old_i))
  }

  skipped <- integer(0L)

  ## CTMC generator list: reuse the one add_predictions() already stored, and
  ## build it only once (it does not depend on the tag), falling back to a fresh
  ## computation only if predictions were not run
  mstar <- if (engine == 2L) {
    if (!is.null(fit$pred$mstar)) fit$pred$mstar else calc_mstar(fit)
  } else NULL

  for (idx in i) {

    tag <- tags[[idx]]
    ind <- which(apply(!is.na(tag[, 1:3]), 1, all))

    if (length(ind) < 2L) {
      warning("Tag ", idx, " has fewer than 2 non-missing (t, x, y) observations; skipping.",
              call. = FALSE)
      skipped <- c(skipped, idx)
      next
    }

    tag <- tag[ind, ]
    xrel <- if (!is.null(xrel0)) xrel0 else tag[1, 2]
    yrel <- if (!is.null(yrel0)) yrel0 else tag[1, 3]
    trel <- tag[1, 1]
    trec <- tag[nrow(tag), 1]

    if (engine == 2L) {

      tall <- dat$pred$time

      itrel <- as.integer(cut(trel, tall, include.lowest = TRUE))
      itrec <- as.integer(cut(trec, tall, include.lowest = TRUE))
      tind <- sapply(tag$t, function(tt) which.min(abs(tall[itrel:itrec] - tt)))
      nt <- max(itrec) - itrel + 1L

      dist_prob <- matrix(0, nt + 1L, nrow(dat$pred$grid$xygrid))
      icrel <- dat$pred$grid$celltable[cbind(cut(xrel, dat$pred$grid$xgr),
                                             cut(yrel, dat$pred$grid$ygr))]
      dist_prob[1L, icrel] <- 1

      for (k in seq_len(nt)) {
        m <- as.matrix(Matrix::expm(
          mstar[[itrel + k - 1L]] * diff(dat$pred$time)[itrel + k - 1L]))
        dist_prob[k + 1L, ] <- as.vector(dist_prob[k, ] %*% m)
      }

      dens_list <- vector("list", nrow(tag))
      for (k in seq_along(tind)) {
        ct <- dat$pred$grid$celltable
        ct[which(!is.na(ct))] <- dist_prob[tind[k], ]
        dens_list[[k]] <- ct
      }

      tag_dist_list[[as.character(idx)]] <- list(
        engine = engine,
        tag = tag,
        i = idx,
        dens_list = dens_list,
        xg = x_centers(dat$pred$grid),
        yg = y_centers(dat$pred$grid),
        xrange = dat$grid$xrange,
        yrange = dat$grid$yrange
      )

    } else {

      funcs <- default_sim_funcs(dat, conf, fit$pl)
      dt_min <- min(dat$min_dt, median(diff(sort(tag$t))))

      out <- build_time(tag$t, mode = "fixed_dt",
                        dt_min = dt_min, dt = dt, eps = 1)
      ts <- out$ts
      dts <- out$dts
      nts <- out$nts
      observed <- out$observed

      if (nts == 1L) {
        warning("Tag ", idx, " is recaptured within the first time step; skipping. ",
                "Use a smaller dt or pick a different tag.", call. = FALSE)
        skipped <- c(skipped, idx)
        next
      }

      kappa <- exp(fit$pl$logKappa)

      traj <- matrix(NA_real_, nts, 4L)
      colnames(traj) <- c("x0", "y0", "v1", "v2")
      xy0 <- matrix(c(xrel, yrel), 1L, 2L)
      P <- c(0, 0)
      traj[1L, ] <- c(xy0[1L], xy0[2L], 0, 0)

      for (t in 2:nts) {
        dt_t <- dts[t - 1L]
        moveT0 <- if (conf$use_taxis)     kappa * funcs$tax(xy0, ts[t - 1L]) * dt_t else c(0, 0)
        moveA0 <- if (conf$use_advection) funcs$adv(xy0, ts[t - 1L]) * dt_t         else c(0, 0)
        D0 <- exp(funcs$dif(xy0, ts[t - 1L]))

        xy0 <- xy0 + moveT0 + moveA0
        PP <- P + 2 * D0 * dt_t

        if (t %in% observed) {
          ind_obs <- which(observed == t) + 1L
          for (j in seq_along(ind_obs)) {
            F <- PP
            P <- PP - PP / F * PP
            obs_xy <- c(tag$x[ind_obs[j]], tag$y[ind_obs[j]])
            xy0 <- xy0 + PP / F * (obs_xy - xy0)
          }
        } else {
          P <- F <- PP
        }

        traj[t, ] <- c(xy0[, 1L], xy0[, 2L], F[1L], F[2L])
      }

      ind.track <- sapply(tag[, 1L], function(tt) which.min(abs(ts - tt)))
      xrange <- range(traj[, 1L], tag[, 2L], na.rm = TRUE)
      yrange <- range(traj[, 2L], tag[, 3L], na.rm = TRUE)

      tag_dist_list[[as.character(idx)]] <- list(
        engine = engine,
        tag = tag,
        i = idx,
        traj = traj,
        ind.track = ind.track,
        xrange = xrange,
        yrange = yrange
      )
    }
  }

  if (length(skipped) > 0L)
    message(length(skipped), " tag(s) skipped: ", .format_ids(skipped))

  fit$tag_dist <- tag_dist_list
  fit
}


##' Summarise a fitted admove model
##'
##' @description
##' Summarises the main results of a fitted \code{admove} model, including
##' parameter estimates and, where available, associated uncertainty measures.
##'
##' @param object A fitted model object of class \code{"admove"}, as returned by
##'   [admove()].
##' @param CI Numeric scalar giving the confidence level used for confidence
##'   intervals. Default is \code{0.95}.
##' @param ... Additional arguments passed to internal summary methods.
##'
##' @details
##' The parameter table has one row per estimated parameter. Elements that are
##' held equal by the parameter map (e.g. the \eqn{x}- and \eqn{y}-direction
##' advection coefficients, which [default_map()] couples) are a single
##' estimated parameter and are therefore shown once, with a row name listing
##' all element indices they cover, such as \code{gamma3,6}. Fixed elements
##' (\code{NA} in the map) are omitted, so the number of rows matches the model
##' degrees of freedom reported by [logLik()].
##'
##' The taxis scaling parameter \code{kappa} is fixed rather than estimated (see
##' [default_par()]) and is reported separately above the table.
##'
##' @return
##' A summary object, typically printed for inspection.
##'
##' @name summarise_fit
##' @export
summarise_fit <- function(object, CI = 0.95, ...) {
  x <- object

  .check_class(x, "admove")

  if (!"digits" %in% names(list(...))) digits <- 7
  ndigits <- digits # Present values with this number of digits after the dot.

  if(CI > 1 || CI < 0) stop("CI has to be between 0 and 1!")
  zscore <- qnorm(CI + (1 - CI)/2)

  cat("<admove>\n")

  cat(paste(' Convergence: ', x$opt$convergence,
            '  MSG: ', x$opt$message, '\n', sep=''))
  if (x$opt$convergence > 0) {
    cat('WARNING: Model did not obtain proper convergence! Estimates and uncertainties are most likely invalid and cannot be trusted.\n')
  }

  ## if('sderr' %in% names(x)) cat('WARNING: Could not calculate all standard deviations. The optimum found may be invalid. Proceed with caution.\n')
  if (x$opt$convergence > 0) {
    txtobj <- 'Objective function: '
  } else {
    txtobj <- 'Objective function at optimum: '
  }
  cat(paste0(" ", txtobj, round(x$obj$fn(), ndigits), '\n'))

  tags <- x$dat$tags
  if (inherits(tags, "list")) {
    tags <- do.call(rbind, tags)
  }
  dims <- dim(tags)

  tags_split <- split(tags, tags$tag_type)
  tags_split2 <- lapply(tags_split, function(x) split(x, x$id))
  n_tags <- lapply(tags_split2, function(x) ifelse(length(x) == 0, NA, length(x)))

  labw <- 22

  cat("\n")

  for (i in 1:3) {
    if (is.na(n_tags[[i]]) || n_tags[[i]] == 0) next
    cat(sprintf(paste0("  %-", labw, "s %s\n"), paste0(c("data-storage","mark-resight","mark-recapture")[i]," tags:"),
                n_tags[[i]]))
  }

  cat("\n")

  ## Seasonality
  ss <- x$conf$seasonal_spline
  sc <- x$conf$seasonal_cov
  if (isTRUE(any(ss, na.rm = TRUE)) || isTRUE(any(sc, na.rm = TRUE))) {

    per  <- tryCatch(period(x), error = function(e) NA_real_)
    unit <- tryCatch(tref(x)$units, error = function(e) NA_character_)

    per_txt <- if (length(per) == 1L && !is.na(per) && is.finite(per)) {
      paste0(signif(per, ndigits),
             if (!is.null(unit) && !is.na(unit) && nzchar(unit)) paste0(" ", unit) else "")
    } else {
      "unspecified"
    }
    cat(sprintf(paste0("  %-", labw, "s %s\n"), "seasonal period:", per_txt))

    ## conf$n_seasons is the specification; fall back to the array dimensions for
    ## models whose breakpoints were set by hand
    nsea_of <- function(a) if (!is.null(a) && length(dim(a)) >= 3L) dim(a)[3L] else NA_integer_
    if (!is.null(x$conf$n_seasons) && any(x$conf$n_seasons > 1L)) {
      cn <- names(x$dat$cov)
      nsea_vals <- x$conf$n_seasons
      names(nsea_vals) <- if (!is.null(cn) && length(cn) == length(nsea_vals)) {
        cn
      } else {
        paste0("cov", seq_along(nsea_vals))
      }
    } else {
      nsea_vals <- c(taxis = nsea_of(x$par$alpha),
                     diffusion = nsea_of(x$par$beta),
                     advection = nsea_of(x$par$gamma))
    }
    nsea_vals <- nsea_vals[!is.na(nsea_vals)]
    if (length(nsea_vals) > 0L) {
      nsea_txt <- if (length(unique(nsea_vals)) <= 1L) {
        as.character(nsea_vals[1L])
      } else {
        paste(paste0(names(nsea_vals), "=", nsea_vals), collapse = ", ")
      }
      cat(sprintf(paste0("  %-", labw, "s %s\n"), "number of seasons:", nsea_txt))
    }

    ## beta carries a seasonal dimension whenever a seasonal basis is used, but
    ## its coefficients are coupled across seasons unless conf$seasonal_dif is
    ## set, so report what is actually estimated rather than the array dimension.
    cat(sprintf(paste0("  %-", labw, "s %s\n"), "seasonal diffusion:",
                if (isTRUE(x$conf$seasonal_dif)) "yes" else "no (coupled)"))

    cov_names <- names(x$dat$cov)
    lab_cov <- function(flag) {
      wi <- which(isTRUE(flag) | flag)
      if (length(wi) == 0L) return("none")
      if (!is.null(cov_names) && length(cov_names) >= max(wi)) {
        paste(cov_names[wi], collapse = ", ")
      } else {
        paste(wi, collapse = ", ")
      }
    }
    if (isTRUE(any(ss, na.rm = TRUE)))
      cat(sprintf(paste0("  %-", labw, "s %s\n"), "seasonal spline (cov):", lab_cov(ss)))
    if (isTRUE(any(sc, na.rm = TRUE)))
      cat(sprintf(paste0("  %-", labw, "s %s\n"), "seasonal cov:", lab_cov(sc)))

    cat("\n")
  }

  ## kappa is a fixed scale factor (only kappa * alpha is identifiable), so it
  ## is never part of the estimates below -- report it separately so that it is
  ## visible which scale the taxis parameters refer to.
  kappa_fixed <- !is.null(x$map$logKappa) && all(is.na(x$map$logKappa))
  if (kappa_fixed && !is.null(x$pl$logKappa) && isTRUE(x$conf$use_taxis)) {

    us <- tryCatch(units_space(x$dat), error = function(e) NA_character_)
    ut <- tryCatch(units_time(x$dat), error = function(e) NA_character_)

    kappa_unit <- if (length(us) == 1L && !is.na(us) && nzchar(us) &&
                        length(ut) == 1L && !is.na(ut) && nzchar(ut)) {
      paste0(" ", us, "^2/", ut)
    } else {
      ""
    }

    cat(sprintf(paste0("  %-", labw, "s %s\n"), "kappa (fixed scale):",
                paste0(signif(exp(x$pl$logKappa), ndigits), kappa_unit)))
    cat("\n")
  }

  sel <- .select_estimated_par(x$pl, x$map)
  keep <- sel$keep

  pl <- unlist(x$pl, use.names = FALSE)[keep]
  names(pl) <- sel$labels[keep]

  plsd_all <- unlist(x$plsd, use.names = FALSE)
  plsd <- if (length(plsd_all) == length(keep)) plsd_all[keep] else numeric(0)

  pllow <- pl - zscore * plsd
  plup <- pl + zscore * plsd

  cat(' Model parameter estimates w ', CI * 100, '% CI \n', sep = "")

  res <- round(cbind(pl, pllow, plup, plsd), ndigits)
  rownames(res) <- names(pl)
  if (length(pllow) > 0) {
    colnames(res) <- c("estimate","cilow","ciupp","sd")
  } else {
    colnames(res) <- c("estimate")
  }

  cat('',paste(capture.output(res), '\n'), '\n')

  invisible(NULL)
}


##' Summary plots for a fitted `admove` object
##'
##' @description
##' `plot_fit()` creates one or several diagnostic summary plots for a fitted
##' object of class `admove`. Depending on the selected `quantity`, the function
##' visualises habitat preference functions, taxis, diffusion, or parameter
##' estimates. Multiple quantities are arranged automatically in a multi-panel
##' layout.
##'
##' @param x A fitted object of class `admove`, as returned by [admove()].
##' @param quantity Character vector specifying which quantities to plot.
##'   Available options are:
##'   \describe{
##'     \item{`"pref"`}{Taxis habitat-preference function vs the covariate(s).}
##'     \item{`"taxis"`}{Taxis (movement direction and magnitude) in space.}
##'     \item{`"advection"`}{Advection (direction and magnitude) in space, one
##'       panel per season. Skipped automatically when the model carries no
##'       advection (fitted with `conf$use_advection = FALSE`, or every `gamma`
##'       coefficient zero), since the whole field would be zero.}
##'     \item{`"dif"`}{Diffusion in space.}
##'     \item{`"pref_dif"`}{Diffusion as a function of the covariate(s). Only
##'       informative when diffusion has more than one knot (`nknots_dif > 1`);
##'       with a single knot diffusion is covariate-independent and this panel is
##'       skipped automatically.}
##'     \item{`"par"`}{Estimated model parameters.}
##'   }
##'   Multiple quantities can be selected. The default includes `"pref_dif"`, so
##'   a covariate-dependent diffusion spline is shown automatically when present.
##'
##'   For the covariate-preference quantities (`"pref"` and `"pref_dif"`), a
##'   covariate whose spline coefficients are all fixed in the parameter map (all
##'   `NA`, i.e. nothing estimated for it) has no fitted preference relationship
##'   and its panel is omitted. This is the case, for example, for the current
##'   covariates when advection is used, whose taxis (`alpha`) and diffusion
##'   (`beta`) coefficients are mapped off. If every covariate is fixed, the
##'   corresponding quantity contributes no panels.
##' @param plot_land Logical; if `TRUE`, land masses are added to spatial plots
##'   using [plot_land()]. Default: `FALSE`.
##' @param auto_layout Logical; if `TRUE`, graphical parameters are set and
##'   restored automatically, and plots are arranged in a multi-panel layout.
##'   Default: `TRUE`.
##' @param col Colours used in the plots. Defaults to `.admove_cols(10)`.
##' @param cor_tax Optional scaling factor for taxis arrows. If `NULL`,
##'   a default scaling is used internally.
##' @param cor_dif Optional scaling factor for diffusion symbols. If `NULL`,
##'   a default scaling is used internally.
##' @param cor_adv Optional scaling factor for advection arrows. If `NULL`,
##'   a default scaling is used internally.
##' @param asp Positive numeric value specifying the target aspect ratio
##'   (columns / rows) for the plot layout. Default: `2`.
##' @param plot.legend Integer controlling legend placement. If `1` (default),
##'   a shared legend is drawn in a separate panel below the plots. If `2`,
##'   the legend is added to the final plot panel.
##' @param bg Optional background colour for the plotting device. If `NULL`
##'   (default), the current background setting is used.
##' @param ... Additional arguments
##'
##' @details
##' The function is a wrapper around [plot_compare_one()] and is designed for
##' quick visual inspection of fitted models. If `auto_layout = TRUE`, the
##' plotting layout is determined automatically using [n2mfrow()], and graphical
##' parameters are reset after plotting.
##'
##' @return
##' Invisibly returns `NULL`. Called for its side effect of producing plots.
##'
##' @seealso [plot_compare_one()]
##'
##' @name plot_fit
##' @export
plot_fit <- function(x,
                     quantity = c("pref","taxis","advection",
                                  "pref_dif","dif",
                                  "par"),
                     plot_land = FALSE,
                     auto_layout = TRUE,
                     col = .admove_cols(10),
                     cor_tax = NULL,
                     cor_dif = NULL,
                     cor_adv = NULL,
                     asp = 2,
                     plot.legend = 1,
                     bg = NULL,
                     ...) {

  fit <- x

  fitlist <- list(fit)

  quantity <- match.arg(quantity, several.ok = TRUE)
  nq <- length(quantity)

  if (!is.null(bg)) {
    par(bg = bg)
  }

  ncov <- if (!is.null(fit$dat$cov)) length(fit$dat$cov) else 1L
  nsea_fit <- if (!is.null(fit$par$alpha)) dim(fit$par$alpha)[3L] else 1L
  nsea_adv <- if (!is.null(fit$par$gamma)) dim(fit$par$gamma)[3L] else 1L
  ## Covariates whose spline coefficients are all fixed (mapped NA) contribute
  ## no fitted preference relationship, so their preference panels are dropped.
  ## This happens e.g. for the current covariates under advection, where alpha
  ## (taxis) and beta (diffusion) are mapped off entirely. A covariate is "active"
  ## if any of its coefficients across knots/seasons is estimated (non-NA in map).
  active_cov <- function(map_par, par_arr) {
    if (is.null(par_arr)) return(seq_len(ncov))
    d <- dim(par_arr)                       ## [nknots, ncov, nsea]
    if (is.null(map_par)) return(seq_len(d[2L]))
    m <- array(as.integer(map_par), dim = d)
    which(apply(m, 2L, function(z) any(!is.na(z))))
  }
  sel_tax <- active_cov(fit$map$alpha, fit$par$alpha)
  ## diffusion varies with the covariate only when it has more than one knot; with
  ## a single knot it is an estimated constant and the "pref_dif" panel is
  ## uninformative, so pref_dif is shown only for multi-knot, non-fixed covariates
  nknots_dif <- if (!is.null(fit$par$beta)) dim(fit$par$beta)[1L] else 1L
  sel_dif <- if (nknots_dif > 1L) active_cov(fit$map$beta, fit$par$beta) else integer(0)
  ## the advection panels are per season rather than per covariate, so
  ## active_cov() does not apply; they are dropped when the model carries no
  ## advection at all, in which case every arrow and the whole field is zero
  adv_active <- .adv_active(fit)
  panels_per_q <- vapply(quantity, function(q) {
    if (q == "pref") length(sel_tax)
    else if (q == "taxis") nsea_fit
    else if (q == "advection") if (adv_active) nsea_adv else 0L
    else if (q == "pref_dif") length(sel_dif)
    else 1L
  }, integer(1L))
  total_panels <- sum(panels_per_q)

  if (total_panels == 0L) {
    warning("Nothing to plot: all requested quantities have zero panels ",
            "(e.g. only preference quantities were selected but every covariate ",
            "is fixed/mapped, or only advection was selected for a model fitted ",
            "without it). Nothing drawn.")
    return(invisible(NULL))
  }

  if (auto_layout) {
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar))
    mfrow <- n2mfrow(total_panels, asp = asp)
    par(mar = c(4.5,4,1,1)+0.1, oma = c(1,1,1,1))
    if(as.integer(plot.legend) == 1){
      layout(rbind(matrix(seq_len(max(total_panels, prod(mfrow))),
                          nrow = mfrow[1],
                          ncol = mfrow[2],
                          byrow = TRUE),
                   rep(total_panels + 1L, mfrow[2])),
             heights = c(rep(1, mfrow[1]), 0.15))
    }else{
      ## pad to the full grid so an odd panel count does not warn/recycle;
      ## unused regions (> total_panels) are simply left blank
      layout(matrix(seq_len(prod(mfrow)),
                    nrow = mfrow[1],
                    ncol = mfrow[2],
                    byrow = TRUE))
    }
  }

  panel_start <- cumsum(c(0L, panels_per_q[-nq]))
  all_labs <- if (nq > 1L) LETTERS[seq_len(total_panels)] else NULL
  for(i in 1:nq){
    if (panels_per_q[i] == 0L) next   ## all covariates fixed (e.g. pref_dif with a single knot)
    q_labs <- if (!is.null(all_labs)) all_labs[panel_start[i] + seq_len(panels_per_q[i])] else NULL
    q_select <- if (quantity[i] == "pref") sel_tax
                else if (quantity[i] == "pref_dif") sel_dif
                else NULL
    plot_compare_one(fitlist,
                     quantity = quantity[i],
                     col = col,
                     plot.legend = as.integer(plot.legend) == 2 && i == nq,
                     plot_land = plot_land,
                     auto_layout = TRUE,
                     panel_lab = q_labs,
                     select = q_select,
                     cor_tax = cor_tax,
                     cor_dif = cor_dif,
                     cor_adv = cor_adv,
                     bg = bg)
  }
}




## Internal functions -----------------------------------------------------------------

## TRUE when a fitted object carries a non-zero advection field: it was fitted
## with advection on and at least one gamma coefficient is non-zero. Fixed
## (mapped NA) but non-zero coefficients still produce a real field, so the test
## is on the values rather than on the map. With advection off, or every gamma
## at zero, the field is identically zero everywhere and there is nothing to
## draw -- plot_fit() drops the panels in that case.
.adv_active <- function(fit) {

  if (!isTRUE(fit$conf$use_advection)) return(FALSE)

  gamma_est <- if (!is.null(fit$pl$gamma)) fit$pl$gamma else fit$par$gamma
  if (is.null(gamma_est)) return(FALSE)

  any(gamma_est != 0, na.rm = TRUE)
}


## kappa enters the likelihood only as kappa * grad(h), and h is linear in
## alpha, so rescaling (kappa, alpha) -> (c * kappa, alpha / c) leaves the
## objective unchanged: the two are exactly confounded and the Hessian is
## singular in that direction. The default map fixes logKappa; warn if a
## user-supplied map frees it again.
##
## The one exception is a map that holds some alpha coefficient at a non-zero
## value. That pins the scale of alpha, kappa becomes identifiable, and the
## warning would be wrong -- so it is skipped in that case.
.check_kappa_map <- function(par, map, conf) {

  if (!isTRUE(conf$use_taxis)) return(invisible(NULL))
  if (is.null(par$logKappa)) return(invisible(NULL))

  kappa_free <- is.null(map$logKappa) || !all(is.na(map$logKappa))
  if (!kappa_free) return(invisible(NULL))

  alpha_fixed <- if (is.null(map$alpha)) {
    rep(FALSE, length(par$alpha))
  } else {
    is.na(map$alpha)
  }
  anchored <- any(alpha_fixed & as.vector(par$alpha) != 0)
  if (anchored) return(invisible(NULL))

  warning("logKappa is not fixed in 'map', but kappa is confounded with alpha: ",
          "only the product kappa * alpha is identifiable, so the likelihood is ",
          "flat along that direction and the Hessian is singular. Estimates and ",
          "standard errors cannot be trusted. Fix it with ",
          "map$logKappa <- factor(NA) (the default_map() behaviour) and set the ",
          "scale via par$logKappa instead. See ?default_par.",
          call. = FALSE)

  invisible(NULL)
}

## Which elements of the parameter list are estimated, and how to label them.
##
## Elements of a parameter that share a map level are a single estimated
## parameter -- they contribute one entry to opt$par and always carry identical
## estimates and standard errors -- so only the first element of each level is
## kept. The label lists all element indices the level covers (e.g. "gamma3,6"),
## following the naming of unlist(pl): no index for length-one parameters.
##
## Returns a list with 'keep' (logical, over unlist(pl)) and 'labels'
## (character, same length).
.select_estimated_par <- function(pl, map) {

  keep <- vector("list", length(pl))
  labels <- vector("list", length(pl))
  names(keep) <- names(labels) <- names(pl)

  for (nm in names(pl)) {

    n <- length(unlist(pl[[nm]]))
    mp <- map[[nm]]

    lab <- function(ii) paste0(nm, if (n > 1) paste(ii, collapse = ",") else "")

    ## no map entry: every element is estimated on its own
    if (is.null(mp)) {
      keep[[nm]] <- rep(TRUE, n)
      labels[[nm]] <- vapply(seq_len(n), lab, character(1))
      next
    }

    ## empty or malformed map entry: nothing to report
    if (length(mp) != n) {
      keep[[nm]] <- rep(FALSE, n)
      labels[[nm]] <- rep(NA_character_, n)
      next
    }

    mp <- as.character(mp)
    k <- rep(FALSE, n)
    l <- rep(NA_character_, n)

    for (lev in unique(mp[!is.na(mp)])) {
      ii <- which(mp == lev)
      k[ii[1L]] <- TRUE
      l[ii[1L]] <- lab(ii)
    }

    keep[[nm]] <- k
    labels[[nm]] <- l
  }

  list(keep = unlist(keep, use.names = FALSE),
       labels = unlist(labels, use.names = FALSE))
}


## Display labels for the estimated parameters of a fitted or simulated object,
## looked up by the element names of unlist(pl) (e.g. "gamma3" -> "gamma3,6").
## Keeps summary() and plot_compare(quantity = "par") referring to coupled
## parameters in the same way. Keys without a label are returned unchanged.
.par_display_labels <- function(x, keys) {

  pl <- if (inherits(x, "admove_sim")) x$par_sim else x$pl
  if (is.null(pl) || is.null(x$map)) return(keys)

  sel <- .select_estimated_par(pl, x$map)

  lookup <- sel$labels[sel$keep]
  names(lookup) <- names(unlist(pl))[sel$keep]

  out <- unname(lookup[keys])
  out[is.na(out)] <- keys[is.na(out)]
  out
}


.get_pl_from_opt <- function(par, map, opt) {

  pl <- par

  ind <- which(names(opt$par) == "alpha")
  pl$alpha[which(!is.na(map$alpha))] <- opt$par[ind]

  ind <- which(names(opt$par) == "beta")
  pl$beta[which(!is.na(map$beta))] <- opt$par[ind]

  ind <- which(names(opt$par) == "gamma")
  pl$gamma[which(!is.na(map$gamma))] <- opt$par[ind]

  ind <- which(names(opt$par) == "logSdO")
  pl$logSdO[which(!is.na(map$logSdO))] <- opt$par[ind]


  pl
}

## Exact AD Hessian of the fixed effects, for RTMB::sdreport(hessian.fixed=).
## Returns NULL whenever it cannot be produced, so the caller falls back to
## sdreport()'s own finite-difference default.
##
## The evaluation point is derived exactly as TMB::sdreport() derives par.fixed,
## so the two paths are comparable: obj$env$last.par.best, minus the random
## effects if there are any. admove never uses random effects, but the check
## costs nothing and keeps this correct if that ever changes -- obj$he() is the
## Hessian of the inner objective, not of the Laplace approximation, so it must
## not be used in that case.
.get_ad_hessian <- function(obj) {

  if (is.null(obj) || !is.function(obj$he)) return(NULL)

  rand <- obj$env$random
  if (!is.null(rand) && length(rand) > 0) return(NULL)

  par_fixed <- obj$env$last.par.best
  if (is.null(par_fixed)) return(NULL)

  hess <- try(obj$he(par_fixed), silent = TRUE)
  if (inherits(hess, "try-error")) return(NULL)
  if (!is.matrix(hess)) return(NULL)
  if (nrow(hess) != length(par_fixed) || ncol(hess) != length(par_fixed)) return(NULL)

  hess
}


.get_lower_bounds <- function(par){
  lower <- lapply(par, function(z) {
    z[] <- -Inf
    z
  })
  lower
}

.get_upper_bounds <- function(par){
  upper <- lapply(par, function(z) {
    z[] <- Inf
    z
  })
  upper
}



## s3 methods -------------------------------------------------------------------------

##' @rdname plot_fit
##' @export
plot.admove <- function(x, ...) {
  plot_fit(x, ...)
  return(invisible(NULL))
}

##' @rdname print-admove
##' @method print admove
##' @export
print.admove <- function(x, ...) {
  tmp <- x
  attributes(tmp) <- NULL
  NextMethod("print", tmp, ...)
}


##' @method summary admove
##' @rdname summarise_fit
##' @export
summary.admove <- function(object, ...) {
  summarise_fit(object, ...)
}


##' Log-likelihood of a fitted admove model
##'
##' @description
##' Extracts the log-likelihood from a fitted \code{admove} model. The result
##' can be passed to [stats::AIC()] or [stats::BIC()], or used for likelihood
##' ratio tests.
##'
##' @param object A fitted model object of class \code{"admove"}.
##' @param ... Currently unused.
##'
##' @return An object of class \code{"logLik"} with attributes \code{df}
##'   (number of estimated parameters) and \code{nobs} (number of used
##'   observations).
##'
##' @examples
##' fit <- admove(skjepo$sim, do_sdreport = FALSE)
##' logLik(fit)
##' AIC(fit)
##' BIC(fit)
##'
##' @importFrom stats logLik
##' @method logLik admove
##' @export
logLik.admove <- function(object, ...) {
  ll <- -object$opt$objective
  k <- length(object$opt$par)
  tags <- object$dat$tags
  if (is.data.frame(tags)) {
    nobs <- sum(tags$use == 1L, na.rm = TRUE)
  } else {
    nobs <- sum(vapply(tags, function(t) sum(t$use == 1L, na.rm = TRUE), integer(1L)))
  }
  structure(ll, df = k, nobs = nobs, class = "logLik")
}


##' Description of model
##'
##' @description Writes a string to install the version of the package which was
##'     used to run the model.
##'
##' @param fit A fitted admove list of class `admove` as returned by the
##'     function [admove].
##' @param ... Additional parameters to be passed to [message].
##'
##' @export
model_version_info <-function(fit, ...){
    .check_class(fit, "admove")

    ret <- c(
        '# The fit was run with a specific version of admove package.',
        '# If in the mean time version on your system has been updated',
        '# you can revert back to the version used by inserting this:',
        '',
        paste0('devtools::install_github("tokami/admove@',
               attr(fit, "RemoteSha"),'")'),
        '',
        '# right before the admove package is loaded'
    )

    message(ret, ...)
}
