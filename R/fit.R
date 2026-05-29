
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
  if(is.null(par)) par <- default_par(dat, conf)
  if(is.null(map)) map <- default_map(dat, conf, par)

  conf$engine <- .get_engine_integer(conf$engine)

  ## check and clean tags
  dat$tags <- check_tags(dat$tags, dat$grid, dat, conf, TRUE, verbose)

  ## period for seasonality
  dat$period <- period(dat)


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
##' @return
##' An updated object of class \code{"admove"} with sdreport results added.
##'
##' @export
add_sdreport <- function(fit, save_covariance = FALSE) {

  .check_class(fit, "admove")

  res <- fit

  t1 <- Sys.time()
  sdrep <- RTMB::sdreport(obj = fit$obj)
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
##'
##' @details
##' This function evaluates the model prediction step and stores the resulting
##' predicted quantities in the fitted object.
##'
##' @return
##' An updated object of class \code{"admove"} with model predictions added.
##'
##' @export
add_predictions <- function(fit) {

  .check_class(fit, "admove")

  res <- fit
  dat <- fit$dat
  conf <- fit$conf


  ## dimensions
  ncov <- length(dat$cov)
  nt <- length(dat$time)
  nc <- nrow(dat$grid$igrid)
  ntp <- length(dat$pred$time)
  ncp <- nrow(dat$pred$grid$igrid)


  t1 <- Sys.time()
  par_est <- get_par_est(fit$par, fit$map, fit$opt)
  kappa <- exp(par_est$logKappa)


  ## Make preference functions --------------------------
  pref_funcs <- .make_pref_funcs(par_est$alpha, par_est$beta, par_est$gamma,
                                 dat$knots_tax, dat$knots_dif)


  ## Local interpolation --------------------------------
  liv <- .get_liv(dat$cov)


  ## Setup habi objects ---------------------------------
  habi_tax <- .make_habi(liv, dat$xrange_cov,
                        dat$yrange_cov, dat$time_cov,
                        pref_funcs$tax, pref_funcs$dtax,
                        dat$time_spline, period(fit),
                        conf$seasonal_cov,
                        conf$seasonal_spline)
  habi_dif <- .make_habi(liv, dat$xrange_cov,
                        dat$yrange_cov, dat$time_cov,
                        pref_funcs$dif, pref_funcs$ddif,
                        dat$time_spline, period(fit),
                        conf$seasonal_cov,
                        conf$seasonal_spline)
  habi_adv_x <- .make_habi(liv, dat$xrange_cov,
                          dat$yrange_cov, dat$time_cov,
                          pref_funcs$adv_x, pref_funcs$dadv_x,
                          dat$time_spline, period(fit),
                          conf$seasonal_cov,
                          conf$seasonal_spline)
  habi_adv_y <- .make_habi(liv, dat$xrange_cov,
                          dat$yrange_cov, dat$time_cov,
                          pref_funcs$adv_y, pref_funcs$dadv_y,
                          dat$time_spline, period(fit),
                          conf$seasonal_cov,
                          conf$seasonal_spline)


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


##' Compute predicted location distributions for a single tag
##'
##' @description
##' Propagates the model's predicted spatial location distribution for one
##' archival tag from its release to recovery time, using either the CTMC
##' forward-pass (engine 2) or repeated Kalman-filter track simulations
##' (engine 1). The result is stored in `fit$tag_dist` and consumed by
##' [plot_tag_dist()]. Separating the expensive computation from rendering
##' means plot aesthetics can be changed without re-running the model.
##'
##' @param fit A fitted object of class `admove`, as returned by [admove()].
##' @param i Integer index of the tag to use. Default is `1`.
##' @param dt Time step used when simulating tracks (engine 1 only). Default
##'   is `0.5`.
##' @param engine Optional integer overriding the engine stored in
##'   `fit$conf$engine`. `1` = Kalman filter, `2` = CTMC.
##' @param xrel0,yrel0 Optional coordinates overriding the release location
##'   for CTMC-based predictions.
##'
##' @return
##' A copy of `fit` with the additional component `$tag_dist`, a list
##' containing the precomputed distributions and metadata required by
##' [plot_tag_dist()].
##'
##' @seealso [plot_tag_dist()]
##'
##' @export
add_tag_dist <- function(fit, i = 1, dt = 0.5,
                         engine = NULL, xrel0 = NULL, yrel0 = NULL) {

  .check_class(fit, "admove")

  dat <- fit$dat
  conf <- fit$conf

  if (inherits(dat$tags, "data.frame")) {
    tags <- split(dat$tags, dat$tags$id)
  } else {
    tags <- dat$tags
  }

  if (i < 1L || i > length(tags))
    stop("Tag index 'i' (", i, ") is out of range [1, ", length(tags), "].")

  tag <- tags[[i]]
  ind <- which(apply(!is.na(tag[, 1:3]), 1, all))

  if (length(ind) < 2L)
    stop("Tag ", i, " has fewer than 2 non-missing (t, x, y) observations; ",
         "only archival tags with full position records are supported.")

  tag <- tag[ind, ]
  xrel <- if (!is.null(xrel0)) xrel0 else tag[1, 2]
  yrel <- if (!is.null(yrel0)) yrel0 else tag[1, 3]
  trel <- tag[1, 1]
  trec <- tag[nrow(tag), 1]

  if (is.null(engine)) engine <- conf$engine

  if (engine == 2L) {

    tall <- dat$pred$time
    mstar <- calc_mstar(fit)

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
        mstar[,, itrel + k - 1L] * diff(dat$pred$time)[k]))
      dist_prob[k + 1L, ] <- as.vector(dist_prob[k, ] %*% m)
    }

    dens_list <- vector("list", nrow(tag))
    for (k in seq_along(tind)) {
      ct <- dat$pred$grid$celltable
      ct[which(!is.na(ct))] <- dist_prob[tind[k], ]
      dens_list[[k]] <- ct
    }

    tag_dist <- list(
      engine = engine,
      tag = tag,
      i = i,
      dens_list = dens_list,
      xg = x_centers(dat$pred$grid),
      yg = y_centers(dat$pred$grid),
      xrange = dat$grid$xrange,
      yrange = dat$grid$yrange
    )

  } else {

    funcs <- default_sim_funcs(dat, conf, fit$pl)
    dt_min <- min(dat$min_dt, median(diff(tag$t)))

    out <- build_time(tag$t, mode = "fixed_dt",
                           dt_min = dt_min, dt = dt, eps = 1)
    ts <- out$ts
    dts <- out$dts
    nts <- out$nts
    observed <- out$observed

    kappa <- exp(fit$pl$logKappa)

    ## Single deterministic forward pass: propagate mean (xy0) and variance
    ## (PP) only. No random draws are needed — the stochastic xy track was
    ## never used for anything other than evaluating taxis/diffusion at a
    ## noisy position, which is well-approximated by evaluating at the mean.
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

    tag_dist <- list(
      engine = engine,
      tag = tag,
      i = i,
      traj = traj,
      ind.track = ind.track,
      xrange = xrange,
      yrange = yrange
    )
  }

  fit$tag_dist <- tag_dist
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

  idx <- lapply(x$map, function(x) if(length(x) == 0) FALSE else !is.na(x))
  tmp <- x$pl[!(names(x$pl) %in% names(x$map))]
  idx <- c(idx, lapply(tmp, function(x) rep(TRUE, length(x))))
  pl <- unlist(x$pl)[unlist(idx[match(names(x$pl), names(idx))])]
  plsd <- unlist(x$plsd)[unlist(idx[match(names(x$plsd), names(idx))])]
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
##'     \item{`"pref"`}{Habitat preference functions.}
##'     \item{`"taxis"`}{Taxis (movement direction and magnitude).}
##'     \item{`"dif"`}{Diffusion.}
##'     \item{`"par"`}{Estimated model parameters.}
##'   }
##'   Multiple quantities can be selected.
##' @param plot_land Logical; if `TRUE`, land masses are added to spatial plots
##'   using [plot_land()]. Default: `FALSE`.
##' @param auto_layout Logical; if `TRUE`, graphical parameters are set and
##'   restored automatically, and plots are arranged in a multi-panel layout.
##'   Default: `TRUE`.
##' @param col Colours used in the plots. Defaults to `.admove_cols(10)`.
##' @param cor_dif Optional scaling factor for diffusion symbols. If `NULL`,
##'   a default scaling is used internally.
##' @param cor_tax Optional scaling factor for taxis arrows. If `NULL`,
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
                     quantity = c("pref","taxis",
                                  "dif","par"),
                     plot_land = FALSE,
                     auto_layout = TRUE,
                     col = .admove_cols(10),
                     cor_dif = NULL,
                     cor_tax = NULL,
                     asp = 2,
                     plot.legend = 1,
                     bg = NULL,
                     ...) {

  fit <- x

  fitlist <- list(fit)

  quantity <- match.arg(quantity, several.ok = TRUE)
  nq <- length(quantity)

  if(!is.null(bg)){
    par(bg = bg)
  }

  ncov <- if (!is.null(fit$dat$cov)) length(fit$dat$cov) else 1L
  panels_per_q <- vapply(quantity, function(q) {
    if (q == "pref") ncov else 1L
  }, integer(1L))
  total_panels <- sum(panels_per_q)

  if(auto_layout){
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
      layout(matrix(seq_len(total_panels),
                    nrow = mfrow[1],
                    ncol = mfrow[2],
                    byrow = TRUE))
    }
  }

  panel_start <- cumsum(c(0L, panels_per_q[-nq]))
  all_labs <- if (nq > 1L) LETTERS[seq_len(total_panels)] else NULL
  for(i in 1:nq){
    q_labs <- if (!is.null(all_labs)) all_labs[panel_start[i] + seq_len(panels_per_q[i])] else NULL
    plot_compare_one(fitlist,
                     quantity = quantity[i],
                     col = col,
                     plot.legend = as.integer(plot.legend) == 2 && i == nq,
                     plot_land = plot_land,
                     auto_layout = TRUE,
                     panel_lab = q_labs,
                     cor_dif = cor_dif,
                     cor_tax = cor_tax,
                     bg = bg)
  }
}




## Internal functions -----------------------------------------------------------------

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
