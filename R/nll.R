##' Negative log-likelihood for the admove model
##'
##' @description
##' Computes the negative log-likelihood of the admove movement model given
##' parameter values and input data.
##'
##' @param par A named list of model parameters. The structure should match the
##'   output of [default_par()].
##' @param dat A list containing model input data and configuration settings,
##'   as returned by [setup_data()] and [default_conf()], respectively.
##'
##' @details
##' Two estimation engines are currently implemented:
##'
##' \itemize{
##'   \item \strong{Kalman filter (KF):} continuous space, discrete time formulation
##'   \item \strong{Continuous-time Markov chain (CTMC):} discrete space, continuous time formulation
##' }
##'
##' Both engines operate on the same input data structure. The choice of engine
##' is controlled via \code{conf$engine}, where \code{1} selects the Kalman filter
##' and \code{2} selects the CTMC approach.
##'
##' @return A numeric value representing the negative log-likelihood.
##'
##' @importFrom Matrix expm
##'
nll <- function(par, dat) {


  ## R's JIT does not preserve a few specific definitions made by RTMB
  "c" <- RTMB::ADoverload("c")
  "[<-" <- RTMB::ADoverload("[<-")
  "diag<-" <- RTMB::ADoverload("diag<-")


  nll <- 0
  loglik_tags <- rep(0, length(dat$tags))


  ## Conversions
  sdO <- exp(par$logSdO)
  varO <- sdO^2
  kappa <- exp(par$logKappa)


  ## Dimensions
  ncov <- length(dat$cov)
  nt <- length(dat$time_cont)
  nc <- nrow(dat$grid$igrid)


  ## Variables
  move0 <- matrix(0, 1, 2)
  xygrid <- dat$grid$xygrid
  cs <- dat$grid$cellsize
  nextTo <- get_neighbours(dat$grid)
  next_dist <- c(dat$grid$cellsize[1], dat$grid$cellsize[1],
                 dat$grid$cellsize[2], dat$grid$cellsize[2])
  ntags <- length(dat$tags)
  boundary_excess <- rep(0, ntags)

  ## Predicted observation distribution, one row per row of the tag data in the
  ## order the tags are split (see tag_predictions()): the mean and variance the
  ## likelihood evaluates each observation against, i.e. before its update and
  ## including observation error. pred_set marks the rows that were filled, and
  ## is plain numeric because only non-AD code writes it.
  nobs_tags <- vapply(dat$tags, nrow, integer(1))
  obs_offset <- c(0L, cumsum(nobs_tags))
  nobs_all <- sum(nobs_tags)
  pred_x <- pred_y <- pred_var_x <- pred_var_y <- rep(0, nobs_all)
  pred_set <- rep(0, nobs_all)


  ## Covariate field limits for the KF. The interpolants are only defined
  ## between the outermost cell centres and return NaN beyond, so a predicted
  ## position that leaves the field (e.g. a mark-recapture track driven by a
  ## strong taxis during the first optimizer steps) turns the whole likelihood
  ## NaN and the optimizer stalls. With kf_boundary = "clamp" predicted
  ## positions are held at the edge instead. As long as no position leaves the
  ## field this is the identity, so such likelihoods are unchanged.
  kf_clamp <- dat$engine == 1 && !is.null(dat$xrange_cov) &&
    !identical(dat$kf_boundary, "none")
  if (kf_clamp) {
    ## RTMB's interpolant is also NaN within ~1e-4 cell widths of the centre of
    ## an NA cell, and a position clamped in both directions sits exactly on the
    ## centre of a corner cell, which is often land. So hold positions a thousandth
    ## of a cell inside the outermost cell centres.
    cs_field <- min(sapply(seq_along(dat$cov), function(k) {
      c(diff(dat$xrange_cov[k, ]) / max(1, dim(dat$cov[[k]])[1] - 1),
        diff(dat$yrange_cov[k, ]) / max(1, dim(dat$cov[[k]])[2] - 1))
    }))
    eps_field <- 1e-3 * cs_field
    xlim_field <- c(max(dat$xrange_cov[, 1]), min(dat$xrange_cov[, 2])) +
      c(eps_field, -eps_field)
    ylim_field <- c(max(dat$yrange_cov[, 1]), min(dat$yrange_cov[, 2])) +
      c(eps_field, -eps_field)
    clamp_xy <- function(xy) {
      out <- xy
      out[, 1] <- .clamp_ad(xy[, 1], xlim_field[1], xlim_field[2])
      out[, 2] <- .clamp_ad(xy[, 2], ylim_field[1], ylim_field[2])
      out
    }
  }


  ## testing
  liv <- .get_liv(dat$cov)


  ## Make preference functions --------------------------
  pref_funcs <- .make_pref_funcs(par$alpha, par$beta, par$gamma,
                                dat$knots_tax, dat$knots_dif,
                                method = dat$smooth_method)

  ## Setup habi objects ---------------------------------
  habi_dif <- .make_habi(liv, dat$xrange_cov,
                        dat$yrange_cov, dat$time_cov,
                        pref_funcs$dif, pref_funcs$ddif,
                        dat$time_spline, dat$period,
                        dat$seasonal_cov,
                        dat$seasonal_spline)
  habi_tax <- .make_habi(liv, dat$xrange_cov,
                        dat$yrange_cov, dat$time_cov,
                        pref_funcs$tax, pref_funcs$dtax,
                        dat$time_spline, dat$period,
                        dat$seasonal_cov,
                        dat$seasonal_spline)
  habi_adv_x <- .make_habi(liv, dat$xrange_cov,
                          dat$yrange_cov, dat$time_cov,
                          pref_funcs$adv_x, pref_funcs$dadv_x,
                          dat$time_spline, dat$period,
                          dat$seasonal_cov,
                          dat$seasonal_spline)
  habi_adv_y <- .make_habi(liv, dat$xrange_cov,
                          dat$yrange_cov, dat$time_cov,
                          pref_funcs$adv_y, pref_funcs$dadv_y,
                          dat$time_spline, dat$period,
                          dat$seasonal_cov,
                          dat$seasonal_spline)


  ## Estimate movement ------------------------------------
  if (dat$engine == 1) { ## KF

    time_mode <- ifelse(is.null(dat$dt) || is.na(dat$dt),
                        "fill_gaps", "fixed_dt")

    for (i in seq_len(ntags)) {

      tag <- dat$tags[[i]]
      last_xy <- as.matrix(tag[1,2:3,drop = FALSE])

      ind_tag_type <- as.integer(tag$tag_type)

      ## Ambiguous final observation: several candidate positions, exactly one
      ## of which is the true one, with known probabilities. Their log-densities
      ## are collected here and combined into a single mixture term after the
      ## time loop.
      ev <- .tag_events(tag)
      last_ev <- max(ev)
      amb <- sum(ev == last_ev) > 1L
      amb_lw <- NULL

      ## unique(): candidate positions may share a time, and a repeated time
      ## would make the median step 0 and silently skip the whole tag below
      dt_min <- min(dat$min_dt, median(diff(sort(unique(tag$t))), na.rm = TRUE), na.rm = TRUE)
      if (time_mode == "fill_gaps" && (!is.finite(dt_min) || dt_min <= 0)) next()

      out <- build_time(tag$t,
                        mode = time_mode,
                        dt_min = dt_min,
                        dt = dat$dt,
                        eps = 0.1)

      ts <- out$ts
      dts <- out$dts
      nts <- out$nts
      observed <- out$observed
      if (length(observed) == 0) next()

      if (nts < 2) stop("Something went wrong (nts < 2).")

      ## initial uncertainty
      P <- dat$p_init

      for (t in 2:nts) {

        dt <- dts[t-1]

        moveT <- moveA <- move0

        ## the updated position can lie just outside the field when an
        ## observation does (outer half of a boundary cell)
        if (kf_clamp) last_xy <- clamp_xy(last_xy)

        ## diffusion
        D <- exp(habi_dif$val(last_xy, ts[t-1]))

        ## taxis
        if (dat$use_taxis) {
          moveT <- kappa * habi_tax$grad(last_xy, ts[t-1]) * dt
        }

        ## advection
        if (dat$use_advection) {
          moveA <- c(habi_adv_x$val(last_xy, ts[t-1]),
                     habi_adv_y$val(last_xy, ts[t-1])) * dt
        }

        pred_xy <- last_xy + moveT + moveA
        if (kf_clamp) {
          pred_in <- clamp_xy(pred_xy)
          boundary_excess[i] <- boundary_excess[i] +
            sum(abs(pred_xy - pred_in))
          pred_xy <- pred_in
        }
        PP <- P + (2 * D * dt)

        if (t %in% observed) {

          ind_obs <- which(observed == t) + 1

          ## multiple observation in same time
          for (j in seq_along(ind_obs)) {
            ind_obs_j <- ind_obs[j]
            ind_tt <- ind_tag_type[ind_obs_j]

            if (is.na(tag$ic[ind_obs_j]) || tag$use[ind_obs_j] == 0) {
              last_xy <- pred_xy
              P <- PP
              next()
            }

            ## default
            F <- PP

            ## obs uncertainty
            ## obs_var_type 1 means "all but the last observation", so the
            ## comparison is against the last observation EVENT of this tag --
            ## not against nts, which counts the integration time points and is
            ## larger than nrow(tag) whenever gaps are filled (the default);
            ## comparing to nts made type 1 behave like type 2. Comparing
            ## events rather than rows keeps every candidate position of an
            ## ambiguous final observation on the same footing.
            if ((dat$obs_var_type[ind_tt] == 1 && ev[ind_obs_j] != last_ev) ||
                  dat$obs_var_type[ind_tt] == 2 ||
                  dat$obs_var_type[ind_tt] == 3) {
              if (dat$obs_var_type[ind_tt] == 3) {
                pxy <- PP + c(tag$sdx[ind_obs_j], tag$sdy[ind_obs_j])^2
              } else {
                pxy <- PP + varO[1:2,ind_tt]
              }

              if (all(!is.na(pxy))) {
                F <- pxy
              }
            }

            this_xy <- c(tag$x[ind_obs_j], tag$y[ind_obs_j])
            w <- this_xy - pred_xy

            ind_pred <- obs_offset[i] + ind_obs_j
            pred_x[ind_pred] <- pred_xy[1]
            pred_y[ind_pred] <- pred_xy[2]
            pred_var_x[ind_pred] <- F[1]
            pred_var_y[ind_pred] <- F[2]
            pred_set[ind_pred] <- 1

            ld <- RTMB::dnorm(w[1], 0, sqrt(F[1]), TRUE) +
              RTMB::dnorm(w[2], 0, sqrt(F[2]), TRUE)

            if (amb && ev[ind_obs_j] == last_ev) {

              ## One of the candidates is the true position, so their densities
              ## are summed (weighted), not multiplied. The state is deliberately
              ## NOT updated: the ambiguous event is the last one, so nothing is
              ## propagated past it and every candidate is evaluated against the
              ## same release-conditioned prediction. That keeps the mixture
              ## exact -- there is no Gaussian-mixture posterior to collapse.
              ##
              ## RTMB's c() coerces every argument with advector(), which
              ## rejects NULL, so the first term has to seed the vector.
              term <- log(tag$prob[ind_obs_j]) + ld
              amb_lw <- if (is.null(amb_lw)) term else c(amb_lw, term)
              last_xy <- pred_xy
              P <- PP

            } else {

              ## likelihood
              loglik_tags[i] <- loglik_tags[i] + ld

              ## update
              if (isTRUE(dat$do_update[ind_tt])) {
                last_xy <- pred_xy + PP / F * w
                P <- PP - PP / F * PP
              } else {
                last_xy <- pred_xy
                P <- PP
              }
            }

          }
        } else {
          last_xy <- pred_xy
          P <- PP
        }
      }

      ## log sum_k prob_k * density_k, accumulated stably (see .logsumexp_ad)
      if (amb && !is.null(amb_lw)) {
        loglik_tags[i] <- loglik_tags[i] + .logsumexp_ad(amb_lw)
      }
    }

  } else if (dat$engine == 2) {  ## CTMC

    time_mode <- ifelse(is.null(dat$dt) || is.na(dat$dt),
                        "fill_gaps", "fixed_dt")


    if (identical(dat$ctmc_method, "expav")) {
      mstar_template <- make_mstar_template(nextTo, ad = TRUE)
    }

    for (i in seq_len(ntags)) {

      tag <- dat$tags[[i]]

      ## ambiguous final observation -- see the KF branch above
      ev <- .tag_events(tag)
      last_ev <- max(ev)
      amb <- sum(ev == last_ev) > 1L
      amb_lw <- NULL

      ## unique(): candidate positions may share a time
      dt_min <- min(dat$min_dt, median(diff(sort(unique(tag$t)))))
      if (time_mode == "fill_gaps" && (!is.finite(dt_min) || dt_min <= 0)) next()

      out <- build_time(tag$t,
                        mode = time_mode,
                        dt_min = dt_min,
                        dt = dat$dt,
                        eps = 0.1)

      ts <- out$ts
      dts <- out$dts
      nts <- out$nts
      observed <- out$observed
      if (length(observed) == 0) next()

      ind_tag_type <- as.integer(tag$tag_type)

      if (tag$use[1] == 0) next()

      ## Distribution probability
      last_dist <- rep(0, nc)
      last_dist[tag$ic[1]] <- 1

      if (nts < 2) stop("Something went wrong (nts < 2).")

      ## Loop over time
      for (t in 2:nts) {

        dt <- dts[t-1]

        ## Set to zero
        if (identical(dat$ctmc_method, "expav")) {
          Zstar <- Astar <- Dstar <- mstar_template
          Zstar@x[] <- Astar@x[] <- Dstar@x[] <- 0
        } else {
          Zstar <- Astar <- Dstar <- RTMB::matrix(0, nc, nc)
        }

        ## taxis
        if (dat$use_taxis) {
          move <- kappa * habi_tax$grad(xygrid, ts[t-1]) * dt  ## distance
          Zstar <- fill_inst_mat(Zstar, move, nextTo, next_dist, dat$drift_scheme)
        }

        ## advection
        if (dat$use_advection) {
          move <- cbind(habi_adv_x$val(xygrid, ts[t-1]),
                        habi_adv_y$val(xygrid, ts[t-1])) * dt  ## distance
          Astar <- fill_inst_mat(Astar, move, nextTo, next_dist, dat$drift_scheme)
        }

        ## diffusion
        D <- exp(habi_dif$val(xygrid, ts[t-1])) ## distance^2 / time
        hD <- D * dt  ## (distance^2)
        for (k in 1:4) {
          j <- k + 1
          ind <- which(!is.na(nextTo[, j]))
          Dstar[cbind(ind, nextTo[ind, j])] <- hD[ind] / next_dist[k]^2
          ## (distance^2) / (distance) = distance
        }

        ## Movement rates
        Mstar <- Zstar + Astar + Dstar

        ## Mass balance
        Mstar[cbind(1:nc, 1:nc)] <- 0
        Mstar[cbind(1:nc, 1:nc)] <- -RTMB::rowSums(Mstar)

        ## dist prob after move
        if (identical(dat$ctmc_method, "expav")) {

          pred_dist <- as.vector(RTMB::expAv(Mstar,
                                     last_dist,
                                     transpose = TRUE,
                                     uniformization = TRUE,
                                     rescale_freq = 1,
                                     trace = FALSE))

        } else {

          M <- Matrix::expm(Mstar)
          pred_dist <- as.vector(RTMB::matrix(last_dist, 1, nc) %*% M)

        }

        if (t %in% observed) {

          ind_obs <- which(observed == t) + 1

          ## multiple observation in same time
          for (j in seq_along(ind_obs)) {
            ind_obs_j <- ind_obs[j]
            ind_tt <- ind_tag_type[ind_obs_j]

            if (is.na(tag$ic[ind_obs_j]) || tag$use[ind_obs_j] == 0) {
              last_dist <- pred_dist
              next()
            }

            ## default
            this_dist <- rep(0, nrow(xygrid))
            this_dist[tag$ic[ind_obs_j]] <- 1

            ## obs uncertainty
            ## see the note in the KF branch: "all but the last observation" is
            ## a comparison against the last observation event, not against the
            ## number of integration time points
            if ((dat$obs_var_type[ind_tt] == 1 && ev[ind_obs_j] != last_ev) ||
                  dat$obs_var_type[ind_tt] == 2 ||
                  dat$obs_var_type[ind_tt] == 3) {

              xLo <- xygrid[,1] - cs[1] / 2
              xUp <- xygrid[,1] + cs[1] / 2
              yLo <- xygrid[,2] - cs[2] / 2
              yUp <- xygrid[,2] + cs[2] / 2

              xObs <- tag$x[ind_obs_j]
              yObs <- tag$y[ind_obs_j]
              if (dat$obs_var_type[ind_tt] == 3) {
                sdx <- tag$sdx[ind_obs_j]
                sdy <- tag$sdy[ind_obs_j]
              } else {
                sdx <- sdO[1,ind_tt]
                sdy <- sdO[2,ind_tt]
              }
              px <- RTMB::pnorm(xUp, mean = xObs, sd = sdx) -
                RTMB::pnorm(xLo, mean = xObs, sd = sdx)
              py <- RTMB::pnorm(yUp, mean = yObs, sd = sdy) -
                RTMB::pnorm(yLo, mean = yObs, sd = sdy)
              pxy <- px * py
              pxy <- pxy / sum(pxy)

              if (all(!is.na(pxy))) {
                this_dist <- pxy
              }
            }

            update_dist <- pred_dist * this_dist

            if (amb && ev[ind_obs_j] == last_ev) {

              ## mixture over the candidate positions; no update, exactly as in
              ## the KF branch (the ambiguous event is the last one)
              term <- log(tag$prob[ind_obs_j]) + log(sum(update_dist))
              amb_lw <- if (is.null(amb_lw)) term else c(amb_lw, term)
              last_dist <- pred_dist

            } else {

              ## likelihood
              loglik_tags[i] <- loglik_tags[i] + log(sum(update_dist))

              ## update
              if (isTRUE(dat$do_update[ind_tt])) {
                last_dist <- update_dist / sum(update_dist)
              } else {
                last_dist <- pred_dist
              }
            }
          }
        } else {
          last_dist <- pred_dist
        }
      }

      ## log sum_k prob_k * density_k, accumulated stably (see .logsumexp_ad)
      if (amb && !is.null(amb_lw)) {
        loglik_tags[i] <- loglik_tags[i] + .logsumexp_ad(amb_lw)
      }
    }

  } else {

    stop("This engine is not yet implemented. Select \"kf\" for the Kalman filter or \"ctmc\" for the CTMC.")

  }

  nll <- nll - sum(loglik_tags)

  REPORT(loglik_tags)


  ## Predicted habi preference functions ------------------------------
  pref_taxis_pred <- habi_tax$cov2val(dat$pred$cov)
  pref_dif_pred <- habi_dif$cov2val(dat$pred$cov)


  REPORT(pref_taxis_pred)
  REPORT(pref_dif_pred)
  ADREPORT(pref_taxis_pred)
  ADREPORT(pref_dif_pred)
  ## per tag: summed distance by which predicted positions were moved back to
  ## the edge of the covariate field (KF, conf$kf_boundary = "clamp")
  REPORT(boundary_excess)

  ## predicted observation distribution per tag row (KF only; see above)
  REPORT(pred_x)
  REPORT(pred_y)
  REPORT(pred_var_x)
  REPORT(pred_var_y)
  REPORT(pred_set)


  return(nll)
}
