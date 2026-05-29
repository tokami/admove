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


  ## testing
  liv <- .get_liv(dat$cov)


  ## Make preference functions --------------------------
  pref_funcs <- .make_pref_funcs(par$alpha, par$beta, par$gamma,
                                dat$knots_tax, dat$knots_dif)

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

      dt_min <- min(dat$min_dt, median(diff(tag$t), na.rm = TRUE), na.rm = TRUE)
      if (is.na(dt_min)) browser()
      if (time_mode == "fill_gaps" && (dt_min == 0 || is.na(dt_min))) next()

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
            if ((dat$obs_var_type[ind_tt] == 1 && ind_obs_j != nts) ||
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

            ## likelihood
            loglik_tags[i] <- loglik_tags[i] + RTMB::dnorm(w[1], 0, sqrt(F[1]), TRUE)
            loglik_tags[i] <- loglik_tags[i] + RTMB::dnorm(w[2], 0, sqrt(F[2]), TRUE)

            ## update
            if (isTRUE(dat$do_update[ind_tt])) {
              last_xy <- pred_xy + PP / F * w
              P <- PP - PP / F * PP
            } else {
              last_xy <- pred_xy
              P <- PP
            }

          }
        } else {
          last_xy <- pred_xy
          P <- PP
        }
      }
    }

  } else if (dat$engine == 2) {  ## CTMC

    flag_expm_uni <- ifelse(dat$ctmc_method == 2, TRUE, FALSE)

    time_mode <- ifelse(is.null(dat$dt) || is.na(dat$dt),
                        "fill_gaps", "fixed_dt")


    if (dat$ctmc_method > 0) {
      mstar_template <- make_mstar_template(nextTo, ad = TRUE)
    }

    for (i in seq_len(ntags)) {

      tag <- dat$tags[[i]]
      nobs <- nrow(tag)

      dt_min <- min(dat$min_dt, median(diff(tag$t)))
      if (time_mode == "fill_gaps" && (dt_min == 0 || is.na(dt_min))) next()

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
        if (dat$ctmc_method > 0) {
          Zstar <- Astar <- Dstar <- mstar_template
          Zstar@x[] <- Astar@x[] <- Dstar@x[] <- 0
        } else {
          Zstar <- Astar <- Dstar <- RTMB::matrix(0, nc, nc)
        }

        ## taxis
        if (dat$use_taxis) {
          move <- kappa * habi_tax$grad(xygrid, ts[t-1]) * dt  ## distance
          Zstar <- fill_inst_mat(Zstar, move, nextTo, next_dist)
        }

        ## advection
        if (dat$use_advection) {
          move <- cbind(habi_adv_x$val(xygrid, ts[t-1]),
                        habi_adv_y$val(xygrid, ts[t-1])) * dt  ## distance
          Astar <- fill_inst_mat(Astar, move, nextTo, next_dist)
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
        if (dat$ctmc_method > 0) {

          pred_dist <- as.vector(RTMB::expAv(Mstar,
                                     last_dist,
                                     transpose = TRUE,
                                     uniformization = flag_expm_uni,
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
            if ((dat$obs_var_type[ind_tt] == 1 && ind_obs_j != nts) ||
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

            ## likelihood
            update_dist <- pred_dist * this_dist
            loglik_tags[i] <- loglik_tags[i] + log(sum(update_dist))

            ## update
            if (isTRUE(dat$do_update[ind_tt])) {
              last_dist <- update_dist / sum(update_dist)
            } else {
              last_dist <- pred_dist
            }
          }
        } else {
          last_dist <- pred_dist
        }
      }
    }

  } else {

    stop("This engine is not yet implemented. Select 1 for Kalman filter and 2 for CTMC.")

  }

  nll <- nll - sum(loglik_tags)

  REPORT(loglik_tags)


  ## Predicted habi preference functions ------------------------------
  pref_taxis_pred <- habi_tax$cov2val(dat$pred$cov)
  pref_dif_pred <- habi_dif$cov2val(dat$pred$cov)


  ADREPORT(pref_taxis_pred)
  ADREPORT(pref_dif_pred)


  return(nll)
}
