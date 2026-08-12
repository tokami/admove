
.dxfield <- function(field, xr) {
  nr <- nrow(field)
  nc <- ncol(field)
  ret <- matrix(NA_real_, nr, nc)

  ## xr = ranges correspond to cell midpoints (thus nr-1)
  onedx <- (xr[2] - xr[1]) / (nr - 1)

  for (i in 1:nr) {
    for (j in 1:nc) {
      val <- field[i, j]
      if (is.na(val)) next

      left <- if (i > 1) field[i - 1, j] else NA_real_
      right <- if (i < nr) field[i + 1, j] else NA_real_

      if (!is.na(left) && !is.na(right)) {
        ret[i, j] <- (right - left) / (2 * onedx)
      } else if (!is.na(right)) {
        ret[i, j] <- (right - val) / onedx
      } else if (!is.na(left)) {
        ret[i, j] <- (val - left) / onedx
      }
    }
  }

  ret[is.na(ret)] <- 0

  return(ret)
}



.dyfield <- function(field, yr){
  t(.dxfield(t(field), yr))
}


.get_liv <- function(fields, r = 1) {

  ncov <- length(fields)
  nts <- sapply(fields, function(x) dim(x)[3])

  xyrange <- .get_cov_xyrange(fields)
  xr <- xyrange$xr
  yr <- xyrange$yr

  liv <- lapply(seq_len(ncov), function(i) {
    lapply(seq_len(nts[i]), function(j) {
      RTMB::interpol2Dfun(
        fields[[i]][, , j],
        xlim = xr[i, ],
        ylim = yr[i, ],
        R = r
      )
    })
  })

  liv_dx <- lapply(seq_len(ncov), function(i) {
    lapply(seq_len(nts[i]), function(j) {
      RTMB::interpol2Dfun(
        .dxfield(fields[[i]][, , j], xr[i, ]),
        xlim = xr[i, ],
        ylim = yr[i, ],
        R = r
      )
    })
  })

  liv_dy <- lapply(seq_len(ncov), function(i) {
    lapply(seq_len(nts[i]), function(j) {
      RTMB::interpol2Dfun(
        .dyfield(fields[[i]][, , j], yr[i, ]),
        xlim = xr[i, ],
        ylim = yr[i, ],
        R = r
      )
    })
  })

  list(
    liv = liv,
    liv_dx = liv_dx,
    liv_dy = liv_dy
  )
}


.make_habi <- function(liv_obj, xr, yr, time_cov,
                       s, ds, time_spline,
                       seasonal_period = NULL,
                       seasonal_cov = NULL,
                       seasonal_spline = NULL) {
  "c" <- ADoverload("c")
  "[<-" <- ADoverload("[<-")

  liv <- liv_obj$liv
  liv_dx <- liv_obj$liv_dx
  liv_dy <- liv_obj$liv_dy

  ncov <- length(liv)
  if (is.null(seasonal_period)) seasonal_period <- 1
  if (is.null(seasonal_cov)) seasonal_cov <- rep(FALSE, ncov)
  if (is.null(seasonal_spline)) seasonal_spline <- rep(FALSE, ncov)

  val <- function(xy, t){
    h <- rep(0, nrow(xy))
    for(i in 1:ncov){
      it <- t2index(t, time_cov[[i]], period = seasonal_period,
                    seasonal = seasonal_cov[i])
      is <- t2index(t, time_spline[[i]], period = seasonal_period,
                    seasonal = seasonal_spline[i])
      if(it > 0 && is > 0 && !.is_empty(s[[i]]) && !.is_empty(s[[i]][[is]])){
        h <- h + s[[i]][[is]](liv[[i]][[it]](xy[,1], xy[,2]))
      }
    }
    return(h)
  }

  grad <- function(xy, t){
    dh <- dxytmp <- RTMB::matrix(0, nrow(xy), 2)
    for(i in 1:ncov){
      it <- t2index(t, time_cov[[i]], period = seasonal_period,
                    seasonal = seasonal_cov[i])
      is <- t2index(t, time_spline[[i]], period = seasonal_period,
                    seasonal = seasonal_spline[i])
      if(it > 0 && is > 0 && !.is_empty(ds[[i]]) && !.is_empty(ds[[i]][[is]])){
        dxytmp[,1] <- liv_dx[[i]][[it]](xy[,1], xy[,2])
        dxytmp[,2] <- liv_dy[[i]][[it]](xy[,1], xy[,2])
        dh <- dh + ds[[i]][[is]](liv[[i]][[it]](xy[,1], xy[,2])) * dxytmp
      }
    }
    return(dh)
  }

  valF <- function(xy, t){
    h <- 0
    for(i in 1:ncov){
      it <- t2index(t, time_cov[[i]], period = seasonal_period,
                    seasonal = seasonal_cov[i])
      h <- h + liv[[i]][[it]](xy[,1], xy[,2])
    }
    return(h)
  }

  cov2val <- function(cov, combine = FALSE){
    nsea <- sapply(time_spline, length)
    h <- array(0, dim = c(nrow(cov), ncol(cov), max(nsea)))
    for (i in 1:ncol(cov)) {
      for (j in 1:nsea[i]) {
        is <- j
        if(!.is_empty(s[[i]]) && !.is_empty(s[[i]][[is]])){
          h[,i,j] <- h[,i,j] + s[[i]][[is]](cov[,i])
        }
      }
      if(combine){
        h <- apply(h, c(1,3), sum)
      }
    }
    return(h)
  }

  res <- list(
    val = val,
    grad = grad,
    valF = valF,
    cov2val = cov2val
  )
  return(res)
}


## Build the preference functions and the four habi objects (taxis, diffusion,
## advection in x and y) from data, configuration and parameter estimates.
## Shared by add_predictions() and .get_habi() so the two always construct them
## the same way.
.build_habi <- function(dat, conf, par_est, per) {

  pref_funcs <- .make_pref_funcs(par_est$alpha, par_est$beta, par_est$gamma,
                                 dat$knots_tax, dat$knots_dif)

  liv <- .get_liv(dat$cov)

  mk <- function(s, ds) .make_habi(liv, dat$xrange_cov,
                                   dat$yrange_cov, dat$time_cov,
                                   s, ds,
                                   dat$time_spline, per,
                                   conf$seasonal_cov,
                                   conf$seasonal_spline)

  list(pref_funcs = pref_funcs,
       habi = list(tax = mk(pref_funcs$tax, pref_funcs$dtax),
                   dif = mk(pref_funcs$dif, pref_funcs$ddif),
                   adv_x = mk(pref_funcs$adv_x, pref_funcs$dadv_x),
                   adv_y = mk(pref_funcs$adv_y, pref_funcs$dadv_y)))
}


## The habi objects are closures over RTMB::interpol2Dfun() interpolants, and
## those hold an external pointer to a C++ object. R does not serialise external
## pointers: after save()/load() (or saveRDS()/readRDS()) they come back as nil
## and any call fails with "external pointer is not valid". Probe the stored
## objects with a single evaluation and, if they are stale, rebuild them from
## the data and estimates the fit carries. Everything needed is stored, so this
## is always possible -- it just costs the interpolant setup again.
.get_habi <- function(x) {

  habi <- x$pred$habi

  if (.habi_usable(habi, x)) return(habi)

  if (is.null(x$dat$cov))
    stop("The fitted object has no covariates, so habitat fields are not ",
         "available.", call. = FALSE)

  par_est <- get_par_est(x$par, x$map, x$opt)

  .build_habi(x$dat, x$conf, par_est, period(x))$habi
}


## TRUE when the stored habi objects can still be evaluated. The probe uses the
## first prediction cell and time so that at least one covariate interpolant is
## actually hit (a time outside every covariate's range would short-circuit the
## loop in val() and hide a dead pointer).
.habi_usable <- function(habi, x) {

  if (is.null(habi) || !is.function(habi$tax$val)) return(FALSE)

  xy <- x$dat$pred$grid$xygrid
  tt <- x$dat$pred$time
  if (is.null(xy) || nrow(xy) < 1L || length(tt) < 1L) return(FALSE)

  isTRUE(tryCatch({
    habi$tax$val(xy[1L, , drop = FALSE], tt[1L])
    TRUE
  }, error = function(e) FALSE))
}
