
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
      if(it > 0 && !.is_empty(s[[i]]) && !.is_empty(s[[i]][[is]])){
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
      if(it > 0 && !.is_empty(ds[[i]]) && !.is_empty(ds[[i]][[is]])){
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
