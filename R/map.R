
## Main functions ---------------------------------------------------------------------

##' Default parameter map
##'
##' @description
##' Construct the default `map` list used to control which parameters are
##' estimated, fixed, or coupled in TMB. The returned list has the same overall
##' structure as the model parameter list and is intended to be passed to
##' `TMB::MakeADFun()` via its `map` argument.
##'
##' By default:
##' \itemize{
##'   \item taxis spline coefficients (`alpha`) are mapped using
##'     `.make_alpha_map()`, with the first row fixed and remaining coefficients
##'     estimated independently;
##'   \item the taxis scaling parameter (`logKappa`) is fixed;
##'   \item diffusion spline coefficients (`beta`) are mapped using
##'     `.make_beta_map()`;
##'   \item advection coefficients (`gamma`) are either fixed or, if advection is
##'     enabled, coupled between the \(x\)- and \(y\)-directions within each
##'     covariate;
##'   \item observation-error parameters (`logSdO`) are fixed unless estimation
##'     is enabled for the corresponding tag type, in which case \(x\)- and
##'     \(y\)-direction standard deviations are coupled by default.
##' }
##'
##' @param dat A data list as produced by [setup_data()].
##' @param conf A configuration list as produced by [default_conf()].
##' @param par A parameter list with initial values as produced by
##'   [default_par()].
##'
##' @return
##' A named list of factors with elements `alpha`, `logKappa`, `beta`, `gamma`,
##' and `logSdO`. Entries with `NA` are fixed, while equal factor levels are
##' estimated as the same parameter.
##'
##' @examples
##' map <- with(skjepo$sim, default_map(dat, conf, par))
##'
##' @export
default_map <- function(dat, conf, par){

  map <- list()

  ## Taxis ----------------------------------------------------
  map$alpha <- .make_alpha_map(par$alpha)

  ## Taxis scaling --------------------------------------------
  map$logKappa <- factor(NA)

  ## Diffusion ------------------------------------------------
  map$beta <- .make_beta_map(par$beta)

  ## Advection ------------------------------------------------
  if(conf$use_advection){
    ## Link x and y direction by default
    map$gamma <- factor(sapply(1:ncol(par$gamma),
                               function(x) rep(x, nrow(par$gamma))))
  }else{
    map$gamma <- factor(rep(NA, length(par$gamma)))
  }

  ## Observation error ----------------------------------------
  logSdO <- rep(NA, 6)

  if(conf$use_dtags && conf$obs_var_type[1] %in% c(1,2)){
    logSdO[1:2] <- c(1, 1) ## assume equal var in x,y
  }

  if(conf$use_stags && conf$obs_var_type[2] %in% c(1,2)){
    mini <- min(logSdO[1:2])
    if (is.na(mini)) mini <- 0
    logSdO[3:4] <- c(mini + 1, mini + 1) ## assume equal var in x,y
  }

  map$logSdO <- factor(logSdO)

  ## return
  map
}


## Internal functions ---------------------------------------------------------------

.make_alpha_map <- function(alpha) {

  dims <- dim(alpha)
  stopifnot(length(dims) == 3)

  nknot <- dims[1]
  ncov <- dims[2]
  nsea <- dims[3]

  map_id <- array(NA_integer_, dim = dims, dimnames = dimnames(alpha))
  next_id <- 1L

  if (nknot > 1) {
    for (s in seq_len(nsea)) {
      for (j in seq_len(ncov)) {
        for (k in 2:nknot) {
          map_id[k, j, s] <- next_id
          next_id <- next_id + 1L
        }
      }
    }
  }

  factor(map_id)
}

.make_beta_map <- function(beta) {

  dims <- dim(beta)
  stopifnot(length(dims) == 3)

  nknot <- dims[1]
  ncov <- dims[2]
  nsea <- dims[3]

  map_id <- array(NA_integer_, dim = dims, dimnames = dimnames(beta))
  next_id <- 1L

  ## intercept: first covariate only, coupled across seasons
  map_id[1, 1, ] <- next_id
  next_id <- next_id + 1L

  ## remaining spline coefficients: independent
  if (nknot > 1) {
    for (s in seq_len(nsea)) {
      for (j in seq_len(ncov)) {
        for (k in 2:nknot) {
          map_id[k, j, s] <- next_id
          next_id <- next_id + 1L
        }
      }
    }
  }

  factor(map_id)
}
