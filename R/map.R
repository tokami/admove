
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
##'     estimated independently for each season;
##'   \item the taxis scaling parameter (`logKappa`) is fixed;
##'   \item diffusion spline coefficients (`beta`) are mapped using
##'     `.make_beta_map()` and are **coupled across seasons**, so diffusion is
##'     constant over the seasonal cycle unless `conf$seasonal_dif` is `TRUE`;
##'   \item advection coefficients (`gamma`, one per advection field, see
##'     [prep_adv()]) are fixed when `conf$use_advection` is `FALSE`. Otherwise
##'     there is one estimated coefficient per field and advection season,
##'     shared by the \(x\)- and \(y\)-components (`conf$adv_gamma =
##'     "shared"`), or one per component (`"xy"`);
##'   \item the constant drift (`adv_const`) is estimated, separately in \(x\)
##'     and \(y\) and per advection season, only when `conf$adv_const` is
##'     `TRUE`, and fixed otherwise;
##'   \item observation-error parameters (`logSdO`) are fixed unless estimation
##'     is enabled for the corresponding tag type via `conf$obs_var_type`, in
##'     which case \(x\)- and \(y\)-direction standard deviations are coupled by
##'     default. All three tag types are supported, mark-recapture tags
##'     included.
##' }
##'
##' Seasonality is opt-in per model component. Setting `conf$seasonal_spline`
##' gives `alpha` and `beta` a seasonal third dimension, but only the taxis
##' coefficients are estimated season by season. The seasons of the advection
##' coefficients are set separately by `conf$n_seasons_adv`. Diffusion
##' is a second-moment quantity: estimating it separately per season divides the
##' information available per season and is easily traded off against a seasonal
##' taxis acting on the same covariate. Set `conf$seasonal_dif <- TRUE` to
##' estimate season-specific diffusion instead. That frees the diffusion
##' intercept (`beta[1, 1, ]`, the overall diffusion level) as well as the
##' covariate-dependent coefficients, so it also has an effect when diffusion is
##' a single estimated constant.
##'
##' Because the seasonal dimension is shared by all covariates but each
##' covariate is indexed through its own `dat$time_spline[[i]]`, a covariate with
##' fewer breaks than the maximum never evaluates its trailing slices. Those
##' coefficients are fixed rather than estimated, since they would otherwise be
##' free parameters that do not enter the likelihood.
##'
##' @param dat A data list as produced by [setup_data()].
##' @param conf A configuration list as produced by [default_conf()].
##' @param par A parameter list with initial values as produced by
##'   [default_par()].
##'
##' @return
##' A named list of factors with elements `alpha`, `logKappa`, `beta`, `gamma`
##' (when the data have advection fields), `adv_const` and `logSdO`. Entries with `NA` are fixed, while equal factor levels are
##' estimated as the same parameter.
##'
##' @examples
##' map <- with(skjepo$sim, default_map(dat, conf, par))
##'
##' @export
default_map <- function(dat, conf, par){

  map <- list()

  ## Seasonal structure is specified in conf; derive the breakpoints it implies
  res <- .resolve_seasons(dat, conf)
  dat <- res$dat
  conf <- res$conf

  ## Spline slices actually used by each covariate
  nsea <- .get_nsea(dat)

  ## Taxis ----------------------------------------------------
  map$alpha <- .make_alpha_map(par$alpha, nsea)

  ## Taxis scaling --------------------------------------------
  map$logKappa <- factor(NA)

  ## Diffusion ------------------------------------------------
  map$beta <- .make_beta_map(par$beta, nsea,
                             seasonal = isTRUE(conf$seasonal_dif))

  ## Advection ------------------------------------------------
  conf <- .adv_conf(conf)
  if (!is.null(par$gamma)) {
    map$gamma <- if (conf$use_advection) {
      .make_gamma_map(par$gamma, conf$adv_gamma)
    } else factor(rep(NA, length(par$gamma)))
  }
  if (!is.null(par$adv_const)) {
    map$adv_const <- if (conf$use_advection && conf$adv_const) {
      factor(seq_along(par$adv_const))
    } else factor(rep(NA, length(par$adv_const)))
  }

  ## Observation error ----------------------------------------
  ## par$logSdO is a 2 x 3 matrix: rows are the x and y direction, columns are
  ## the tag types in the order dtags, stags, ctags. x and y are coupled within
  ## a tag type, and every tag type that estimates its own observation error
  ## gets its own level. Mark-recapture tags are included: the reported
  ## recapture position carries error too (vessel position, rounding to a
  ## fishing block), and nll.R implements obs_var_type for all three types.
  logSdO <- rep(NA_integer_, 6)
  use_type <- c(isTRUE(conf$use_dtags), isTRUE(conf$use_stags),
                isTRUE(conf$use_ctags))
  next_id <- 1L

  for (i in seq_len(3)) {
    if (use_type[i] &&
          .get_obs_var_type_integer(conf$obs_var_type)[i] %in% c(1, 2)) {
      logSdO[(2 * i - 1):(2 * i)] <- next_id  ## assume equal var in x,y
      next_id <- next_id + 1L
    }
  }

  map$logSdO <- factor(logSdO)

  ## return
  map
}


## Internal functions ---------------------------------------------------------------

## Per-covariate number of spline slices, recycled to the covariate dimension of
## a parameter array and capped at that array's seasonal dimension. Falls back to
## the full seasonal dimension when the two are inconsistent, so a mismatch never
## silently fixes coefficients that the likelihood does use.
.expand_nsea <- function(nsea, ncov, nsea_max) {
  if (is.null(nsea)) nsea <- nsea_max
  nsea <- as.integer(nsea)
  if (length(nsea) == 1L) nsea <- rep(nsea, ncov)
  if (length(nsea) != ncov) nsea <- rep(max(nsea, nsea_max), ncov)
  pmin(pmax(nsea, 1L), nsea_max)
}

.make_alpha_map <- function(alpha, nsea = NULL) {

  dims <- dim(alpha)
  stopifnot(length(dims) == 3)

  nknot <- dims[1]
  ncov <- dims[2]
  nsea <- .expand_nsea(nsea, ncov, dims[3])

  map_id <- array(NA_integer_, dim = dims, dimnames = dimnames(alpha))
  next_id <- 1L

  if (nknot > 1) {
    for (j in seq_len(ncov)) {
      for (s in seq_len(nsea[j])) {
        for (k in 2:nknot) {
          map_id[k, j, s] <- next_id
          next_id <- next_id + 1L
        }
      }
    }
  }

  factor(map_id)
}

.make_beta_map <- function(beta, nsea = NULL, seasonal = FALSE) {

  dims <- dim(beta)
  stopifnot(length(dims) == 3)

  nknot <- dims[1]
  ncov <- dims[2]
  nsea <- .expand_nsea(nsea, ncov, dims[3])

  map_id <- array(NA_integer_, dim = dims, dimnames = dimnames(beta))
  next_id <- 1L

  ## Intercept: carried by a single covariate, since one constant per covariate
  ## is not identifiable. Under seasonal diffusion the intercept is the overall
  ## diffusion level in each season -- the part of a seasonal diffusion that is
  ## still present when diffusion does not vary with the covariate at all -- so
  ## it is carried by the covariate with the most slices. That is covariate 1
  ## whenever all covariates have the same number of slices.
  j0 <- which.max(nsea)
  if (seasonal) {
    for (s in seq_len(nsea[j0])) {
      map_id[1, j0, s] <- next_id
      next_id <- next_id + 1L
    }
  } else {
    map_id[1, j0, seq_len(nsea[j0])] <- next_id
    next_id <- next_id + 1L
  }

  ## remaining spline coefficients: independent across knots and covariates,
  ## and coupled across seasons unless seasonal diffusion is requested
  if (nknot > 1) {
    for (j in seq_len(ncov)) {
      for (k in 2:nknot) {
        if (seasonal) {
          for (s in seq_len(nsea[j])) {
            map_id[k, j, s] <- next_id
            next_id <- next_id + 1L
          }
        } else {
          map_id[k, j, seq_len(nsea[j])] <- next_id
          next_id <- next_id + 1L
        }
      }
    }
  }

  factor(map_id)
}

## gamma is [direction (x, y), field, season]. "shared": one level per field and
## season for both directions; "xy": one per direction as well.
.make_gamma_map <- function(gamma, type = c("shared", "xy")) {

  type <- match.arg(type)
  dims <- dim(gamma)
  stopifnot(length(dims) == 3, dims[1] == 2)

  map_id <- array(NA_integer_, dim = dims, dimnames = dimnames(gamma))
  next_id <- 1L

  for (s in seq_len(dims[3])) {
    for (f in seq_len(dims[2])) {
      if (type == "shared") {
        map_id[, f, s] <- next_id
        next_id <- next_id + 1L
      } else {
        map_id[, f, s] <- next_id + 0:1
        next_id <- next_id + 2L
      }
    }
  }

  factor(map_id)
}
