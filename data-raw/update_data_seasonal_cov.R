## Move seasonal_cov from the configuration to the data in the stored example
## objects (2026-09-27): whether a covariate field repeats every period is now
## set with setup_data(seasonal_cov = ), and check_conf() rejects
## conf$seasonal_cov. No stored example has a seasonal covariate, so nothing is
## re-simulated or refitted and the data stay identical. make_skjepo.R and
## make_montagus_harrier.R already produce this structure when rerun.

devtools::load_all()

fix <- function(x) {
  if (is.null(x)) return(x)
  ncov <- length(x$dat$cov)
  sc <- x$conf$seasonal_cov
  x$conf$seasonal_cov <- NULL
  if (!is.null(x$dat)) {
    x$dat$seasonal_cov <- if (ncov == 0L) logical(0) else
      rep_len(if (is.null(sc)) FALSE else sc, ncov)
  }
  x
}

skjepo$sim <- fix(skjepo$sim)
skjepo$fit <- fix(skjepo$fit)
montagus_harrier <- fix(montagus_harrier)

usethis::use_data(skjepo, overwrite = TRUE)
usethis::use_data(montagus_harrier, overwrite = TRUE)
