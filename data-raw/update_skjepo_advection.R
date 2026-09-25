## Bring the stored skjepo objects to the parameter structure of advection
## fields (2026-09-25): gamma only exists with an advection field
## (setup_data(adv = )), and every model has a constant drift adv_const. The
## example has no advection, so gamma is dropped and adv_const added at 0,
## fixed. Nothing is re-simulated or refitted, so the data stay identical.
## make_skjepo.R already produces this structure when rerun.

devtools::load_all()

fix_par <- function(p) {
  if (is.null(p)) return(p)
  p$gamma <- NULL
  if (is.null(p$adv_const)) p$adv_const <- matrix(0, 2, 1)
  p
}
fix_map <- function(m) {
  if (is.null(m)) return(m)
  m$gamma <- NULL
  m$adv_const <- factor(c(NA, NA))
  m
}
fix_conf <- function(conf) .adv_conf(conf)

skjepo$sim$par <- fix_par(skjepo$sim$par)
skjepo$sim$par_true <- fix_par(skjepo$sim$par_true)
skjepo$sim$map <- fix_map(skjepo$sim$map)
skjepo$sim$conf <- fix_conf(skjepo$sim$conf)

skjepo$fit$par <- fix_par(skjepo$fit$par)
skjepo$fit$pl <- fix_par(skjepo$fit$pl)
if (!is.null(skjepo$fit$plsd)) {
  skjepo$fit$plsd$gamma <- NULL
  skjepo$fit$plsd$adv_const <- matrix(0, 2, 1)
}
skjepo$fit$map <- fix_map(skjepo$fit$map)
## the stored habi objects are rebuilt on use (their pointers do not survive
## saving anyway); drop the old per-direction advection entries
skjepo$fit$pred$habi$adv_x <- NULL
skjepo$fit$pred$habi$adv_y <- NULL
skjepo$fit$conf <- fix_conf(skjepo$fit$conf)

usethis::use_data(skjepo, overwrite = TRUE)
