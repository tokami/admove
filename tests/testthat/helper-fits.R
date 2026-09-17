## Shared fixtures for tests that need *a* fitted model but do not test the
## estimates themselves (names, labels, knots, summary output, convergence
## bookkeeping). Building the RTMB tape dominates the cost of a fit, and it
## scales with the number of tag positions, so these tests use a small simulated
## data set, and each variant is built once per test run and reused.
##
## Tests that modify a fit work on their own copy: an admove object is a list,
## so e.g. `fit$pl <- NULL` does not affect the cached object. The RTMB 'obj'
## inside is shared, so tests must not change its state.

.admove_test_cache <- new.env(parent = emptyenv())

.cached <- function(key, expr) {
  if (!exists(key, envir = .admove_test_cache, inherits = FALSE)) {
    assign(key, force(expr), envir = .admove_test_cache)
  }
  get(key, envir = .admove_test_cache, inherits = FALSE)
}


## Small simulation on the skjepo example: same grid, covariate and settings as
## skjepo$sim, but 3 archival and 40 conventional tags instead of 20 and 200.
small_sim <- function() {
  .cached("sim", withr::with_seed(1, suppressMessages(suppressWarnings(
    sim_data(skjepo$sim, n_dtags = 3, n_ctags = 40, verbose = FALSE)
  ))))
}


## Fit of small_sim() without predictions. 'sdreport' and 'report' add the
## corresponding post-processing, exactly as admove(do_sdreport = TRUE,
## do_report = TRUE) would.
small_fit <- function(sdreport = FALSE, report = FALSE) {

  fit <- .cached("fit", suppressWarnings(
    admove(small_sim(), do_sdreport = FALSE, do_predictions = FALSE,
           do_report = FALSE, verbose = FALSE)
  ))

  if (sdreport) {
    fit <- .cached("fit_sd", suppressWarnings(add_sdreport(fit)))
  }
  if (report) {
    key <- if (sdreport) "fit_sd_report" else "fit_report"
    fit <- .cached(key, suppressWarnings(add_report(fit)))
  }

  fit
}
