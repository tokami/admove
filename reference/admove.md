# Fit an admove movement model

Fits the `admove` movement model to tagging data to estimate movement
processes and, where applicable, habitat preference and advection
effects. The model supports different types of tagging data and can be
fitted using either a Kalman filter (continuous space, discrete time) or
a continuous-time Markov chain (CTMC; discrete space, continuous time)
formulation.

## Usage

``` r
admove(
  dat,
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
  ...
)
```

## Arguments

- dat:

  A data list containing model input data, as produced by
  [`setup_data()`](https://tokami.github.io/admove/reference/setup_data.md).

- conf:

  An optional configuration list, typically created by
  [`default_conf()`](https://tokami.github.io/admove/reference/default_conf.md).
  If `NULL`, a default configuration is generated from `dat`.

- par:

  An optional named list of initial parameter values, typically created
  by
  [`default_par()`](https://tokami.github.io/admove/reference/default_par.md).
  If `NULL`, default initial values are generated from `dat` and `conf`.

- map:

  An optional parameter map, typically created by
  [`default_map()`](https://tokami.github.io/admove/reference/default_map.md).
  If `NULL`, a default map is generated from `dat`, `conf`, and `par`.

- engine:

  Optional integer to override `conf$engine`. Use `1` for the Kalman
  filter and `2` for the CTMC formulation. If `NULL`, the value in
  `conf` is used.

- run:

  Logical; if `TRUE` (default), the model is optimized. If `FALSE`, only
  the RTMB objective object is constructed and returned.

- lower:

  Optional lower bounds for optimization. If `NULL`, no explicit lower
  bounds are supplied.

- upper:

  Optional upper bounds for optimization. If `NULL`, no explicit upper
  bounds are supplied.

- rel_tol:

  Relative convergence tolerance passed to
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html). Default is
  `1e-10`.

- do_predictions:

  Logical; if `TRUE` (default), model predictions are computed after
  fitting. If `FALSE`, prediction-related outputs are skipped, and some
  plotting methods may not be available.

- do_sdreport:

  Logical; if `TRUE` (default), `RTMB::sdreport()` is run to obtain
  uncertainty estimates for model parameters and derived quantities.

- do_report:

  Logical; if `TRUE` (default), `obj$report()` is run to extract
  reported RTMB quantities.

- save_covariance:

  Logical; if `TRUE`, the covariance matrix from `RTMB::sdreport()` is
  retained. This may substantially increase memory use.

- dbg:

  Logical; if `TRUE`, the function is run in debugging mode. Default is
  `FALSE`.

- control:

  An optional named list of control settings passed to the optimizer.

- verbose:

  Logical; if `TRUE`, progress messages are printed.

- ...:

  Additional arguments passed to `RTMB::MakeADFun()`.

## Value

A fitted model object of class `"admove"`. Depending on the function
arguments, the returned object may include the RTMB objective function,
the optimization output, reported quantities, predictions, and
uncertainty estimates.

## Details

This is the main model-fitting function in `admove`. It combines the
supplied data, configuration, parameter values, and parameter mapping,
constructs an RTMB objective function, and optionally optimizes it using
[`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html).

If configuration settings, initial parameter values, or parameter maps
are not supplied, they are generated automatically using
[`default_conf()`](https://tokami.github.io/admove/reference/default_conf.md),
[`default_par()`](https://tokami.github.io/admove/reference/default_par.md),
and
[`default_map()`](https://tokami.github.io/admove/reference/default_map.md),
respectively.

## Examples

``` r
fit <- admove(skjepo$sim, do_sdreport = FALSE)
#> Setting tref$origin on object (was NA).
#> Building the model, that can take a few minutes.
#> Model built (0.27min). Minimizing neg. loglik.
#>   0: 2.6621148e+08:  0.00000  0.00000  0.00000
#>   1:     97824176.: 0.0304759 0.0479782 0.998383
#>   2:     58120629.: 0.144985 -0.724747  1.62271
#>   3:     20987545.: 0.0501387 -0.525837  2.59813
#>   4:     11698081.: 0.804236 -0.390874  3.24087
#>   5:     4176101.9: 0.653388 -0.172935  4.20511
#>   6:     2315402.5:  1.38283 -0.0565109  4.87917
#>   7:     899135.28: 0.881767 0.0650439  5.73601
#>   8:     477695.10:  1.00349 -0.518613  6.53883
#>   9:     182171.19: 0.849998 -0.184989  7.46896
#>  10:     113092.53:  1.49966 0.328658  8.02941
#>  11:     58857.934:  1.09352 0.155624  8.92669
#>  12:     44857.604:  1.75497 -0.0553574  9.64640
#>  13:     33307.738:  1.51547 0.211375  10.5799
#>  14:     31866.793:  1.93611  1.08067  10.8395
#>  15:     30696.247:  1.49547 0.690490  11.6480
#> Warning: NA/NaN function evaluation
#>  16:     30669.312:  1.42363 0.758147  11.6641
#>  17:     30619.129:  1.02488 0.728827  11.6525
#>  18:     30612.294: 0.685632 0.525310  11.5934
#>  19:     30610.207: 0.683988 0.526026  11.6453
#>  20:     30610.119: 0.667122 0.508135  11.6376
#>  21:     30610.106: 0.664376 0.512619  11.6350
#>  22:     30610.102: 0.669039 0.516039  11.6361
#>  23:     30610.102: 0.668728 0.515873  11.6366
#>  24:     30610.102: 0.669268 0.516236  11.6365
#>  25:     30610.102: 0.669468 0.516370  11.6365
#>  26:     30610.102: 0.669485 0.516383  11.6365
#> Minimisation done (0.025min). Model converged.
#> Predicting movement rates.
#> Predictions done (0.0011min).
#> Reporting variables.
#> Reporting done (0.086min).
```
