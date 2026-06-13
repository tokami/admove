# Log-likelihood of a fitted admove model

Extracts the log-likelihood from a fitted `admove` model. The result can
be passed to [`stats::AIC()`](https://rdrr.io/r/stats/AIC.html) or
[`stats::BIC()`](https://rdrr.io/r/stats/AIC.html), or used for
likelihood ratio tests.

## Usage

``` r
# S3 method for class 'admove'
logLik(object, ...)
```

## Arguments

- object:

  A fitted model object of class `"admove"`.

- ...:

  Currently unused.

## Value

An object of class `"logLik"` with attributes `df` (number of estimated
parameters) and `nobs` (number of used observations).

## Examples

``` r
fit <- admove(skjepo$sim, do_sdreport = FALSE)
#> Setting tref$origin on object (was NA).
#> Building the model, that can take a few minutes.
#> Model built (0.25min). Minimizing neg. loglik.
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
#> Minimisation done (0.024min). Model converged.
#> Predicting movement rates.
#> Predictions done (0.00058min).
#> Reporting variables.
#> Reporting done (0.077min).
logLik(fit)
#> 'log Lik.' -30610.1 (df=3)
AIC(fit)
#> [1] 61226.2
BIC(fit)
#> [1] 61243.71
```
