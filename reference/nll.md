# Negative log-likelihood for the admove model

Computes the negative log-likelihood of the admove movement model given
parameter values and input data.

## Usage

``` r
nll(par, dat)
```

## Arguments

- par:

  A named list of model parameters. The structure should match the
  output of
  [`default_par()`](https://tokami.github.io/admove/reference/default_par.md).

- dat:

  A list containing model input data and configuration settings, as
  returned by
  [`setup_data()`](https://tokami.github.io/admove/reference/setup_data.md)
  and
  [`default_conf()`](https://tokami.github.io/admove/reference/default_conf.md),
  respectively.

## Value

A numeric value representing the negative log-likelihood.

## Details

Two estimation engines are currently implemented:

- **Kalman filter (KF):** continuous space, discrete time formulation

- **Continuous-time Markov chain (CTMC):** discrete space, continuous
  time formulation

Both engines operate on the same input data structure. The choice of
engine is controlled via `conf$engine`, where `1` selects the Kalman
filter and `2` selects the CTMC approach.
