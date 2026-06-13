# Default simulation functions for admove

Constructs default function objects used for simulation in `admove`,
including functions for taxis, diffusion, the spatial gradient of
diffusion, and advection.

## Usage

``` r
default_sim_funcs(dat, conf, par, funcs = NULL)
```

## Arguments

- dat:

  An `admove_data` object, as produced by
  [`setup_data()`](https://tokami.github.io/admove/reference/setup_data.md).

- conf:

  A configuration list, typically produced by
  [`default_conf()`](https://tokami.github.io/admove/reference/default_conf.md).

- par:

  A parameter list, typically produced by
  [`default_par()`](https://tokami.github.io/admove/reference/default_par.md)
  or
  [`default_sim_par()`](https://tokami.github.io/admove/reference/default_sim_par.md).

- funcs:

  An optional named list of user-supplied functions that overwrite the
  corresponding defaults. Allowed names are `"tax"`, `"dif"`, `"ddif"`,
  and `"adv"`.

## Value

A named list of simulation functions with elements `tax`, `dif`, `ddif`,
and `adv`.

## Details

If covariate data are available in `dat`, the default functions are
constructed from the parameter values and interpolated covariate fields.
This allows taxis, diffusion, and advection to vary in space and time.

If no covariate data are available, simple default functions are
returned: diffusion is constant, and taxis, advection, and the gradient
of diffusion are set to zero.

User-supplied functions in `funcs` replace the corresponding default
functions in the returned list.
