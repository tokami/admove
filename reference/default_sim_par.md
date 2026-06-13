# Default parameter values for simulation

Creates a named list of simulation parameter values for `admove`. The
defaults are chosen to give movement rates and observation error that
are scaled to the spatial and temporal extent of the supplied data,
where available.

## Usage

``` r
default_sim_par(
  par = NULL,
  dat = NULL,
  time_unit = NULL,
  alpha_template = c(0, 5, 4),
  target_tax_per_time = NULL,
  target_tax_frac = 1/10,
  target_dif_per_time = NULL,
  target_dif_frac = 1/500,
  target_sdO_frac = 1/30
)
```

## Arguments

- par:

  An optional named list of parameter values used to overwrite the
  corresponding defaults.

- dat:

  An optional data object of class `"admove_data"` used to scale the
  default simulation parameters to the spatial domain, time range, and
  available covariates.

- time_unit:

  Optional time unit label attached to the returned parameter object
  when `dat` is not provided.

- alpha_template:

  Numeric vector giving the template spline coefficients used to
  construct default taxis parameters.

- target_tax_per_time:

  Optional target taxis speed in distance units per time unit. If
  `NULL`, `target_tax_frac` is used instead.

- target_tax_frac:

  Target taxis speed as a fraction of the characteristic spatial scale
  per unit time. Used only if `target_tax_per_time` is `NULL`.

- target_dif_per_time:

  Optional target diffusion coefficient in distance squared per time
  unit. If `NULL`, `target_dif_frac` is used instead.

- target_dif_frac:

  Target diffusion coefficient as a fraction of the squared
  characteristic spatial scale per unit time. Used only if
  `target_dif_per_time` is `NULL`.

- target_sdO_frac:

  Target observation error as a fraction of the characteristic spatial
  scale.

## Value

A named list of simulation parameter values. The returned object also
carries a `"units"` attribute containing the distance and time units,
when available.

## Details

If `dat` is provided, default values are scaled to the spatial extent of
the grid and the total time range of the data. If covariates are
available, the taxis spline coefficients are additionally scaled so that
the resulting taxis strength is approximately consistent with the
requested target movement rate.

The returned parameter list includes defaults for taxis (`alpha`),
diffusion (`beta`), advection (`gamma`), taxis scaling (`logKappa`), and
observation error (`logSdO`). User-supplied values in `par` overwrite
the corresponding defaults.
