# Plot spatial habitat preference surfaces

Plot habitat preference in space for taxis or diffusion as a raster-like
surface over the model grid. Preference surfaces can be shown separately
for selected covariates and seasons, or combined across covariates
and/or seasons.

## Usage

``` r
plot_pref_grid(
  x,
  type = "taxis",
  select_cov = NULL,
  select.y = NULL,
  select.sea = NULL,
  combine_cov = FALSE,
  combine.sea = FALSE,
  main = NULL,
  col = hcl.colors(14, "YlOrRd", rev = TRUE),
  ci = 0.95,
  plot_land = FALSE,
  auto_layout = TRUE,
  add = FALSE,
  xlab = "x",
  ylab = "y",
  bg = NULL,
  asp = 2,
  ...
)
```

## Arguments

- x:

  An object of class `admove` or `admove_sim`.

- type:

  Character string specifying which preference surface to plot:
  `"taxis"` (default) or `"diffusion"`.

- select_cov:

  Optional index vector specifying which covariates to plot. If `NULL`,
  all available covariates for the selected `type` are used.

- select.y:

  Optional time slice to use from the covariate field. If `NULL`, the
  first available layer is used.

- select.sea:

  Optional index vector specifying which seasonal components to plot. If
  `NULL`, all available seasons are used.

- combine_cov:

  Logical; if `TRUE`, preference surfaces are summed across selected
  covariates before plotting. Default is `FALSE`.

- combine.sea:

  Logical; if `TRUE`, preference surfaces are summed across selected
  seasonal components before plotting. Default is `FALSE`.

- main:

  Optional main title for the plot panels.

- col:

  Colours for the preference surface. Defaults to
  `hcl.colors(14, "YlOrRd", rev = TRUE)`.

- ci:

  Confidence level for pointwise confidence intervals. Default is
  `0.95`.

- plot_land:

  Logical; if `TRUE`, land masses are added to the plot. Default is
  `FALSE`.

- auto_layout:

  Logical; if `TRUE`, the plotting layout is set automatically. Default
  is `TRUE`.

- add:

  Logical; if `TRUE`, the preference surface is added to an existing
  plot. If `FALSE` (default), a new plot is created.

- xlab:

  Label for the x-axis. Default is `"x"`.

- ylab:

  Label for the y-axis. Default is `"y"`.

- bg:

  Optional background colour for the plotting device. If `NULL`
  (default), the current background setting is used.

- asp:

  Positive numeric value giving the target aspect ratio (columns / rows)
  for multi-panel plot arrangements. Default is `2`.

- ...:

  Additional arguments passed to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) when a new
  plot is created.

## Value

Invisibly returns `NULL`. Called for its side effect of producing plots.

## Details

For fitted objects of class `admove`, the function evaluates the
estimated preference functions on the spatial covariate fields and
displays the resulting preference surface for each selected covariate
and season.

Preference surfaces can be plotted separately or combined across
covariates (`combine_cov = TRUE`) and/or seasons (`combine.sea = TRUE`).

For simulated objects of class `admove_sim`, the preference surface is
reconstructed from the simulated covariates and parameter values.
