# Plot habitat preference functions

Plot estimated or simulated habitat preference functions against
covariate values for taxis or diffusion. The function can display one or
several covariate-specific preference functions, optionally with
confidence bands for fitted `admove` objects.

## Usage

``` r
plot_pref_func(
  x,
  type = "taxis",
  select = NULL,
  main = NULL,
  cols = .admove_cols(10),
  lwd = 1,
  ci = 0.95,
  auto_layout = TRUE,
  add = FALSE,
  panel_lab = NULL,
  xlab = NULL,
  ylab = "Preference",
  bg = NULL,
  ylim = NULL,
  xlim = NULL,
  return_limits = FALSE,
  data.range = FALSE,
  asp = 2,
  leg_ncol = 1,
  ...
)
```

## Arguments

- x:

  An object of class `admove` or `admove_sim`.

- type:

  Character string specifying which preference function to plot:
  `"taxis"` (default) or `"diffusion"`.

- select:

  Optional index vector specifying which covariates to plot. If `NULL`,
  all available covariates for the selected `type` are shown.

- main:

  Optional main title. Can be a single character string or a character
  vector with one title per panel. If `NULL`, covariate names are used
  where available.

- cols:

  Colours used for the plotted preference functions. Defaults to
  `admove:::.admove_cols(10)`.

- lwd:

  Line width. Default is `1`.

- ci:

  Confidence level for pointwise confidence intervals. Default is
  `0.95`.

- auto_layout:

  Logical; if `TRUE`, the plotting layout is set automatically. Default
  is `TRUE`.

- add:

  Logical; if `TRUE`, the preference function is added to an existing
  plot. If `FALSE` (default), a new plot is created.

- panel_lab:

  Optional character string added as a panel label in the top-left
  corner of the plot. Default is `NULL`.

- xlab:

  Label for the x-axis. If `NULL` (default), the covariate name from
  `dat$cov` is used for each panel; falls back to `"Covariate"` when no
  name is available.

- ylab:

  Label for the y-axis. Default is `"Preference"`.

- bg:

  Optional background colour for the plotting device. If `NULL`
  (default), the current background setting is used.

- ylim:

  Optional y-axis limits.

- xlim:

  Optional x-axis limits.

- return_limits:

  Logical; if `TRUE`, no plot is produced and a list with `xlim` and
  `ylim` is returned instead. Default is `FALSE`.

- data.range:

  Logical; if `TRUE`, x-axis limits are based on the range of the
  observed covariate data rather than the prediction grid. Default is
  `FALSE`.

- asp:

  Positive numeric value giving the target aspect ratio (columns / rows)
  for multi-panel plot arrangements. Default is `2`.

- leg_ncol:

  Number of columns in the legend when seasonal curves are shown.
  Default is `1`.

- ...:

  Additional arguments passed to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) when a new
  plot is created.

## Value

Invisibly returns `NULL` when plotting. If `return_limits = TRUE`,
returns a list with components `xlim` and `ylim`.

## Details

For fitted objects of class `admove`, the function plots estimated
preference functions based on predicted values stored in the fitted
object. If standard deviations are available, pointwise confidence bands
are drawn.

For simulated objects of class `admove_sim`, the function reconstructs
the corresponding preference function directly from the simulated spline
knots and coefficients.

If multiple seasonal curves are available for a covariate, they are
shown as separate line types, and a legend is added.
