# Plot taxis on a spatial grid

Plot the taxis component as arrows over the spatial prediction grid for
a fitted or simulated `admove` object. The function can display taxis at
selected time steps or the average taxis across multiple time steps.

## Usage

``` r
plot_taxis(
  x,
  select = NULL,
  average = TRUE,
  cor = NULL,
  col = "black",
  alpha = 0.5,
  lwd = 1,
  main = "Taxis",
  plot_land = FALSE,
  image_bg = TRUE,
  auto_layout = TRUE,
  add = FALSE,
  xlab = "x",
  ylab = "y",
  xaxt = "s",
  yaxt = "s",
  bg = NULL,
  ...
)
```

## Arguments

- x:

  An object of class `admove` or `admove_sim`.

- select:

  Optional index vector specifying which prediction time steps to plot.
  If `NULL`, all available prediction time steps are used.

- average:

  Logical; if `TRUE` (default), the taxis vectors are averaged over the
  selected time steps. If `FALSE`, taxis is plotted separately for each
  selected time step.

- cor:

  Optional scaling factor for arrow lengths. If `NULL`, the longest
  arrow is automatically scaled to one grid cell width.

- col:

  Colour of the arrows. Default is `"black"`.

- alpha:

  Transparency value. Currently not used directly in the plotting call.
  Default is `0.5`.

- lwd:

  Line width of the arrows. Default is `1`.

- main:

  Main title of the plot. Default is `"Taxis"`.

- plot_land:

  Logical; if `TRUE`, land masses are added using
  [`plot_land()`](https://tokami.github.io/admove/reference/plot_land.md).
  Default is `FALSE`.

- image_bg:

  Logical; if `TRUE` (default), a colour image of taxis magnitude is
  drawn underneath the arrows.

- auto_layout:

  Logical; if `TRUE`, the plotting layout is set automatically. If
  multiple time steps are plotted and `average = FALSE`, panels are
  arranged using
  [`n2mfrow()`](https://rdrr.io/r/grDevices/n2mfrow.html). Default is
  `TRUE`.

- add:

  Logical; if `TRUE`, taxis arrows are added to an existing plot. If
  `FALSE` (default), a new plot is created.

- xlab:

  Label for the x-axis. Default is `"x"`.

- ylab:

  Label for the y-axis. Default is `"y"`.

- xaxt:

  A character specifying the x-axis type, passed to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html). Default is
  `"s"`.

- yaxt:

  A character specifying the y-axis type, passed to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html). Default is
  `"s"`.

- bg:

  Optional background colour for the plotting device. If `NULL`
  (default), the current background setting is used.

- ...:

  Additional arguments passed to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) when a new
  plot is created.

## Value

Invisibly returns `NULL`. Called for its side effect of producing a
plot.

## Details

For objects of class `admove`, the function plots predicted taxis from
`x$pred$hTdx` and `x$pred$hTdy`. For objects of class `admove_sim`,
taxis is recomputed from the simulated covariates and parameter values.

If `average = TRUE`, the mean taxis over the selected time steps is
plotted. Otherwise, one panel per selected time step is produced unless
`add = TRUE`.
