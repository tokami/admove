# Plot diffusion on a spatial grid

`plot_diffusion()` plots the diffusion component over the spatial
prediction grid for a fitted or simulated `admove` object.

For fitted `admove` objects, the function plots the mean predicted
diffusion across prediction times. For `admove_sim` objects, diffusion
is reconstructed from the simulated covariates and parameter values and
then averaged over prediction times.

## Usage

``` r
plot_diffusion(
  x,
  cor = NULL,
  col = "black",
  alpha = 0.5,
  lwd = 1,
  main = "Diffusion",
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

- cor:

  Optional scaling factor controlling the size of the diffusion symbols.
  If `NULL`, the largest circle is automatically scaled so that its
  diameter equals one grid cell width.

- col:

  Colour of the plotted diffusion symbols. Default: `"black"`.

- alpha:

  Transparency value. Currently not used directly in the plotting call.
  Default: `0.5`.

- lwd:

  Line width used for the plotted symbols. Default: `1`.

- main:

  Main title of the plot. Default: `"Diffusion"`.

- plot_land:

  Logical; if `TRUE`, land masses are added using
  [`plot_land()`](https://tokami.github.io/admove/reference/plot_land.md).
  Default: `FALSE`.

- image_bg:

  Logical; if `TRUE` (default), a colour image of diffusion intensity is
  drawn underneath the circles.

- auto_layout:

  Logical; if `TRUE`, graphical parameters are set and restored
  automatically. Default: `TRUE`.

- add:

  Logical; if `TRUE`, diffusion is added to an existing plot. If `FALSE`
  (default), a new plot is created.

- xlab:

  Label for the x-axis. Default: `"x"`.

- ylab:

  Label for the y-axis. Default: `"y"`.

- xaxt:

  A character specifying the x-axis type, passed to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html). Default:
  `"s"`.

- yaxt:

  A character specifying the y-axis type, passed to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html). Default:
  `"s"`.

- bg:

  Optional background colour for the plotting device. If `NULL`
  (default), the current background setting is used.

- ...:

  Additional graphical arguments passed to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) when a new
  plot is created.

## Value

Invisibly returns `NULL`. Called for its side effect of producing a
plot.

## Details

Diffusion is represented by point sizes on the prediction grid, with
larger symbols indicating higher diffusion. For fitted objects,
diffusion is based on `x$pred$hD`. For simulated objects, diffusion is
reconstructed from the simulation setup using
[`default_sim_funcs()`](https://tokami.github.io/admove/reference/default_sim_funcs.md).
