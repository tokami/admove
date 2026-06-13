# Plot an admove grid

Plots the spatial grid of an `admove_grid` object or of an `admove`
object containing a grid.

## Usage

``` r
plot_grid(
  x,
  main = "Grid",
  labels = TRUE,
  plot_grid = TRUE,
  plot_land = FALSE,
  plot_bg = TRUE,
  auto_layout = TRUE,
  xlab = "x",
  ylab = "y",
  bg = NULL,
  ...
)

# S3 method for class 'admove_grid'
plot(x, ...)
```

## Arguments

- x:

  A grid object of class `"admove_grid"` or an object containing a grid,
  such as `"admove_data"`, `"admove_sim"`, or `"admove"`.

- main:

  Main title for the plot. Default is `"Grid"`.

- labels:

  Logical; if `TRUE`, cell numbers are plotted at the cell centres.

- plot_grid:

  Logical; if `TRUE`, grid lines are added.

- plot_land:

  Logical; if `TRUE`, land masses are added to the plot.

- plot_bg:

  Logical; if `TRUE`, active cells are shaded in the background.

- auto_layout:

  Logical; if `TRUE`, graphical parameters are set automatically and
  restored on exit.

- xlab:

  Label for the x-axis.

- ylab:

  Label for the y-axis.

- bg:

  Optional background colour for the plotting device.

- ...:

  Additional arguments passed to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

Invisibly returns `NULL`. The function is called for its plotting side
effects.

## Examples

``` r
plot_grid(skjepo$grid)

```
