# Plot covariate fields

`plot_cov()` plots one or more covariate fields stored in an
`admove_cov` object or in a higher-level *admove* object that contains
covariate data.

Each time slice is plotted as a separate panel, optionally with contour
lines and land added to the plot.

## Usage

``` r
plot_cov(
  x,
  i = 1,
  select = NULL,
  main = "Covariate fields",
  labels = TRUE,
  plot_land = FALSE,
  auto_layout = TRUE,
  xlab = "x",
  ylab = "y",
  bg = NULL,
  plot_contour = TRUE,
  ...
)

# S3 method for class 'admove_cov'
plot(x, ...)
```

## Arguments

- x:

  An object of class `admove_cov`, or an object containing covariate
  data such as `admove_data`, `admove_sim`, or `admove`.

- i:

  Optional index used when `x` is a list-like covariate object. Default:
  `1`.

- select:

  Optional vector of time-step indices to plot. Default: `NULL`, in
  which case all time steps are plotted.

- main:

  Main title of the plot. Default: `"Covariate fields"`.

- labels:

  Logical; currently reserved for plotting cell labels. Default: `TRUE`.

- plot_land:

  Logical; if `TRUE`, add land to the plot. Default: `FALSE`.

- auto_layout:

  Logical; if `TRUE`, plotting parameters are set automatically and
  restored afterwards. Default: `TRUE`.

- xlab:

  Label for the x-axis. Default: `"x"`.

- ylab:

  Label for the y-axis. Default: `"y"`.

- bg:

  Optional background colour for the plot. Default: `NULL`.

- plot_contour:

  Logical; if `TRUE`, add contour lines. Default: `TRUE`.

- ...:

  Additional graphical arguments passed to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

Invisibly returns `NULL`.

## Details

The function uses the x and y dimension names, if present, as plotting
coordinates. Otherwise, row and column indices are used. Each selected
time slice is displayed with
[`graphics::image()`](https://rdrr.io/r/graphics/image.html), and
optional contours are added with
[`graphics::contour()`](https://rdrr.io/r/graphics/contour.html).

## Examples

``` r
plot_cov(skjepo$cov)

```
