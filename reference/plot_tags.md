# Plot tagging data

`plot_tags()` plots tagging data in space, showing release locations,
intermediate observations where available, and final recovery or resight
locations.

Depending on the tag type, tags are drawn either as trajectories
connecting successive observations or as straight lines between release
and recovery.

## Usage

``` r
plot_tags(
  x,
  main = "Tags",
  plot_land = FALSE,
  auto_layout = TRUE,
  xlim = NULL,
  ylim = NULL,
  add = FALSE,
  xlab = "x",
  ylab = "y",
  leg_pos = "topright",
  labels = FALSE,
  bg = NULL,
  by_tag_type = TRUE,
  by_tag = FALSE,
  col = c(adjustcolor("grey60", 0.3), .admove_cols(2)),
  pch = c(1, 0, 16),
  cex = 0.8,
  ...
)

# S3 method for class 'admove_tags'
plot(x, ...)
```

## Arguments

- x:

  An object of class `admove_tags`, or an object containing tagging data
  such as `admove_data`, `admove_sim`, or `admove`.

- main:

  Main title of the plot. Default: `"Tags"`.

- plot_land:

  Logical; if `TRUE`, add land to the plot. Default: `FALSE`.

- auto_layout:

  Logical; if `TRUE`, plotting parameters are set automatically and
  restored afterwards. Default: `TRUE`.

- xlim:

  Optional x-axis limits.

- ylim:

  Optional y-axis limits.

- add:

  Logical; if `TRUE`, add to an existing plot. Default: `FALSE`.

- xlab:

  Label for the x-axis. Default: `"x"`.

- ylab:

  Label for the y-axis. Default: `"y"`.

- leg_pos:

  Position of the legend. Default: `"topright"`.

- labels:

  Logical; if `TRUE`, label observations by time instead of plotting
  intermediate points. Default: `FALSE`.

- bg:

  Optional background colour for the plot. Default: `NULL`.

- by_tag_type:

  Logical; if `TRUE`, create separate panels by tag type. Default:
  `TRUE`.

- by_tag:

  Logical; if `TRUE`, create separate panels for individual tags.
  Default: `FALSE`.

- col:

  Character vector of length 1 to 3 giving colours for tag paths,
  release positions, and recovery or final observation positions.

- pch:

  Integer vector of length 1 to 3 giving plotting symbols for
  intermediate observations, release positions, and recovery or final
  observation positions. Default: `c(1, 0, 16)`.

- cex:

  Numeric character expansion factor for plotted points. Default: `0.8`.

- ...:

  Additional graphical arguments passed to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

No return value. Called for its side effect of producing a plot.

## Details

Data-storage and mark-resight tags with intermediate observations are
plotted as trajectories through space. Conventional mark-recapture tags
are plotted as straight-line segments between release and recovery
positions.

When `by_tag_type = TRUE`, separate panels are created for each tag type
present in the data. When `by_tag = TRUE`, each tag is shown in a
separate panel.

## Examples

``` r
plot_tags(skjepo$sim$tags)


```
