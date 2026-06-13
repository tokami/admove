# Plot predicted location distributions for a single tag

Displays the model's predicted spatial location distribution alongside
the observed track for a single archival tag. The predicted
distributions must be precomputed by
[`add_tag_dist()`](https://tokami.github.io/admove/reference/add_tag_dist.md)
before calling this function; an informative error is raised otherwise.

A multi-panel layout shows a subset of `n` evenly spaced observations.
In each panel the full observed track is drawn in grey, the highlighted
observation in blue, and the predicted density as a colour image.

## Usage

``` r
plot_tag_dist(
  x,
  n = 16,
  asp = 2,
  plot_land = FALSE,
  plot_contour = FALSE,
  xlab = "x",
  ylab = "y"
)
```

## Arguments

- x:

  A fitted object of class `admove` with `$tag_dist` added by
  [`add_tag_dist()`](https://tokami.github.io/admove/reference/add_tag_dist.md).

- n:

  Number of tag observations to display as panels. Default is `16`.

- asp:

  Positive numeric value giving the target aspect ratio (columns / rows)
  for the multi-panel plot arrangement. Default is `2`.

- plot_land:

  Logical; if `TRUE`, land masses are added. Default is `FALSE`.

- plot_contour:

  Logical; if `TRUE`, contour lines are added on top of the predicted
  density image. Default is `FALSE`.

- xlab:

  Label for the x-axis. Default is `"x"`.

- ylab:

  Label for the y-axis. Default is `"y"`.

## Value

Invisibly returns `NULL`. Called for its side effect of producing plots.

## See also

[`add_tag_dist()`](https://tokami.github.io/admove/reference/add_tag_dist.md)
