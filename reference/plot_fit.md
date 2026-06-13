# Summary plots for a fitted `admove` object

`plot_fit()` creates one or several diagnostic summary plots for a
fitted object of class `admove`. Depending on the selected `quantity`,
the function visualises habitat preference functions, taxis, diffusion,
or parameter estimates. Multiple quantities are arranged automatically
in a multi-panel layout.

## Usage

``` r
plot_fit(
  x,
  quantity = c("pref", "taxis", "dif", "par"),
  plot_land = FALSE,
  auto_layout = TRUE,
  col = .admove_cols(10),
  cor_dif = NULL,
  cor_tax = NULL,
  asp = 2,
  plot.legend = 1,
  bg = NULL,
  ...
)

# S3 method for class 'admove'
plot(x, ...)
```

## Arguments

- x:

  A fitted object of class `admove`, as returned by
  [`admove()`](https://tokami.github.io/admove/reference/admove.md).

- quantity:

  Character vector specifying which quantities to plot. Available
  options are:

  `"pref"`

  :   Habitat preference functions.

  `"taxis"`

  :   Taxis (movement direction and magnitude).

  `"dif"`

  :   Diffusion.

  `"par"`

  :   Estimated model parameters.

  Multiple quantities can be selected.

- plot_land:

  Logical; if `TRUE`, land masses are added to spatial plots using
  [`plot_land()`](https://tokami.github.io/admove/reference/plot_land.md).
  Default: `FALSE`.

- auto_layout:

  Logical; if `TRUE`, graphical parameters are set and restored
  automatically, and plots are arranged in a multi-panel layout.
  Default: `TRUE`.

- col:

  Colours used in the plots. Defaults to `.admove_cols(10)`.

- cor_dif:

  Optional scaling factor for diffusion symbols. If `NULL`, a default
  scaling is used internally.

- cor_tax:

  Optional scaling factor for taxis arrows. If `NULL`, a default scaling
  is used internally.

- asp:

  Positive numeric value specifying the target aspect ratio (columns /
  rows) for the plot layout. Default: `2`.

- plot.legend:

  Integer controlling legend placement. If `1` (default), a shared
  legend is drawn in a separate panel below the plots. If `2`, the
  legend is added to the final plot panel.

- bg:

  Optional background colour for the plotting device. If `NULL`
  (default), the current background setting is used.

- ...:

  Additional arguments

## Value

Invisibly returns `NULL`. Called for its side effect of producing plots.

## Details

The function is a wrapper around
[`plot_compare_one()`](https://tokami.github.io/admove/reference/plot_compare_one.md)
and is designed for quick visual inspection of fitted models. If
`auto_layout = TRUE`, the plotting layout is determined automatically
using [`n2mfrow()`](https://rdrr.io/r/grDevices/n2mfrow.html), and
graphical parameters are reset after plotting.

## See also

[`plot_compare_one()`](https://tokami.github.io/admove/reference/plot_compare_one.md)
