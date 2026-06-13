# Compare a single quantity across fitted or simulated admove objects

`plot_compare_one()` provides a compact plotting interface for comparing
a single quantity across multiple fitted or simulated *admove* objects.
Depending on the selected `quantity`, the function compares habitat
preference functions, taxis, diffusion, or parameter estimates.

The input can be supplied either as individual `admove` / `admove_sim`
objects or as a list of such objects.

## Usage

``` r
plot_compare_one(
  fit,
  ...,
  quantity = c("pref", "taxis", "dif", "par"),
  plot_land = FALSE,
  auto_layout = TRUE,
  asp = 2,
  col = .admove_cols(10),
  lty = 1:10,
  cor_tax = NULL,
  cor_dif = NULL,
  plot.legend = 1,
  panel_lab = NULL,
  bg = NULL
)
```

## Arguments

- fit:

  An object of class `admove` or `admove_sim`, or a list of such
  objects.

- ...:

  Additional `admove` or `admove_sim` objects to compare.

- quantity:

  Character string specifying which quantity to compare. Currently
  implemented options are:

  - `"pref"` for habitat preference functions,

  - `"taxis"` for taxis,

  - `"dif"` for diffusion, and

  - `"par"` for parameter estimates.

- plot_land:

  Logical; if `TRUE`, add land masses to spatial plots using
  [`plot_land()`](https://tokami.github.io/admove/reference/plot_land.md).
  Default: `FALSE`.

- auto_layout:

  Logical; if `TRUE`, restore graphical parameters on exit. Default:
  `TRUE`.

- asp:

  Positive numeric value giving the target aspect ratio (columns / rows)
  of the plot arrangement. Default: `2`.

- col:

  Vector of colours used for the different objects being compared. By
  default, colours are taken from `.admove_cols(10)`.

- lty:

  Vector of line types used for the different objects being compared.
  Default: `1:10`.

- cor_tax:

  Optional scaling factor for taxis arrows. If `NULL` (default), the
  longest arrow is automatically scaled to one grid cell width.

- cor_dif:

  Optional scaling factor for diffusion symbols. If `NULL` (default),
  the largest circle is automatically scaled to one grid cell width.

- plot.legend:

  Logical or integer indicating whether, or which, legend should be
  plotted. Default: `1`.

- panel_lab:

  Optional character string added as a panel label in the top-left
  corner of the plot. Default: `NULL`.

- bg:

  Optional background colour for the plotting device. Default: `NULL`.

## Value

No return value. Called for its side effect of producing a comparison
plot.

## Details

For `quantity = "pref"`, habitat preference curves are overlaid in a
single panel. For `quantity = "taxis"` and `quantity = "dif"`, spatial
movement components are compared on the prediction grid. For
`quantity = "par"`, fitted parameter estimates are shown, with
confidence intervals for fitted `admove` objects where available.

Simulated objects of class `admove_sim` can be included in the
comparison, allowing direct visual comparison between simulated and
fitted quantities.
