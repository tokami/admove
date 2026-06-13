# Compare fitted and simulated admove objects

Create comparison plots for one or more fitted or simulated `admove`
objects. Depending on `quantity`, the function can compare habitat
preference, taxis, diffusion, or parameter estimates across objects.
Multiple requested quantities are arranged automatically in a
multi-panel layout.

## Usage

``` r
plot_compare(
  fit,
  ...,
  quantity = c("pref", "taxis", "dif", "par"),
  plot_land = FALSE,
  auto_layout = TRUE,
  col = .admove_cols(10),
  cor_dif = NULL,
  cor_tax = NULL,
  asp = 2,
  plot.legend = 1,
  bg = NULL
)
```

## Arguments

- fit:

  Either a single object of class `admove` or `admove_sim`, or a list of
  such objects. If a named list is supplied, names are used in the
  legend.

- ...:

  Additional `admove` or `admove_sim` objects to compare.

- quantity:

  Character vector specifying which quantities to compare. Implemented
  options are:

  `"pref"`

  :   Habitat preference.

  `"taxis"`

  :   Taxis.

  `"dif"`

  :   Diffusion.

  `"par"`

  :   Parameter estimates.

  Multiple quantities can be selected.

- plot_land:

  Logical; if `TRUE`, land masses are added to spatial plots using
  [`plot_land()`](https://tokami.github.io/admove/reference/plot_land.md).
  Default is `FALSE`.

- auto_layout:

  Logical; if `TRUE`, the plot layout and graphical parameters are set
  automatically. Default is `TRUE`.

- col:

  Colours used for the different objects being compared. Defaults to
  `admove:::.admove_cols(10)`.

- cor_dif:

  Optional scaling factor for diffusion symbols. If `NULL`, the default
  internal scaling is used.

- cor_tax:

  Optional scaling factor for taxis arrows. If `NULL`, the default
  internal scaling is used.

- asp:

  Positive numeric value giving the target aspect ratio (columns / rows)
  for the plot arrangement. Default is `2`.

- plot.legend:

  Logical or integer controlling legend placement. If set to `1`
  (default), a shared legend is drawn in a separate layout panel. If set
  to `2`, the legend is added within the final plot panel.

- bg:

  Optional background colour for the plotting device. If `NULL`
  (default), the current background setting is used.

## Value

Invisibly returns `NULL`. Called for its side effect of producing plots.

## Details

If `auto_layout = TRUE`, the function arranges the requested comparison
plots automatically. For spatial quantities, land can optionally be
added via
[`plot_land()`](https://tokami.github.io/admove/reference/plot_land.md).
Simulated objects of class `admove_sim` can be compared directly with
fitted objects of class `admove`.

When `plot.legend = 1`, a shared legend is drawn below the plots. If the
input objects are unnamed, fitted objects are labelled sequentially and
simulated objects are labelled `"Sim"`.
