# Plot components of an `admove_data` object

Create summary plots for the main components of an object of class
`admove_data`. Depending on which components are present, the function
plots the spatial grid, covariate fields, and tag data in a multi-panel
layout.

## Usage

``` r
plot_data(x, auto_layout = TRUE, ...)

# S3 method for class 'admove_data'
plot(x, ...)
```

## Arguments

- x:

  An object of class `admove_data`, as returned by
  [`setup_data()`](https://tokami.github.io/admove/reference/setup_data.md).

- auto_layout:

  Logical; if `TRUE`, the plotting layout and graphical parameters are
  set automatically. Default is `TRUE`.

- ...:

  Additional arguments passed to the underlying plotting functions,
  including
  [`plot_grid()`](https://tokami.github.io/admove/reference/plot_grid.md),
  [`plot_cov()`](https://tokami.github.io/admove/reference/plot_cov.md),
  and
  [`plot_tags()`](https://tokami.github.io/admove/reference/plot_tags.md).

## Value

Invisibly returns `NULL`. Called for its side effect of producing plots.

## Details

The function inspects `x` and plots all available data components. If
present, the grid is plotted first, followed by each covariate field,
and then tag data split by tag type:

- `"d"`:

  Archival tags.

- `"s"`:

  Mark-resight tags.

- `"c"`:

  Conventional tags.

If `auto_layout = TRUE`, panels are arranged automatically using
[`n2mfrow()`](https://rdrr.io/r/grDevices/n2mfrow.html). Panel labels
are added with
[`add_lab()`](https://tokami.github.io/admove/reference/add_lab.md).
