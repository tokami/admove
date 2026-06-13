# Plot a summary of simulated admove data

Produces a set of summary plots for an object of class `"admove_sim"`.
Depending on the contents of the simulation object, the plots may
include a covariate field, habitat preference function, taxis,
diffusion, and the simulated tag tracks.

## Usage

``` r
plot_sim(
  x,
  auto_layout = TRUE,
  plot_land = FALSE,
  cor_taxis = NULL,
  cor_diffusion = NULL,
  asp = 2,
  by_tag_type = TRUE,
  ...
)

# S3 method for class 'admove_sim'
plot(x, ...)
```

## Arguments

- x:

  An object of class `"admove_sim"`, as returned by
  [`sim_data()`](https://tokami.github.io/admove/reference/sim_data.md).

- auto_layout:

  Logical; if `TRUE`, the function automatically sets and restores
  graphical parameters.

- plot_land:

  Logical; if `TRUE`, land masses are added to spatial plots.

- cor_taxis:

  Optional scaling factor for taxis arrows.

- cor_diffusion:

  Optional scaling factor for diffusion arrows.

- asp:

  Positive numeric value giving the target aspect ratio (columns / rows)
  of the plot arrangement.

- by_tag_type:

  Logical; if `TRUE`, simulated tracks are plotted separately for each
  tag type. If `FALSE`, all tracks are shown in a single panel.

- ...:

  Additional arguments passed to lower-level plotting functions.

## Value

Invisibly returns `NULL`. The function is called for its plotting side
effects.

## Details

When `auto_layout = TRUE`, the function arranges multiple panels in a
suitable plotting layout. The exact set of panels depends on the
simulation object and whether simulated tag data are available.

If simulated tags are present and `by_tag_type = TRUE`, a separate plot
is produced for each tag type.
