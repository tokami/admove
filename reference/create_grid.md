# Create or modify a spatial grid for admove

Creates a new spatial grid or derives one from an existing object for
use in `admove`. Grids are used for spatial prediction and are required
for model formulations based on a discrete spatial domain, such as the
continuous-time Markov chain (CTMC) approach.

## Usage

``` r
create_grid(
  x = NULL,
  cellsize = NULL,
  xrange = NULL,
  yrange = NULL,
  select = FALSE,
  crs = NULL,
  units = NULL,
  crs_scale = NULL,
  plot_land = FALSE,
  auto_layout = TRUE,
  plot = FALSE,
  force = FALSE,
  verbose = TRUE
)
```

## Arguments

- x:

  An optional object from which to derive the grid. Supported inputs
  include objects of class `"admove_grid"`, `"admove_cov"`,
  `"admove_tags"`, `"admove_data"`, `"admove_sim"`, `"admove"`, as well
  as `sf`/`sfc` objects, raster objects, and matrices. If `NULL`, a new
  grid is created from the supplied `xrange`, `yrange`, and `cellsize`.

- cellsize:

  Numeric vector giving the grid cell size in the x- and y-direction. If
  a single value is supplied, it is used for both directions. If `NULL`,
  a default value is derived from the spatial extent.

- xrange:

  Numeric vector of length 2 giving the range of the x dimension of the
  spatial domain. Ignored if extracted from `x`.

- yrange:

  Numeric vector of length 2 giving the range of the y dimension of the
  spatial domain. Ignored if extracted from `x`.

- select:

  Controls which cells are retained in the grid. If `FALSE` (default),
  all eligible cells are kept. If `TRUE` or an integer value, cells can
  be selected interactively. If a numeric vector of length greater than
  1 is supplied, it is interpreted as cell indices to retain; if all
  supplied indices are negative, those cells are removed.

- crs:

  Optional coordinate reference system for the grid.

- units:

  Optional spatial units for the grid, for example `"degree"`, `"m"`, or
  `"km"`.

- crs_scale:

  Optional scaling between CRS units and the numeric units used in the
  grid.

- plot_land:

  Logical; if `TRUE`, land masses are added to plots, where supported.

- auto_layout:

  Logical; if `TRUE`, plotting methods may adjust graphical parameters
  automatically.

- plot:

  Logical; if `TRUE`, the resulting grid is plotted.

- force:

  Logical; if `TRUE`, the grid is created even when the implied number
  of grid breaks in x or y exceeds 1000. Use with care.

- verbose:

  Logical; if `TRUE`, informative messages are printed.

## Value

An object of class `"admove_grid"` describing the spatial grid.

## Details

The function can construct a grid from scratch, inherit the extent and
missing-cell structure from an existing `admove` object, or derive the
extent from spatial objects such as `sf` geometries, rasters, or
covariate arrays.

When `x` contains missing spatial cells, these are propagated to the new
grid. For `sf` objects, cells can be restricted to locations inside the
supplied geometry. Optional interactive selection can be used to include
or exclude cells manually.

The returned object stores the grid-cell centres, cell indices, lookup
table, grid breaks, and spatial extent, together with the associated
spatial reference.

## Examples

``` r
grid <- create_grid()
dim(grid)
#> nx ny 
#> 10 10 
bbox_grid(grid)
#> xmin xmax ymin ymax 
#>    0    1    0    1 
```
