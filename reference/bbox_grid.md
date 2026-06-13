# Return the bounding box of an admove grid

Returns the spatial extent of an `admove_grid` object as a named vector
containing the minimum and maximum x and y coordinates.

## Usage

``` r
bbox_grid(grid)
```

## Arguments

- grid:

  A grid object of class `"admove_grid"`, as returned by
  [`create_grid()`](https://tokami.github.io/admove/reference/create_grid.md).

## Value

A named numeric vector with elements `xmin`, `xmax`, `ymin`, and `ymax`.
