# Get neighbouring cells of an admove grid

Returns the four-neighbour structure of an `admove_grid`, listing for
each active cell the indices of its neighbouring cells above, below,
left, and right.

## Usage

``` r
get_neighbours(grid)
```

## Arguments

- grid:

  A grid object of class `"admove_grid"`, as returned by
  [`create_grid()`](https://tokami.github.io/admove/reference/create_grid.md).

## Value

A matrix with one row per active cell and columns giving the cell index
and the indices of its neighbouring cells.

## Examples

``` r
grid <- create_grid()
neighbours <- get_neighbours(grid)
```
