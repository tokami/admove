# Add a one-cell buffer around an admove grid

Expands an `admove_grid` by adding a one-cell buffer around its outer
boundary. Missing cells in the original grid are propagated to the
buffered grid, and boundary cells adjacent to missing regions may also
be excluded.

## Usage

``` r
add_buffer(grid, plot = FALSE)
```

## Arguments

- grid:

  A grid object of class `"admove_grid"`, as returned by
  [`create_grid()`](https://tokami.github.io/admove/reference/create_grid.md).

- plot:

  Logical; if `TRUE`, the buffered grid is plotted.

## Value

An object of class `"admove_grid"` with an additional one-cell buffer
around the original spatial domain.

## Examples

``` r
grid <- create_grid()
grid_with_buffer <- add_buffer(grid)
```
