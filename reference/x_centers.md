# Get x-coordinates of grid-cell centres

Returns the x-coordinates of the cell centres of an `admove_grid`.

## Usage

``` r
x_centers(x)
```

## Arguments

- x:

  A grid object of class `"admove_grid"`, as returned by
  [`create_grid()`](https://tokami.github.io/admove/reference/create_grid.md).

## Value

A numeric vector of x-coordinates for the grid-cell centres.

## Examples

``` r
grid <- create_grid()
xcens <- x_centers(grid)
```
