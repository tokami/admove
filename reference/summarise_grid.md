# Summarise an admove grid

Prints a short summary of an `admove_grid` object, including the number
of active cells, grid dimensions, cell size, spatial extent, number of
missing cells, and spatial units.

## Usage

``` r
summarise_grid(object, ...)

# S3 method for class 'admove_grid'
summary(object, ...)
```

## Arguments

- object:

  An object of class `"admove_grid"` or an object containing a grid,
  such as `"admove_data"`, `"admove_sim"`, or `"admove"`.

- ...:

  Additional arguments

## Value

Invisibly returns the extracted `admove_grid` object.

## Details

If `x` is not itself an `admove_grid`, the grid is extracted from the
corresponding component of the supplied object.

## Examples

``` r
summarise_grid(skjepo$grid)
#> <admove_grid>
#>   cells:     89
#>   dims:      11 x 10
#>   cellsize:  750.00 x 750.00
#>   xrange:    [-4295.37, 3954.63]
#>   yrange:    [-3632.01, 3867.99]
#>   NAs:       21
#>   units:     km
```
