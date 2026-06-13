# Summarise admove data

Summarise data of any `admove` object

## Usage

``` r
summarise_data(object, ...)

# S3 method for class 'admove_data'
summary(object, ...)
```

## Arguments

- object:

  an object of class `admove_data` (created by `setup_data`) or an
  object containing such an object (`admove_data`, `admove_sim`, or
  `admove`).

- ...:

  Additional arguments

## Value

Nothing.

## Examples

``` r
summarise_data(skjepo$sim$dat)
#> <admove_grid>
#>   cells:     89
#>   dims:      11 x 10
#>   cellsize:  750.00 x 750.00
#>   xrange:    [-4295.37, 3954.63]
#>   yrange:    [-3632.01, 3867.99]
#>   NAs:       21
#>   units:     km
#> 
#> <admove_cov>
#>   cells:     156
#>   dims:      13 x 12 x 8
#>   cellsize:  750 x 750
#>   xrange:    [-5045.37, 4704.63]
#>   yrange:    [-4382.01, 4617.99]
#>   trange:    [0.00, 21.00]
#>              [2020-01-01 00:00:00
#>                  2021-09-30 18:00:00]
#>   cov range: [20.00, 28.00]
#>   NAs:       35
#>   units:     km x month
#> 
#> <admove_tags>
#>   tags total:     220
#>   ---------------------------------
#>   data-storage tags
#>   n:              20
#>   average over ids:
#>   n obs:          106.60
#>   duration:       10.67
#>   time step:      0.10
#>   x step:         -4.91
#>   y step:         -1.99
#>   ---------------------------------
#>   mark-recapture tags
#>   n:              200
#>   average over ids:
#>   n obs:          2.00
#>   duration:       10.20
#>   time step:      10.20
#>   x step:         -285.87
#>   y step:         -140.05
#>   ---------------------------------
#>   crs:            Azimuthal Equidistant [custom]
#>   datum:          World Geodetic System 1984
#>   crs units:      metre
#>   stored units:   km
#>   crs scale:      1 metre = 0.001 km
#> 
#> 
```
