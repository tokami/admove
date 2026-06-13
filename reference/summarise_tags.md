# Summarise tagging data

`summarise_tags()` prints a compact summary of tagging data stored in an
object of class `admove_tags` or in a higher-level *admove* object that
contains tagging data.

The summary includes the total number of tags, the number of tags by tag
type, and simple averages such as the number of observations per tag,
tag duration, and average step sizes in time and space.

## Usage

``` r
summarise_tags(object, ...)

# S3 method for class 'admove_tags'
summary(object, ...)
```

## Arguments

- object:

  An object of class `admove_tags`, or an object containing tagging data
  such as `admove_data`, `admove_sim`, or `admove`.

- ...:

  Additional arguments

## Value

Invisibly returns the corresponding `admove_tags` object.

## Examples

``` r
summarise_tags(skjepo$sim$dat$tags)
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
