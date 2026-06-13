# Summarise covariate fields

`summarise_cov()` prints a compact summary of one or more covariate
fields. It works on objects of class `admove_cov`, `admove_cov_list`, or
on higher-level *admove* objects that contain covariate data.

The summary includes dimensions, cell size, spatial and temporal ranges,
covariate value range, number of missing cells, and the associated
spatial and temporal units where available.

## Usage

``` r
summarise_cov(object, ...)

# S3 method for class 'admove_cov'
summary(object, ...)

# S3 method for class 'admove_cov'
summary(object, ...)
```

## Arguments

- object:

  An object of class `admove_cov`, `admove_cov_list`, or an object
  containing covariate data such as `admove_data`, `admove_sim`, or
  `admove`.

- ...:

  Additional arguments

## Value

Invisibly returns the corresponding covariate object, coerced internally
to a covariate list if needed.

## Examples

``` r
summarise_cov(skjepo$cov)
#> <admove_cov>
#>   cells:     156
#>   dims:      13 x 12 x 8
#>   cellsize:  750 x 750
#>   xrange:    [-5045.37, 4704.63]
#>   yrange:    [-4382.01, 4617.99]
#>   trange:    [0.00, 7.00]
#>              [2020-01-01 00:00:00
#>                  2021-09-30 18:00:00]
#>   cov range: [20.00, 28.00]
#>   NAs:       35
#>   units:     km x quarter
```
