# Get coordinate reference system (CRS)

Generic function to extract the coordinate reference system (CRS) from a
supported object.

## Usage

``` r
crs(x, ...)

# Default S3 method
crs(x, ...)

# S3 method for class 'admove_sref'
crs(x, ...)

# S3 method for class 'admove_grid'
crs(x, ...)

# S3 method for class 'admove_cov'
crs(x, ...)

# S3 method for class 'admove_tags'
crs(x, ...)

# S3 method for class 'admove_data'
crs(x, ...)
```

## Arguments

- x:

  An object from which to extract a CRS.

- ...:

  Further arguments passed to methods.

## Value

A CRS representation. For admove objects, this is typically whatever is
stored in the corresponding spatial reference object, for example a WKT
string, EPSG code, or another sf-compatible CRS specification.

## See also

[`sf::st_crs()`](https://r-spatial.github.io/sf/reference/st_crs.html)
