# Spatial reference for admove objects

A lightweight container that stores the spatial reference information
used throughout admove. It bundles:

- **crs**: a CRS specification (typically WKT, PROJ string, or EPSG)

- **units**: the units of the coordinates stored in the object

- **crs_scale**: conversion factor from CRS-units to stored units

## Usage

``` r
create_sref(crs = NA, units = NA_character_, crs_scale = 1)
```

## Arguments

- crs:

  CRS specification. Can be WKT/PROJ string, EPSG integer, or an
  [`sf::crs`](https://r-spatial.github.io/sf/reference/coerce-methods.html)
  input.

- units:

  Character string describing units of stored coordinates, e.g. `"m"`,
  `"km"`, `"degree"`.

- crs_scale:

  Numeric scalar conversion factor from CRS units to stored units. For
  example, if CRS is meters and stored coordinates are kilometers,
  `crs_scale = 0.001`.

## Value

An object of class `admove_sref`.

## Details

The intention is that *objects store their coordinates in the chosen
display/model units* (e.g. km), while `crs` still defines the underlying
standard CRS (typically meters). The `crs_scale` links the two.

## Examples

``` r
sp <- create_sref(crs = 32631, units = "km", crs_scale = 0.001)
sp
#> $crs
#> [1] "PROJCRS[\"WGS 84 / UTM zone 31N\",\n    BASEGEOGCRS[\"WGS 84\",\n        ENSEMBLE[\"World Geodetic System 1984 ensemble\",\n            MEMBER[\"World Geodetic System 1984 (Transit)\"],\n            MEMBER[\"World Geodetic System 1984 (G730)\"],\n            MEMBER[\"World Geodetic System 1984 (G873)\"],\n            MEMBER[\"World Geodetic System 1984 (G1150)\"],\n            MEMBER[\"World Geodetic System 1984 (G1674)\"],\n            MEMBER[\"World Geodetic System 1984 (G1762)\"],\n            MEMBER[\"World Geodetic System 1984 (G2139)\"],\n            MEMBER[\"World Geodetic System 1984 (G2296)\"],\n            ELLIPSOID[\"WGS 84\",6378137,298.257223563,\n                LENGTHUNIT[\"metre\",1]],\n            ENSEMBLEACCURACY[2.0]],\n        PRIMEM[\"Greenwich\",0,\n            ANGLEUNIT[\"degree\",0.0174532925199433]],\n        ID[\"EPSG\",4326]],\n    CONVERSION[\"UTM zone 31N\",\n        METHOD[\"Transverse Mercator\",\n            ID[\"EPSG\",9807]],\n        PARAMETER[\"Latitude of natural origin\",0,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8801]],\n        PARAMETER[\"Longitude of natural origin\",3,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8802]],\n        PARAMETER[\"Scale factor at natural origin\",0.9996,\n            SCALEUNIT[\"unity\",1],\n            ID[\"EPSG\",8805]],\n        PARAMETER[\"False easting\",500000,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8806]],\n        PARAMETER[\"False northing\",0,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8807]]],\n    CS[Cartesian,2],\n        AXIS[\"(E)\",east,\n            ORDER[1],\n            LENGTHUNIT[\"metre\",1]],\n        AXIS[\"(N)\",north,\n            ORDER[2],\n            LENGTHUNIT[\"metre\",1]],\n    USAGE[\n        SCOPE[\"Navigation and medium accuracy spatial referencing.\"],\n        AREA[\"Between 0°E and 6°E, northern hemisphere between equator and 84°N, onshore and offshore. Algeria. Andorra. Belgium. Benin. Burkina Faso. Denmark - North Sea. France. Germany - North Sea. Ghana. Luxembourg. Mali. Netherlands. Niger. Nigeria. Norway. Spain. Togo. United Kingdom (UK) - North Sea.\"],\n        BBOX[0,0,84,6]],\n    ID[\"EPSG\",32631]]"
#> 
#> $units
#> [1] "km"
#> 
#> $crs_scale
#> [1] 0.001
#> 
#> attr(,"class")
#> [1] "admove_sref"
```
