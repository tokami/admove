

## Test with lon/lat


## Test with common projection crs


## Test with custom crs without EPSG:
crs_laea <- sf::st_crs("+proj=laea +lat_0=30 +lon_0=-2.5 +datum=WGS84 +units=m +no_defs")
crs_laea$epsg  ## typically NA


## Test with CRS = NA
