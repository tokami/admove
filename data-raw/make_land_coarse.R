## create land for offline plot_land usage

if (!requireNamespace("rnaturalearth", quietly = TRUE)) stop("need rnaturalearth")
if (!requireNamespace("sf", quietly = TRUE)) stop("need sf")

land <- rnaturalearth::ne_download(
  scale = 110, type = "land", category = "physical", returnclass = "sf"
)

land <- sf::st_make_valid(land)

## dissolve to a single geometry
land_geom <- sf::st_union(sf::st_geometry(land))

## create a proper sf object with a geometry column named "geometry"
land <- sf::st_sf(geometry = land_geom, crs = 4326)

land <- sf::st_make_valid(land)

dir.create("inst/extdata", showWarnings = FALSE, recursive = TRUE)
saveRDS(land, "inst/extdata/land_110m.rds")
