## make skjepo-inspired data set

library(admove)
library(lubridate)
library(sf)
library(rnaturalearth)


set.seed(123)


## spatial domain ------------------------------------------
ocean <- try(ne_download(scale = 50, type = "ocean",
                     category = "physical",
                     returnclass = "sf"), silent = TRUE)

## for offline use:
if (inherits(ocean, "try-error")) {
  ocean <- admove:::.get_land()
}


bb <- st_bbox(c(xmin = -150, ymin = -30, xmax = -70, ymax = 30),
              crs = 4326)
bb_sfc <- st_as_sfc(bb)

g <- st_make_valid(st_geometry(ocean))
g <- st_union(g)

epo <- st_crop(g, bb_sfc)

## Azimuthal Equidistant centered at (−110, 0)
crs_aeqd <- "+proj=aeqd +lat_0=0 +lon_0=-110 +datum=WGS84 +units=m +no_defs"

epo_proj <- st_transform(epo, crs_aeqd)

plot(st_geometry(epo_proj))


## do not run (requires manual selection)
if (FALSE) {

  grid <- create_grid(epo_proj,
                      cellsize = 750e3,
                      select = 2,
                      plot = TRUE)
  paste(which(is.na(grid$celltable)), collapse = ",")

} else {

  grid <- create_grid(epo_proj,
                      cellsize = 750e3,
                      select = -c(44,55,66,76,77,85,86,87,88,94,95,96,97,98,99,
                                  105,106,107,108,109,110),
                      plot = TRUE)

}

grid_km <- scale_sref(grid, 0.001)

## or: grid_km2 <- scale_coords(grid, units = "km")

plot(grid_km, plot_land = TRUE)


grid <- grid_km


## time info ---------------------------------------------

origin <- as.POSIXct("2020-01-01", "UTC")
units_cov <- "quarter"
units_model <- "month"

tref_model <- create_tref(origin, units_model)
tref_cov <- create_tref(origin, units_cov)


## covariates ---------------------------------------------

## avoid edge effects
grid_buff <- add_buffer(grid)

## quarterly fields
cov <- sim_cov(grid_buff, nt = 8, rho_t = 0.4,
               simple = FALSE,
               tref = tref_cov)

plot_cov(cov, plot_land = TRUE)



## release events ----------------------------------------
trange_rel <- c(0,5) ## months
xrange_rel <- c(-1000, 1000)
yrange_rel <- c(-1000, 1000)
n_release_events <- 5

release_events <- sim_release_events(grid = grid,
                                     trange_rel = trange_rel,
                                     xrange_rel = xrange_rel,
                                     yrange_rel = yrange_rel,
                                     n_release_events = n_release_events)

head(release_events)
nrow(release_events)


plot(grid, labels = FALSE, auto_layout = FALSE)
points(release_events[,1], release_events[,2])



## tags ---------------------------------------------

## for reference
## par <- list(alpha = array(c(0,40,30), dim = c(3, 1, 1)),
##             beta = array(log(200^2), dim = c(1,1,1)),  ## km²/month
##             logSdO = matrix(log(0.1),2,3))

n_ctags <- 200
n_dtags <- 20
## in months:
trange <- c(0,24)
trange_rec <- c(1,24)

target_tax_frac <- 1/5
target_dif_frac <- 1/50
target_sdO_frac <- 1/100


sim_ctags <- sim_tags("c",
                      grid = grid,
                      cov = cov,
                      n_tags = n_ctags,
                      trange = trange,
                      trange_rec = trange_rec,
                      release_events = release_events,
                      tref = tref_model,
                      target_dif_frac = target_dif_frac,
                      target_tax_frac = target_tax_frac,
                      target_sdO_frac = target_sdO_frac)

sim_ctags$par_sim
ctags <- sim_ctags$tags

plot(ctags)

plot(grid, auto_layout = FALSE, plot_land = TRUE)
plot(ctags, add = TRUE)

## plot(ctags)

sim_dtags <- sim_tags("d",
                grid = grid,
                cov = cov,
                n_tags = n_dtags,
                trange = trange,
                dt_tags = 0.1,
                trange_rec = trange_rec,
                release_events = release_events,
                tref = list(origin = origin,
                            units = units_model),
                target_dif_frac = target_dif_frac,
                target_tax_frac = target_tax_frac,
                target_sdO_frac = target_sdO_frac)

sim_dtags$par_sim
dtags <- sim_dtags$tags

plot(grid, auto_layout = FALSE, plot_land = TRUE)
plot(dtags, add = TRUE)

plot(dtags)

## convert tags to raw format
ctags <- ctags[which(apply(ctags[,c("x","y")], 1, function(x) all(!is.na(x)))),]
ctags0 <- ctags

## create lon,lat to demo prep_tags functionality
sp <- sref(ctags0)
pts <- st_as_sf(
  data.frame(
    x = ctags$x / sp$crs_scale,
    y = ctags$y / sp$crs_scale
  ),
  coords = c("x", "y"),
  crs = sp$crs
)
ll <- st_transform(pts, 4326)
ctags[,c("x","y")] <- st_coordinates(ll)


## create numeric dates to demo prep_tags functionality
start_date <- origin(tref_model)
m <- floor(ctags$t)
date.t <- as.Date(start_date) %m+% months(m) + ddays((ctags$t - m) * 30.4375)
ctags$t <- as.numeric(date.t - as.POSIXct("1899-12-30", tz = "UTC"))

## wide format
ctags <- admove:::ctags_list_2_wide(split(ctags, ctags$id))

skjepo_ctags <- data.frame(fish_id = ctags$id,
                           date_time = ctags$t0,
                           species = rep(111, nrow(ctags)),
                           rel_len = c(rnorm(round(nrow(ctags)/1.5),48,4),
                                       rnorm(nrow(ctags) - round(nrow(ctags)/1.5),65,5)
                                       ),
                           rel_lon = ctags$x0,
                           rel_lat = ctags$y0,
                           date_caught = ctags$t1,
                           recap_lon = ctags$x1,
                           recap_lat = ctags$y1)

head(skjepo_ctags)

ctags <- prep_ctags(skjepo_ctags,
                    names = c(t0 = "date_time",
                              t1 = "date_caught",
                              x0 = "rel_lon",
                              x1 = "recap_lon",
                              y0 = "rel_lat",
                              y1 = "recap_lat",
                              id = "fish_id"),
                    date_origin = "1899-12-30",
                    sref = list(crs = 4326),
                    tref = tref_model)

ctags <- add_sref(ctags, grid, transform_crs = TRUE)

plot(grid, auto_layout = TRUE)
plot_tags(ctags, plot_land = TRUE, add = TRUE)

crs_scale(ctags)

dtags0 <- dtags
sp <- sref(dtags0)
start_date <- origin(dtags0)
dtags2 <- split(dtags, dtags$id)


skjepo_dtags <- lapply(dtags2, function(x){
  x <- x[,1:3]
  x <- x[which(apply(x[,c("x","y")], 1, function(x) all(!is.na(x)))),]
  pts <- st_as_sf(
    data.frame(
      x = x$x / sp$crs_scale,
      y = x$y / sp$crs_scale
    ),
    coords = c("x", "y"),
    crs = sp$crs
  )
  ll <- st_transform(pts, 4326)
  x[,c("x","y")] <- st_coordinates(ll)
  m <- floor(x$t)
  date.t <- as.Date(start_date) %m+% months(m) + ddays((x$t - m) * 30.4375)
  x$t <- as.numeric(date.t - as.POSIXct("1899-12-30", tz = "UTC"))
  x$varlon <- c(0, rbeta(nrow(x)-2, 3, 100), 0)
  x$varlat <- c(0, rbeta(nrow(x)-2, 2, 10), 0)
  colnames(x) <- c("time","mptlon","mptlat","varlon","varlat")
  return(as.data.frame(x))
})

dtags <- prep_dtags(skjepo_dtags,
                    names = c(t = "time", x = "mptlon", y = "mptlat"),
                    date_origin = "1899-12-30",
                    sref = list(crs = 4326))


dtags <- add_sref(dtags, grid, transform_crs = TRUE)
dtags <- add_tref(dtags, tref_model, shift_origin = TRUE)

sref(dtags)


plot_tags(dtags, plot_land = TRUE)

skjepo_cov <- as.array(cov)
class(skjepo_cov) <- "array"


cov <- prep_cov(skjepo_cov)

cov <- add_sref(cov, grid)

plot_cov(cov[,,1:4], plot_land = TRUE,
         xlab = "lon", ylab = "lat")

dat <- setup_data(grid = grid,
                  cov = cov,
                  tags = c(dtags, ctags),
                  tref = tref_model,
                  transform_sref = TRUE,
                  shift_tref = TRUE)

stopifnot(sim_ctags$par_sim$alpha == sim_dtags$par_sim$alpha)
stopifnot(sim_ctags$par_sim$beta == sim_dtags$par_sim$beta)

## simulation list (admove_sim)
sim <- list()
sim$grid <- grid
sim$cov <- cov
sim$par_sim <- sim_ctags$par_sim ## same as: sim_dtags$par_sim
sim$tags <- c(dtags, ctags)
sim$dat <- dat
sim$conf <- default_conf(dat)
sim$par <- default_par(dat, sim$conf)
## copy kappa as it is fixed
sim$par$logKappa <- sim$par_sim$logKappa
sim$map <- default_map(dat, sim$conf, sim$par)

## fix sdO
## sim$par$logSdO <- sim$par_sim$logSdO
## sim$map$logSdO <- factor(rep(NA, 6))

sim <- admove:::.add_class(sim, "admove_sim")
sim <- add_sref(sim, dat)
sim <- add_tref(sim, dat)


plot_sim(sim, cor_diffusion = 0.01, plot_land = TRUE)

sim$tags

## Fit
fit <- admove(sim)

sim$par_sim

summary(fit)

plot_compare(fit, sim, cor_dif = 0.01, cor_tax = 1)


## keep together
skjepo <- list(sim = sim,
               grid = grid,
               cov = cov,
               ctags = skjepo_ctags,
               dtags = skjepo_dtags,
               fit = fit)

## save
usethis::use_data(skjepo, overwrite = TRUE)
