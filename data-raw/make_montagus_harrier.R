## make montagus harrier data set

library(admove)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)


crs_proj <- sf::st_crs(3035)

land <- ne_countries(scale = "medium", returnclass = "sf")
land_proj <- st_transform(land, crs = crs_proj)

is_on_land <- function(x, y, land) {
  pt <- st_sfc(st_point(c(x, y)), crs = st_crs(land))
  any(st_intersects(pt, st_geometry(land), sparse = FALSE))
}

set.seed(123)


## spatial domain ------------------------------------------

## Example for now based on: alerstam_ecology_2018 (Montagu’s harrier)
xlim <- c(-25, 20)
ylim <- c(-5, 65)

xlim <- c(-15, 15)
ylim <- c(5, 60)

bb_ll <- st_bbox(c(xmin = xlim[1], xmax = xlim[2],
                   ymin = ylim[1], ymax = ylim[2]),
                 crs = st_crs(4326))
bb_ll_sfc <- st_as_sfc(bb_ll)



bb_pr <- st_bbox(st_transform(bb_ll_sfc, crs_proj))

xrange_m <- c(bb_pr["xmin"], bb_pr["xmax"])
yrange_m <- c(bb_pr["ymin"], bb_pr["ymax"])


grid <- create_grid(xrange = xrange_m,
                    yrange = yrange_m,
                    cellsize = 5e5,
                    crs = crs_proj)

grid <- scale_sref(grid, scale = 1/1e6)

crs_scale <- sref(grid)$crs_scale


plot_grid(grid, plot_land = TRUE, labels = FALSE)


## time info ---------------------------------------------

origin <- as.POSIXct("2020-01-01", "UTC")
tref <- create_tref(origin, "month")



## covariates ---------------------------------------------

hotspot1 <- function(xy, sigma = 1e8) {
  p_ll <- st_sfc(st_point(c(9, 55)), crs = 4326)
  p_target <- st_transform(p_ll, crs_proj)
  cent <- st_coordinates(p_target) * crs_scale
  res <- exp(-((xy[,1] - cent[1])^2 + (xy[,2] - cent[2])^2) / (2 * sigma^2))
  return(res)
}

hotspot2 <- function(xy, sigma = 1e8) {
  p_ll <- st_sfc(st_point(c(-7.2, 6)), crs = 4326)
  p_target <- st_transform(p_ll, crs_proj)
  cent <- st_coordinates(p_target) * crs_scale
  res <- exp(-((xy[,1] - cent[1])^2 + (xy[,2] - cent[2])^2) / (2 * sigma^2))
  return(res)
}

grid_buff <- add_buffer(grid)

sigma <- 5


bree <- reshape2::acast(data.frame(grid_buff$xygrid,
                                   z = hotspot1(grid_buff$xygrid, sigma = sigma)),
                        x ~ y,
                        value.var = "z")
bree <- admove:::.rescale_cov(bree, zrange = c(1, 20))
ow <- reshape2::acast(data.frame(grid_buff$xygrid,
                                 z = hotspot2(grid_buff$xygrid, sigma = sigma)),
                      x ~ y, value.var = "z")
ow <- admove:::.rescale_cov(ow, zrange = c(1, 20))

image(ow)

cov_list <- list(ow, bree, bree, bree, bree, bree, bree, ow, ow, ow, ow, ow)
dt_cov <- 1/12

env <- prep_cov(c(cov_list, cov_list),
                times = seq(0, 2/dt_cov - dt_cov, 1),
                tref = tref,
                sref = sref(grid))

env

plot_cov(env[,,1:12], plot_land = TRUE)

par <- list(alpha = array(c(0,200,350), dim = c(3, 1, 1)),
            beta = array(log(0.02), dim = c(1,1,1)),
            logKappa = log(0.02),
            logSdO = matrix(log(0.0001),2,3))


knots_tax <- matrix(round(quantile(env, probs = c(0.1,0.5,0.9), na.rm = TRUE)), 3, 1)
knots_dif <- matrix(round(quantile(env, probs = c(0.5), na.rm = TRUE)), 1, 1)


## release locations
df <- data.frame(lon = c(9.009992, 11.762067,
                         10.064680, 9.377790),
                 lat = c(56.345311, 55.517587,
                         57.259584, 55.327421))

pts_sf <- st_as_sf(df, coords = c("lon","lat"), crs = 4326)
pts_proj <- st_transform(pts_sf, crs_proj)
rel_loc <- st_coordinates(pts_proj) * crs_scale

times <- runif(3, 0.4 * 12, 0.5 * 12)

tmp <- expand.grid(t0 = times, i = seq_len(nrow(rel_loc)))

rel_events <- transform(
  tmp,
  x0 = rel_loc[i, 1],
  y0 = rel_loc[i, 2]
)[c("t0", "x0", "y0")]

plot(grid, auto_layout = FALSE, plot_land = TRUE)
points(rel_events[,2:3], col = 2)

sim_list <- sim_tags("d",
                     grid = grid,
                     cov = env,
                     par = par,
                     n_tags = 20,
                     dt_tags = 0.1,
                     trange = c(0,2) * 12,
                     trange_rec = c(1.9,2) * 12,
                     release_events = rel_events,
                     knots_tax = knots_tax,
                     knots_dif = knots_dif,
                     tref = tref,
                     sref = sref(grid))

plot_pref_func(sim_list)

atags0 <- sim_list$tags

summary(atags0)

plot(atags0, plot_land = TRUE)

## plot_tags(atags0, plot_land = TRUE, labels = TRUE, col = "black")


plot(env, auto_layout = FALSE, plot_land = TRUE)
plot(atags0, add = TRUE)

atags <- split(atags0, atags0$id)

atags <- atags[sapply(atags, function(tag) {
  x <- tag[1, "x"] / crs_scale(atags0)
  y <- tag[1, "y"] / crs_scale(atags0)
  is_on_land(x, y, land_proj)
})]

print(paste0("tags on land: ", length(atags)))

atags <- admove:::.add_class(atags, "admove_tags")
atags <- admove:::add_sref(atags, sref(grid))
atags <- admove:::add_tref(atags, tref)


par(mfrow = c(1,1))
plot_tags(atags, leg_pos = "bottomright", plot_land = TRUE,
          auto_layout = FALSE)

range(sapply(atags, function(x) max(x[,3])))
range(sapply(atags, function(x) max(x[,2])))

## Assume mark-resight (ringing for Montagus harrier)
stags <- lapply(atags, function(x) {
  idx <- sort(unique(which(x[,3] >= 2.6)))
  idx <- unique(idx)
  n_sample <- round(runif(1, 10, 15))
  if (length(idx) > n_sample) {
    idx_sample <- sort(sample(idx, n_sample))
  } else {
    idx_sample <- idx
  }
  idx <- sort(unique(which(x[,3] <= -0.5)))
  idx <- unique(idx)
  n_sample <- round(runif(1, 5, 8))
  if (length(idx) > n_sample) {
    idx_sample <- c(idx_sample, sort(sample(idx, n_sample)))
  } else {
    idx_sample <- c(idx_sample, idx)
  }
  x[sort(unique(c(1,idx_sample))), ]
})

stags <- stags[sapply(stags, function(x) nrow(x)) > 1]
stags <- stags[!sapply(stags, function(x) any(is.na(x)))]

## on land
stags <- lapply(stags, function(tag) {
  res <- rep(NA, nrow(tag))
  for(i in 1:nrow(tag)){
    lon <- tag[i, "x"] / crs_scale(atags0)
    lat <- tag[i, "y"] / crs_scale(atags0)
    res[i] <- is_on_land(lon, lat, land_proj)
  }
  tag$tag_type <- rep(admove:::.get_tag_type("s"), nrow(tag))
  tag[res,]
})

stags <- check_tags(stags)
sref(stags) <- sref(grid)
tref(stags) <- tref

sref(stags)
head(stags)

summary(stags)

par(mfrow = c(1,1))
plot_tags(stags, leg_pos = "bottomright",
          plot_land = TRUE, auto_layout = FALSE)

## plot(stags, plot_land = TRUE, labels = TRUE, col = "black")

## testing
## stags <- atags0

head(stags)

rownames(stags) <- NULL

dat <- setup_data(grid = grid,
                  cov = env,
                  tags = stags,
                  knots_tax = knots_tax,
                  knots_dif = knots_dif,
                  transform_sref = TRUE,
                  shift_tref = TRUE)


confi <- default_conf(dat)
confi$obs_var_type[2] <- FALSE
pari <- default_par(dat, confi)
pari$logSdO <- sim_list$par_sim$logSdO
mapi <- default_map(dat, confi, pari)


res <- list()
res$grid <- grid
res$cov <- env
res$par_sim <- sim_list$par_sim
res$tags <- stags
res$dat <- dat
res$conf <- confi
res$par <- pari
## copy kappa as it is fixed
res$par$logKappa <- res$par_sim$logKappa
res$map <- mapi

res <- admove:::.add_class(res, "admove_sim")
res <- admove:::add_sref(res, sref(dat))
res <- admove:::add_tref(res, tref(dat))

montagus_harrier <- res


## save
usethis::use_data(montagus_harrier, overwrite = TRUE)
