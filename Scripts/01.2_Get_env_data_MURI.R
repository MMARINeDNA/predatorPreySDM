####predatorPrey Distribution
####Get environmental data - MURI
####Spring 2026 AVC

library(tidyverse)
library(sdmpredictors)
library(terra)

# set environmental data directory
options(sdmpredictors_datadir = "Data/env-data")

load("ProcessedData/detect_data_muri.Rdata")

#### Filter MURI data for target species ---------------------------------------

detect_data_target <- detect_data_muri %>% 
  filter(BestTaxon %in% c("Lagenorhynchus obliquidens",
                          "Megaptera novaeangliae",
                          "Mirounga angustirostris",
                          "Zalophus californianus",
                          "Phocoena phocoena",
                          "Berardius bairdii"))

detect_data_target %>% group_by(station) %>% n_groups() #177 station
detect_data_target %>% group_by(depth, station) %>% n_groups() #527 station/depths

#### Get environmental data ----------------------------------------------------

#datasets <- list_datasets(terrestrial = FALSE, marine = TRUE)
#bioOracle <- list_layers(datasets[1,1])
#MARSPEC <- list_layers(datasets[2,1])

# MS data
env_dataMS <- load_layers(c("MS_bathy_5m", "MS_biogeo05_dist_shore_5m", 
                            "MS_biogeo06_bathy_slope_5m",
                            "MS_sss07_5m",
                            "MS_sst07_5m"))

names(env_dataMS) <- c("bathy", "distShore", "slope", "SSS07", "SST07")

# BO data
env_dataBO <- load_layers(c("BO2_curvelmean_bdmean",
                            "BO2_chlomean_ss"))

names(env_dataBO) <- c("curVel", "Chla")


# MLD downloaded from de Boyer montégut Clément, https://www.seanoe.org/data/00806/91774/
mld <- rast("Data/MURI/mld_dr003_ref10m_v2023.nc")
mld7 <- mld[["mld_dr003_7"]]

# nepac <- ext(-160, -100, 10, 70)
# mld_np <- crop(mld, nepac)
# plot(mld_np[["mld_dr003_7"]])

#resample to fit resolution
template <- rast(env_dataMS[[1]])
mld7_res <- resample(mld7, template, method = "bilinear")

### Combine all layers ---------------------------------------------------------
env_data <- c(rast(env_dataMS), rast(env_dataBO), mld7_res)

#need upwelling, shelf width, SSH? 
#SSH can get from ROMS, shelf width: calculate distance from shore to 200m isobath for points w/in that space, upwelling mostly covered by SST and salinity?

#### Crop environmental data to survey area

env_data_muri <- raster::crop(env_data, extent(min(detect_data_target$lon_deg),
                                             max(detect_data_target$lon_deg),
                                             min(detect_data_target$lat_deg),
                                             max(detect_data_target$lat_deg)))

data_correlations <- cor(values(env_data_muri), use = "pairwise.complete.obs")
plot_correlation(data_correlations)

#### Merge environmental and detection data using nearest neighbor -------------

env_df_muri <- as.data.frame(env_data_muri, xy=TRUE) %>% 
  mutate(xlon = x, ylat = y) %>% 
  st_as_sf(coords = c("xlon", "ylat"), crs = 4326) %>% 
  st_transform(32610) 

detect_data_merge_muri <- detect_data_target %>% 
  st_join(env_df_muri, join = st_nearest_feature) 

### Save env data --------------------------------------------------------------

save(env_data_muri, env_df_muri, detect_data_target,
     detect_data_merge_muri,
     file = "ProcessedData/detect_and_env_muri.Rdata")

#### Env data research ---------------------------------------------------------
#Summary: shelf width, SST, SSH, MLD, bathymetry, dist to shore, dist to cape,
#dist to freshwater (or SSS?), upwelling index, lat/lon, Chl-a

# From Becker et al. 2020
# Lags
# LON:LAT + shelf + SST + SSH + MLD
# humpy
# LON:LAT + year + depth + SST + MLD
# bairds
# LON:LAT + depth + MLD + SSH
# 
# From Barlow et al. 2026
# phocoena
# depth, distance to coast, shelf width, distance to cape, distance to estuary,
# SST, CUTI (cumulative upwelling index)

## No SDMs for elephant seals?!?! focus on land
# From Robinson et al. tag data:
# Chl-a boundary (TZCF) seems important, 
# and distance to shore/distance to breeding colonies
# latitude (40-50N), associated with gyre boundary
# temperature inversion ()
# narrow persistent density band post-molt (July through November)

# Similar for california sea lions. focus on land
# From Kuhn and Costa tag data
# SST and upwelling index
