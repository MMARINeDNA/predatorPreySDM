####predatorPrey Distribution
####Get environmental data - MURI
####Spring 2026 AVC

library(tidyverse)
library(sdmpredictors)

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
                            "MS_sss08_5m",
                            "MS_sst08_5m"))

names(env_dataMS) <- c("bathy", "distShore", "slope", "SSS08", "SST08")

# BO data
env_dataBO <- load_layers(c("BO2_curvelmean_bdmean",
                            "BO2_chlomean_ss"))

names(env_dataBO) <- c("curVel", "Chla")

#need MLD, upwelling, shelf width, SSH

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
# depth, distace to coast, shelf width, distance to cape, distance to estuary,
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
