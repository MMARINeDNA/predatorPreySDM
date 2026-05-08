###AMBON environmental data
###AVC Jan 2025

### set up environment ---------------------------------------------------------
library(tidyverse)

#spatial modeling
library(sdmpredictors)
library(mregions2)
library(raster)
library(terra)
library(sp)
library(sf)
library(rgdal)
library(viridis)

load("./data products/detect_data.Rdata")

# set environmental data directory
options(sdmpredictors_datadir = "~/env-data")

# Sea ice data 
#downloaded from NCIS Sea Ice Index: https://nsidc.org/data/explore-data

iceShape <- st_union(st_read("./data/metadata/extent_N_202309_polygon/extent_N_202309_polygon_v4.0.shp")) 

### Get environmental data -----------------------------------------------------

shipboard_meta <- read.csv("./data/metadata/FAIRe_noaa-afsc-dbo1.csv") %>% 
  dplyr::select(samp_category, decimalLongitude, decimalLatitude, 
         tot_depth_water_col, temp, chlorophyll, salinity, 
         diss_oxygen, cruise_id, station_id) %>% 
  filter(samp_category == "sample") %>% 
  filter(!grepl("not applicable", station_id)) %>% 
  dplyr::select(-samp_category) %>% 
  mutate(across(1:7, .fns = as.numeric)) %>% 
  group_by(cruise_id,station_id) %>% 
  summarise(across(everything(), mean)) %>% 
  ungroup() %>% 
  mutate(decimalLongitude_360 = ifelse(decimalLongitude < 0,
                                       decimalLongitude + 360,
                                       decimalLongitude))
  
#Also available from ship data: nitrate, nitrite, ammonium, phosphate, pressure, silicate

ship_vect <- vect(shipboard_meta, geom = c("decimalLongitude_360", "decimalLatitude"), 
                  crs = "EPSG:4326")

detect_data_env <- detect_data %>% 
  left_join(shipboard_meta, by = "station_id") %>% 
  mutate(x = decimalLongitude, y = decimalLatitude) %>% 
  st_as_sf(coords = c("x", "y"), crs = 4326) %>% 
  st_transform(st_crs(iceShape))

## we have shipboard data for all stations but ICO2

### Get benthic sp data --------------------------------------------------------

benthic_biomass <- read.csv("./data/metadata/2023 AMBON Epifauna CPUE data.csv",
                            check.names = FALSE) %>% 
  filter(measurementType == "biomass") %>% 
  separate(eventID, into = c(NA, "station_id"), remove = FALSE, sep = "-") %>% 
  dplyr::select(eventDate, eventID, station_id, family,
         measurementValue)

benthic_biomass_family <- benthic_biomass %>% 
  group_by(eventDate, eventID, station_id, 
           family) %>% 
  summarize(measurementValue = mean(measurementValue)) %>% 
  ungroup() %>% 
  pivot_wider(id_cols = c(eventDate, eventID, station_id), names_from = family,
              values_from = measurementValue, values_fill = 0)
  
### Get environmental data -----------------------------------------------------

#datasets <- list_datasets(terrestrial = FALSE, marine = TRUE)
#bioOracle <- list_layers(datasets[1,1])
#MARSPEC <- list_layers(datasets[2,1])

# MS data
env_dataMS <- load_layers(c("MS_bathy_5m", "MS_biogeo05_dist_shore_5m", 
                          "MS_biogeo06_bathy_slope_5m"))

names(env_dataMS) <- c("bathy", "distShore", "slope")

# BO data
env_dataBO <- load_layers(c("BO2_curvelmean_bdmean",
                            "BO2_icethickmean_ss", "BO2_icecovermean_ss"))

names(env_dataBO) <- c("curVel", "iceThick", "iceCov")

# ice raster
iceShape_proj <- st_transform(iceShape, crs(env_dataMS[[1]]))
iceDist_raster <- distance(rast(env_dataMS[[1]]), vect(iceShape_proj))
names(iceDist_raster) <- "iceDist"


# merge with shipboard data
ship_temp <- terra::rasterize(ship_vect, rast(env_dataBO)[[1]], field="temp") 
has_data <- !is.na(ship_temp)
nearest_idx <- terra::which.min(terra::distance(has_data))
###STOP HERE: trying to get nearest neighbor data converting shipboard to raster
ship_temp_filled <- ship_temp
ship_temp_filled[] <- ship_temp[nearest_idx]

ship_temp_filled <- terra::nearest(ship_temp)
ship_chlorophyll <- terra::rasterize(ship_vect, rast(env_dataBO)[[1]], field="chlorophyll")
ship_salinity <- terra::rasterize(ship_vect, rast(env_dataBO)[[1]], field="salinity")
ship_diss_oxygen <- terra::rasterize(ship_vect, rast(env_dataBO)[[1]], field = "diss_oxygen")

env_data <- c(rast(env_dataMS), rast(env_dataBO), 
              iceDist_raster, ship_temp, ship_chlorophyll,
              ship_salinity, ship_diss_oxygen)

names(env_data) <- c("bathy", "distShore", "slope",
                     "curVel", "iceThick", "iceCov", "iceDist",
                     "temp", "chlorophyll", "salinity", "diss_oxygen")

data_extent <- raster::crop(env_data, extent(min(detect_data_env$decimalLongitude),
                                                 max(detect_data_env$decimalLongitude),
                                                 min(detect_data_env$decimalLatitude),
                                                 max(detect_data_env$decimalLatitude)))

data_correlations <- cor(values(data_extent), use = "pairwise.complete.obs")
plot_correlation(data_correlations)

env_df <- as.data.frame(data_extent, xy=TRUE) %>% 
  mutate(xlon = x, ylat = y) %>% 
  st_as_sf(coords = c("xlon", "ylat"), crs = 4326)

iceShape <- iceShape %>% 
  st_crop(st_bbox(detect_data_env))

### Merge satellite and shipboard data -----------------------------------------

env_df_merge <- env_df %>% 
  dplyr::select(-c("temp","chlorophyll","salinity","diss_oxygen")) %>% 
  mutate(lat_factor = as.factor(round(y, digits = 1))) %>% 
  mutate(lon_factor = as.factor(round(x, digits = 1))) %>% 
  group_by(lat_factor, lon_factor) %>% 
  dplyr::select(-geometry) %>% 
  summarise_all(mean) %>% 
  sf::st_drop_geometry()
  

detect_data_merge <- detect_data_env %>% 
  mutate(lat_factor = as.factor(round(decimalLatitude, digits = 1))) %>% 
  mutate(lon_factor = as.factor(round(decimalLongitude, digits = 1))) %>% 
  left_join(env_df_merge, by = c("lat_factor", "lon_factor"))

### Save env data --------------------------------------------------------------

save(data_extent, env_data, env_df, 
     benthic_biomass_family, detect_data_env, detect_data_merge,
     iceShape, file = "env-data/env_data_rasterdf.Rdata")


