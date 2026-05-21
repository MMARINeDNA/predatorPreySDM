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

load("ProcessedData/detect_data_ambon.Rdata")

# set environmental data directory
options(sdmpredictors_datadir = "~/env-data")

# Sea ice data 
#downloaded from NCIS Sea Ice Index: https://nsidc.org/data/explore-data

iceShape <- st_union(st_read("./Data/AMBON/extent_N_202309_polygon/extent_N_202309_polygon_v4.0.shp")) 

### Get environmental data -----------------------------------------------------

shipboard_meta <- read.csv("./Data/AMBON/FAIRe_noaa-afsc-dbo1.csv") %>% 
  dplyr::select(samp_name, samp_category, decimalLongitude, decimalLatitude, 
         tot_depth_water_col, temp, chlorophyll, salinity, 
         diss_oxygen, cruise_id, station_id, eventDate) %>% 
  filter(samp_category == "sample") %>% 
  filter(!grepl("not applicable", station_id)) %>% 
  dplyr::select(-samp_category) %>% 
  mutate(across(2:8, .fns = as.numeric), eventDate = as.Date(
    parse_date_time(
      eventDate,
      orders = c("ymd HMS", "mdy")))) %>% 
  group_by(cruise_id,station_id) %>% 
  summarise(across(where(is.numeric), mean, na.rm = TRUE),
            eventDate = mean(eventDate, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(decimalLongitude_360 = ifelse(decimalLongitude < 0,
                                       decimalLongitude + 360,
                                       decimalLongitude)) %>% 
  mutate(month = month(eventDate)) %>% 
  mutate(year = year(eventDate))
  
#Also available from ship data: nitrate, nitrite, ammonium, phosphate, pressure, silicate

detect_data_env <- detect_data %>% 
  left_join(shipboard_meta, by = "station_id") %>% 
  mutate(x = decimalLongitude, y = decimalLatitude) %>% 
  st_as_sf(coords = c("x", "y"), crs = 4326) %>% 
  st_transform(st_crs(iceShape))

## we have shipboard data for all stations but ICO2

### Get benthic sp data --------------------------------------------------------

# benthic_biomass <- read.csv("Data/AMBON/2023 AMBON Epifauna CPUE data.csv",
#                             check.names = FALSE) %>% 
#   filter(measurementType == "biomass") %>% 
#   separate(eventID, into = c(NA, "station_id"), remove = FALSE, sep = "-") %>% 
#   dplyr::select(eventDate, eventID, station_id, family,
#          measurementValue)
# 
# benthic_biomass_family <- benthic_biomass %>% 
#   group_by(eventDate, eventID, station_id, 
#            family) %>% 
#   summarize(measurementValue = mean(measurementValue)) %>% 
#   ungroup() %>% 
#   pivot_wider(id_cols = c(eventDate, eventID, station_id), names_from = family,
#               values_from = measurementValue, values_fill = 0)
  
### Get environmental data -----------------------------------------------------

#datasets <- list_datasets(terrestrial = FALSE, marine = TRUE)
#bioOracle <- list_layers(datasets[1,1])
#MARSPEC <- list_layers(datasets[2,1])

# MS data
env_dataMS <- load_layers(c("MS_bathy_5m", "MS_biogeo05_dist_shore_5m", 
                          "MS_biogeo06_bathy_slope_5m", "MS_sst08_5m", "MS_sst09_5m",
                          "MS_sss08_5m", "MS_sss09_5m"))

names(env_dataMS) <- c("bathy", "distShore", "slope", "MS_sst8", "MS_sst9",
                       "MS_sss8","MS_sss9")

# BO data
env_dataBO <- load_layers(c("BO2_curvelmean_bdmean",
                            "BO2_icethickmean_ss", "BO2_icecovermean_ss",
                            "BO21_chlomean_ss", "BO_chlomean", "BO_dissox"))

names(env_dataBO) <- c("curVel", "iceThick", "iceCov", "BO_meanChl_ss", "BO_meanChl", "BO_O2")

# ice raster
iceShape_proj <- st_transform(iceShape, crs(env_dataMS[[1]]))
iceDist_raster <- distance(rast(env_dataMS[[1]]), vect(iceShape_proj))
names(iceDist_raster) <- "iceDist"

env_data <- c(rast(env_dataMS), rast(env_dataBO), 
              iceDist_raster)

# merge satellite data with shipboard data (detections and environment)
data_extent <- raster::crop(env_data, extent(min(detect_data_env$decimalLongitude),
                                                 max(detect_data_env$decimalLongitude),
                                                 min(detect_data_env$decimalLatitude),
                                                 max(detect_data_env$decimalLatitude)))

data_correlations <- cor(values(data_extent), use = "pairwise.complete.obs")
plot_correlation(data_correlations)

env_df <- as.data.frame(data_extent, xy=TRUE) %>% 
  mutate(xlon = x, ylat = y) %>% 
  st_as_sf(coords = c("xlon", "ylat"), crs = 4326) %>% 
  st_transform(32610) 

iceShape <- iceShape %>% 
  st_crop(st_bbox(detect_data_env))

### Merge satellite and shipboard data -----------------------------------------

detect_data_merge <- detect_data_env %>% 
  st_transform(32610) %>% 
  st_join(env_df, join = st_nearest_feature) %>% 
  mutate(MS_sst = case_when(month>= 9~MS_sst9,
                            month==8~MS_sst8,
                            TRUE~NA),
         MS_sss = case_when(month>= 9~MS_sss9,
                            month==8~MS_sss8,
                            TRUE~NA))

### Save env data --------------------------------------------------------------

save(data_extent, env_data, env_df, 
     detect_data_env, detect_data_merge,
     iceShape, file = "ProcessedData/detect_and_env_ambon.Rdata")


