#### MURI Distribution w prey variables
#### Summer 2026
#### AVC&MS

#general
library(tidyverse)
library(PNWColors)

#mapping
library(marmap)
library(terra)

#modeling
library(mgcv)

load("./ProcessedData/detect_and_env_muri.Rdata")
load("ProcessedData/detect_env_prey_muri.Rdata")
detect_data_merge_muri <- detect_data_merge_muri %>% mutate(BestTaxon = as.factor(BestTaxon))
detect_data_all_muri <- detect_data_all_muri %>% 
  mutate(BestTaxon = as.factor(BestTaxon)) %>% 
  filter(BestTaxon != "Berardius bairdii")

### Get bathymetry for pred grid -----------------------------------------------

bathy <- getNOAA.bathy(lon1 = min(detect_data_merge_muri$lon_deg), 
                       lon2 = max(detect_data_merge_muri$lon_deg), 
                       lat1 = min(detect_data_merge_muri$lat_deg),  
                       lat2 = max(detect_data_merge_muri$lat_deg),
                       resolution = 1)

bathy_raster <- marmap::as.raster(bathy)
bathy_r <- rast(bathy_raster)

### Q2.0 Detection rate smoothed over depth, env, and prey presence/absence
### with shape and intercept by species ----------------------------------------

detect_data_prey_presence_muri <- detect_data_all_muri %>% 
  mutate(across(Cololabis:Thunnus, ~ceiling(.)))

m2.0 <-
  bam(Detected ~ 
        s(depth, by = BestTaxon, bs= "ts") + 
        #s(bathy, by = BestTaxon, bs="ts") +
        s(distShore, by = BestTaxon, bs = "ts") +
        s(slope, by = BestTaxon, bs="ts") +
        s(Chla, by = BestTaxon, bs = "ts") +
        s(SST07, by = BestTaxon, bs="ts") +
        s(curVel, by = BestTaxon, bs="ts") +
        s(SSS07, by = BestTaxon, bs="ts") +
        #s(Merluccius, by = BestTaxon, bs="ts") +
        s(Engraulis, by = BestTaxon, bs="ts") +
        s(Trachurus, by = BestTaxon, bs="ts"),
      #s(mld_dr003_7, by = BestTaxon, bs = "ts"),
      family = "binomial",
      method = "fREML",
      data = detect_data_prey_presence_muri,
      discrete = TRUE)

summary(m2.0)

AIC(m2.0)
#954

ggplot(detect_data_all_muri, aes(y = Engraulis, x = BestTaxon, 
                                 fill = as.factor(Detected))) +
  geom_boxplot(outliers = FALSE) +
  #facet_wrap("predator", scales = "free", ncol = 2) +
  theme_minimal() +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))
