#### predatorPrey Distribution
#### AMBON sampling map
#### AVC summer 2026

library(tidyverse)
library(ggOceanMaps)
library(sf)

## Get data --------------------------------------------------------------------

load("./ProcessedData/detect_data_ambon.Rdata")


## Get station locations -------------------------------------------------------

station_loc <- detect_data_ambon %>% 
  group_by(station_id) %>% 
  slice_head() %>% 
  ungroup() %>% 
  st_as_sf(coords = c("lon", "lat"),
           crs = 4326,       # WGS84: lon/lat
           remove = FALSE)

## Convex hulls around station groups ------------------------------------------

hulls <- station_loc %>%
  group_by(geo_group) %>%
  summarise(geometry = st_union(geometry),
    .groups = "drop") %>%
  st_concave_hull(ratio = 1,
    allow_holes = FALSE,
    dist = 500)

## Sample map ------------------------------------------------------------------

map <- basemap(limits = c(min(station_loc$lon)-0.2,
                          max(station_loc$lon)+0.2,
                          min(station_loc$lat)-0.1,
                          max(station_loc$lat)+0.1),
               bathy.style = "rcb", crs = 4326,
               rotate = FALSE) +
  geom_sf(data = station_loc, 
                                size = 1.5, alpha = 0.8,
                                aes(color = geo_group)) +
  theme(#axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        #axis.text.y=element_blank(),  
        axis.ticks.y=element_blank(),
        legend.position = "bottom",
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.grid = element_blank()) +
  guides(fill = "none") +
  ggspatial::annotation_scale(location = "bl", width_hint = 0.3)+
  guides(color=guide_legend(override.aes=list(fill=NA))) +
  theme(plot.margin = margin(0, 0, 0, 0))  +
  geom_sf(
    data = hulls,
    
    alpha = 0.3)

map

save(map, file = "./ProcessedData/sampling_map.Rdata")

