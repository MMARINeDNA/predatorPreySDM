#### 3D Distribution
#### Figures for Q3 models
#### October 2025
#### AVC

library(tidyverse)
library(PNWColors)
library(mgcv)
library(patchwork)
library(terra)
library(ggOceanMaps)
library(sf)
library(cowplot)
library(ggspatial)
library(marmap)

#load("./ProcessedData/m3.0models_preds_0.05degree.Rdata")
load("ProcessedData/m3.0models_preds_0.05degree_depthmask.Rdata")
load("./ProcessedData/detect_data.Rdata")
load("./ProcessedData/detect_data_clean.RData")
mmEcoEvo <- read.csv("./Data/MM_metadata.csv")
metadata <- read.csv("./Data/Hake_2019_metadata.csv")

det_colors <- c("500" = "#03051AFF", "300" = "#611F53FF", "150" = "#CB1B4FFF","50" = "#F16445FF","0" = "#F69C73FF")

#"#751F58FF" "#AC185AFF" "#DB2946FF" 
### m3.0c max POD depth map for three species ----------------------------------

# pull depth of max POD
maxPOD_depth <- m3.0c_sePreds %>% 
  group_by(BestTaxon, lat,lon) %>% #3816 groups
  arrange(desc(mu), .by_group = TRUE) %>% 
  mutate(max_mu = first(mu), max_low50 = first(low50), 
         max_high50 = first(high50), max_depth = first(depth)) %>%
  filter(mu > max_low50 & mu < max_high50) %>%
  mutate(depth_min = min(depth), depth_max = max(depth)) %>% 
  slice_head() %>% 
  ungroup() %>% 
  mutate(depthWidth = depth_max - depth_min) %>% 
  mutate(ci95 = high-low) %>% 
  ungroup()
  
# create convex hull study area
study_area <- st_as_sf(metadata, coords = c("lon", "lat"), crs = 4326) %>%
  summarise(geometry = st_union(geometry)) %>%
  st_convex_hull() 

# convert POD to sf
maxPOD_depth_sf <- maxPOD_depth %>% 
  mutate(lon_plain = lon, lat_plain = lat) %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326) 

maxPOD_depth_clipped <- st_join(study_area, maxPOD_depth_sf, left = TRUE) %>% 
  mutate(depth = case_when(max_mu < 0.005~NA,
                           TRUE~depth)) %>% 
  mutate(depthWidth = case_when(max_mu < 0.005~NA,
                           TRUE~depthWidth))

# pull depth of detections
pos_detect <- detect_data %>% 
  filter(BestTaxon %in% c("Lagenorhynchus obliquidens",
                          "Megaptera novaeangliae",
                          "Berardius bairdii")) %>% 
  filter(Detected == 1) %>% 
  mutate(depth = case_when(depth %in% c(48,50)~50,
                           depth %in% c(467,485,495,500)~500,
                           TRUE~depth))
  

### Create coastline shapefile -------------------------------------------------
world <- st_read("Data/ne_10m_land/ne_10m_land.shp")

data_bbox <- st_as_sf(maxPOD_depth, coords = c("lon", "lat"), crs = 4326) %>%
  st_bbox() %>%
  st_as_sfc()  %>%
  st_buffer(dist = 70000)

westcoast_land <- st_crop(world, data_bbox)

### Get bathymetry data --------------------------------------------------------

# lon‐range and lat‐range:
lon1 <- min(maxPOD_depth_clipped$lon_plain); lon2 <- max(maxPOD_depth_clipped$lon_plain) 
lat1 <- min(maxPOD_depth_clipped$lat_plain); lat2 <- max(maxPOD_depth_clipped$lat_plain)

# Download bathymetry
bath <- getNOAA.bathy(lon1 = lon1, lon2 = lon2,
                      lat1 = lat1, lat2 = lat2,
                      resolution = 1)  # “1” ~ 1-minute (~1.8 km) resolution

# Convert bathy to a data.frame for ggplot
bath_df <- fortify.bathy(bath) 

save(maxPOD_depth, pos_detect, maxPOD_depth_clipped, file = "./ProcessedData/H3.0c_pred_flt.Rdata")

#### Max POD plot with matching color scale ------------------------------------

depth_max_detect <- ggplot(westcoast_land) +
  geom_tile(data = maxPOD_depth_clipped, 
            aes(x = lon_plain, y = lat_plain, fill = depth)) +
  geom_sf(fill = "grey50", colour = NA) +
  scale_fill_viridis_c(name = "Depth (m)",
                       option = "mako",
                       trans = "reverse",
                       begin = 0.4, end = 0.9, na.value = "transparent") +
  ggspatial::geom_spatial_point(data = pos_detect,
                                aes(x = lon, y = lat,
                                    color = as.factor(depth)),
                                size = 1,
                                alpha = 0.8,
                                stroke = 1,
                                position = position_jitter(width = 0.05,
                                                           height = 0.05)) +
  scale_color_manual(values = det_colors) +
  facet_wrap(~BestTaxon, labeller = label_wrap_gen(width=10)) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.margin = margin(0, 0, 0, 0),
        legend.position = "right", #c(0.54, 0.45),    # <<-- Adjust to place legend over land
        legend.justification = c("left"),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        strip.text = element_text(face = "italic")) +
  geom_contour(data = bath_df, aes(x = x, y = y, z = z),
               breaks = c(-500, -1000, -2000), color = "grey70", linewidth = 0.3) +
  guides(color = guide_legend("Detection\ndepth (m)"))

depth_max_detect

ci95_plot <- ggplot(westcoast_land) +
  geom_tile(data = maxPOD_depth_clipped, 
            aes(x = lon_plain, y = lat_plain, 
                fill = depthWidth)) +
  theme_minimal() +
  #ggtitle(names(species)[i]) +
  geom_sf(fill = "grey50", colour = NA) +
  scale_fill_viridis_c(name = "50% CI\ndepth range (m)",
                       option = "magma",
                       trans = "reverse",
                       begin = 0.15, end = 1, na.value = "transparent") +
  facet_wrap(~BestTaxon, labeller = label_wrap_gen(width=10)) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.margin = margin(0, 0, 0, 0),
        legend.position = "right", #c(0.54, 0.45),    # <<-- Adjust to place legend over land
        legend.justification = c("left"),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        strip.text = element_text(face = "italic")) +
  geom_contour(data = bath_df, aes(x = x, y = y, z = z),
               breaks = c(-500, -1000, -2000), color = "grey70", linewidth = 0.3)

depth_max_detect / ci95_plot

save(depth_max_detect, ci95_plot, file = "./Figures/H3.0c_map.Rdata")

### Repeat with "clean" dataset ------------------------------------------------

### m3.0c max POD depth map for three species ----------------------------------

# pull depth of max POD
maxPOD_depth_clean <- m3.0c_clean_sePreds %>% 
  group_by(BestTaxon, lat,lon) %>% #3816 groups
  arrange(desc(mu), .by_group = TRUE) %>% 
  mutate(max_mu = first(mu), max_low50 = first(low50), 
         max_high50 = first(high50), max_depth = first(depth)) %>%
  filter(mu > max_low50 & mu < max_high50) %>%
  mutate(depth_min = min(depth), depth_max = max(depth)) %>% 
  slice_head() %>% 
  ungroup() %>% 
  mutate(depthWidth = depth_max - depth_min) %>% 
  mutate(ci95 = high-low) %>% 
  ungroup()

# create convex hull study area
study_area_clean <- st_as_sf(metadata, coords = c("lon", "lat"), crs = 4326) %>%
  summarise(geometry = st_union(geometry)) %>%
  st_convex_hull() 

# convert POD to sf
maxPOD_depth_sf_clean <- maxPOD_depth_clean %>% 
  mutate(lon_plain = lon, lat_plain = lat) %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326) 

maxPOD_depth_clipped_clean <- st_join(study_area_clean, maxPOD_depth_sf_clean, left = TRUE) %>% 
  mutate(depth = case_when(max_mu < 0.016~NA,
                           TRUE~depth)) %>% 
  mutate(depthWidth = case_when(max_mu < 0.016~NA,
                                TRUE~depthWidth))

# pull depth of detections
pos_detect_clean <- detect_data_clean %>% 
  filter(BestTaxon %in% c("Lagenorhynchus obliquidens",
                          "Megaptera novaeangliae",
                          "Berardius bairdii")) %>% 
  filter(Detected == 1) %>% 
  mutate(depth = case_when(depth %in% c(48,50)~50,
                           depth %in% c(467,485,495,500)~500,
                           TRUE~depth))



### Get bathymetry data --------------------------------------------------------

# lon‐range and lat‐range:
lon1 <- min(maxPOD_depth_clipped_clean$lon_plain); lon2 <- max(maxPOD_depth_clipped_clean$lon_plain) 
lat1 <- min(maxPOD_depth_clipped_clean$lat_plain); lat2 <- max(maxPOD_depth_clipped_clean$lat_plain)

# Download bathymetry
bath_clean <- getNOAA.bathy(lon1 = lon1, lon2 = lon2,
                      lat1 = lat1, lat2 = lat2,
                      resolution = 1)  # “1” ~ 1-minute (~1.8 km) resolution

# Convert bathy to a data.frame for ggplot
bath_df_clean <- fortify.bathy(bath_clean) 

#save(maxPOD_depth_clean, pos_detect_clean, maxPOD_depth_clipped_clean, file = "./ProcessedData/H3.0c_pred_flt_clean.Rdata")

#### Max POD plot with matching color scale ------------------------------------
sp_labels <- c("Berardius bairdii" = "northern giant\nbottlenose whale", "Lagenorhynchus obliquidens" = "Pacific\nwhite-sided dolphin", 
               "Megaptera novaeangliae" = "humpback whale")

depth_max_detect_clean <- ggplot(westcoast_land) +
  geom_tile(data = maxPOD_depth_clipped_clean, 
            aes(x = lon_plain, y = lat_plain, fill = depth)) +
  geom_sf(fill = "grey50", colour = NA) +
  scale_fill_viridis_c(name = "Maximized detection\ndepth (m)",
                       option = "mako",
                       trans = "reverse",
                       begin = 0.4, end = 0.9, na.value = "grey90") +
  ggspatial::geom_spatial_point(data = pos_detect_clean,
                                aes(x = lon, y = lat,
                                    color = as.factor(depth)),
                                size = 1,
                                alpha = 0.8,
                                stroke = 1,
                                position = position_jitter(width = 0.05,
                                                           height = 0.05)) +
  scale_color_manual(values = det_colors) +
  facet_wrap(~BestTaxon, labeller = as_labeller(sp_labels), ncol = 1) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.margin = margin(0, 0, 0, 0),
        legend.position = "bottom", #c(0.54, 0.45),    # <<-- Adjust to place legend over land
        legend.justification = c("center"),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9)) +
  geom_contour(data = bath_df_clean, aes(x = x, y = y, z = z),
               breaks = c(-500, -1000, -2000), color = "grey50", linewidth = 0.3) +
  guides(color = guide_legend("Detection\ndepth (m)"))

depth_max_detect_clean

ci95_plot_clean <- ggplot(westcoast_land) +
  geom_tile(data = maxPOD_depth_clipped_clean, 
            aes(x = lon_plain, y = lat_plain, 
                fill = depthWidth)) +
  theme_minimal() +
  #ggtitle(names(species)[i]) +
  geom_sf(fill = "grey50", colour = NA) +
  scale_fill_viridis_c(name = "Model uncertainty:\n50% CI strip width,\ndepth (m)",
                       option = "magma",
                       trans = "reverse",
                       begin = 0.15, end = 1, na.value = "grey90") +
  facet_wrap(~BestTaxon, labeller = as_labeller(sp_labels), ncol = 1) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.margin = margin(0, 0, 0, 0),
        legend.position = "bottom", #c(0.54, 0.45),    # <<-- Adjust to place legend over land
        legend.justification = c("center"),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9)) +
  geom_contour(data = bath_df_clean, aes(x = x, y = y, z = z),
               breaks = c(-500, -1000, -2000), color = "grey50", linewidth = 0.3)

(depth_max_detect_clean | ci95_plot_clean) + 
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

#### Plot all as one grid ####

depth_max_detect_long <- st_drop_geometry(maxPOD_depth_clipped_clean) %>% 
  #select(-geometry) %>% 
  select(BestTaxon, lon_plain, lat_plain, depth, depthWidth) %>% 
  pivot_longer(depth:depthWidth, names_to = "pred_type", values_to = "pred_value")

depth_max_detect_grid <- ggplot(westcoast_land) +
  geom_tile(data = subset(depth_max_detect_long, pred_type == "depth"),
            aes(x = lon_plain, y = lat_plain, fill = pred_value)) +
  scale_fill_viridis_c(name = "Maximized detection\ndepth (m)",
                       option = "mako",
                       trans = "reverse",
                       begin = 0.4, end = 0.9, na.value = "grey90") +
  ggnewscale::new_scale_fill() +
  geom_tile(data = subset(depth_max_detect_long, pred_type == "depthWidth"),
            aes(x = lon_plain, y = lat_plain, fill = pred_value)) +
  scale_fill_viridis_c(name = "Model uncertainty:\n50% CI strip width,\ndepth (m)",
                       option = "magma",
                       trans = "reverse",
                       begin = 0.15, end = 1,
                       na.value = "grey90") +
  geom_sf(fill = "grey50", colour = NA) +
  ggspatial::geom_spatial_point(data = pos_detect_clean,
                                aes(x = lon, y = lat,
                                    color = as.factor(depth)),
                                size = 1,
                                alpha = 0.8,
                                stroke = 1,
                                position = position_jitter(width = 0.05,
                                                           height = 0.05)) +
  scale_color_manual(values = det_colors) +
  facet_grid(BestTaxon~pred_type,
             labeller = labeller(BestTaxon = as_labeller(sp_labels))) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.margin = margin(0, 0, 0, 0),
        legend.position = "right", #c(0.54, 0.45),    # <<-- Adjust to place legend over land
        legend.justification = c("center"),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        strip.text.x = element_blank()) +
  geom_contour(data = bath_df_clean, aes(x = x, y = y, z = z),
               breaks = c(-500, -1000, -2000), color = "grey50", linewidth = 0.3) +
  guides(color = guide_legend("Detection\ndepth (m)"))

depth_max_detect_grid

save(depth_max_detect_clean, ci95_plot_clean, depth_max_detect_grid, file = "./Figures/H3.0c_map_clean.Rdata")

################################################################################
### Max POD with individual color scales ---------------------------------------

### depth of max POD map -------------------------------------------------------
# depth_max_detect <- list()
# 
# species <- c("Lagenorhynchus obliquidens",
#              "Megaptera novaeangliae",
#              "Berardius bairdii")
# names(species) <- c("Lobl","Mnov","Bbai")
# 
# for (i in 1:length(species)){
# depth_max_detect[[i]] <- ggplot(westcoast_land) +
#   geom_tile(data = maxPOD_depth_clipped %>% 
#               filter(BestTaxon == species[i]), 
#             aes(x = lon_plain, y = lat_plain, fill = depth)) +
#   geom_sf(fill = "grey50", colour = NA) +
#   scale_fill_viridis_c(name = "Depth(m)",
#                        option = "mako",
#                        trans = "reverse",
#                        begin = 0.4, end = 0.9, na.value = "transparent") +
#   ggspatial::geom_spatial_point(data = pos_detect %>%
#                                   filter(BestTaxon == species[i]),
#                                 aes(x = lon, y = lat,
#                                     color = as.factor(depth)),
#                                 size = 1,
#                                 alpha = 0.8,
#                                 stroke = 1,
#                                 position = position_jitter(width = 0.05,
#                                                            height = 0.05)) +
#   scale_color_manual(values = det_colors) +
#   #scale_color_viridis_d("Detection depth (m)", option = "rocket",
#   #                      direction = -1, begin = 0.3, end = 0.8, guide = "none") +
#   ggtitle(names(species)[i]) +
#   theme_classic() +
#   theme(axis.text = element_blank(),
#         axis.ticks = element_blank(),
#         axis.title = element_blank(),
#         plot.margin = margin(0, 0, 0, 0),
#         legend.position = c(0.54, 0.45),    # <<-- Adjust to place legend over land
#         legend.justification = c("left"),
#         legend.background = element_blank(),
#         legend.box.background = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         legend.title = element_text(size = 10),
#         legend.text = element_text(size = 9)) +
#   geom_contour(data = bath_df, aes(x = x, y = y, z = z),
#                breaks = c(-500, -1000, -2000), color = "grey70", linewidth = 0.3) +
#   guides(color = "none")
# }
#   
# #depth_max_detect[[1]] + depth_max_detect[[2]] + depth_max_detect[[3]] 
# 
# # "legend-only" plot 
# dummy_df <- data.frame(
#   x = 5, y = 1,
#   depth = factor(unique(pos_detect$depth), 
#                  levels = sort(unique(pos_detect$depth))))
# 
# dummy_plot <- ggplot(dummy_df, aes(x = x, y = y, color = depth)) +
#   geom_point() +
#   scale_color_manual(values = det_colors) +
#   # scale_color_viridis_d(
#   #   name = "Detection depth (m)",
#   #   option = "rocket",
#   #   direction = -1,
#   #   begin = 0.3, end = 0.8) +
#   theme_minimal() +
#   theme(legend.position = "bottom",
#     plot.margin = margin(0, 0, 0, 0),
#     legend.title = element_text(size = 11),
#     legend.text = element_text(size = 10))
# 
# # Extract legend 
# g <- ggplotGrob(dummy_plot)
# leg_index <- which(sapply(g$grobs, function(x) x$name) == "guide-box")
# shared_legend <- g$grobs[[leg_index]]
# 
# # Combine maps and legend
# 
# combined_plot <- cowplot::plot_grid(
#   plot_grid(plotlist = depth_max_detect, nrow = 1),
#   shared_legend,
#   ncol = 1,
#   rel_heights = c(1, 0.08))
# 
# combined_plot
# 
# # 95% CI of max POD
# 
# ci95_plot_list <- list()
# 
# for (i in 1:length(species)){
# ci95_plot_list[[i]] <- ggplot(westcoast_land) +
#   geom_tile(data = maxPOD_depth_clipped %>% 
#               filter(BestTaxon == species[i]), 
#             aes(x = lon_plain, y = lat_plain, 
#                                           fill = depthWidth)) +
#   theme_classic() +
#   #ggtitle(names(species)[i]) +
#   geom_sf(fill = "grey50", colour = NA) +
#   scale_fill_viridis_c(name = "50% CI\ndepth range",
#                        option = "magma",
#                        trans = "reverse",
#                        begin = 0.15, end = 1, na.value = "transparent") +
#   theme(axis.text = element_blank(),
#         axis.ticks = element_blank(),
#         axis.title = element_blank(),
#         plot.margin = margin(0, 0, 0, 0),
#         legend.position = c(0.54, 0.45),    # <<-- Adjust to place legend over land
#         legend.justification = c("left"),
#         legend.background = element_blank(),
#         legend.box.background = element_blank(),
#         legend.title = element_text(size = 10),
#         legend.text = element_text(size = 9)) +
#   geom_contour(data = bath_df, aes(x = x, y = y, z = z),
#                breaks = c(-500, -1000, -2000), color = "grey70", linewidth = 0.3)
# 
# if (i == 1){
#   ci95_plot_list[[i]] <- ci95_plot_list[[i]] +
#     annotation_north_arrow(location = "bl",
#     which_north = "true",
#     style = north_arrow_fancy_orienteering(fill = c("grey30", "white"),
#                                            line_col = "grey30"),
#     height = unit(1, "cm"),
#     width = unit(1, "cm"))
# }
# }
# 
# #ci95_plot_list[[1]] + ci95_plot_list[[2]] + ci95_plot_list[[3]]
# 
# depthPODmax <- plot_grid(
#   plot_grid(plotlist = depth_max_detect, nrow = 1),
#   shared_legend,
#   plot_grid(plotlist = ci95_plot_list, nrow = 1),
#   ncol = 1,
#   rel_heights = c(1, 0.08, 1))
# 
# pdf(height = 15, file = "./Figures/maxPODdepthmap.pdf")
# depthPODmax
# dev.off()

#### Do it again with poop separate --------------------------------------------
# load("./ProcessedData/m3.0c_baleen.Rdata")
# 
# # pull depth of max POD
# maxPOD_depth_baleen <- m3.0c_baleen_sePreds %>% 
#   group_by(BestTaxon, lat,lon) %>% #3816 groups
#   arrange(desc(mu), .by_group = TRUE) %>% 
#   mutate(max_mu = first(mu), max_low50 = first(low50), 
#          max_high50 = first(high50), max_depth = first(depth)) %>%
#   filter(mu > max_low50 & mu < max_high50) %>%
#   mutate(depth_min = min(depth), depth_max = max(depth)) %>% 
#   slice_head() %>% 
#   ungroup() %>% 
#   mutate(depthWidth = depth_max - depth_min) %>% 
#   mutate(ci95 = high-low) %>% 
#   ungroup()
# 
# # convert POD to sf
# maxPOD_depth_sf_baleen <- maxPOD_depth_baleen %>% 
#   mutate(lon_plain = lon, lat_plain = lat) %>% 
#   st_as_sf(coords = c("lon", "lat"), crs = 4326) 
# 
# maxPOD_depth_clipped_baleen <- st_join(study_area, maxPOD_depth_sf_baleen, left = TRUE) %>% 
#   mutate(depth = case_when(max_mu < 0.005~NA,
#                            TRUE~depth)) %>% 
#   mutate(depthWidth = case_when(max_mu < 0.005~NA,
#                                 TRUE~depthWidth))
# 
# # pull depth of detections
# pos_detect_baleen <- detect_data_baleen %>% 
#   filter(BestTaxon %in% c("Megaptera novaeangliae",
#                           "Megaptera novaeangliae_poop")) %>% 
#   filter(Detected == 1) %>% 
#   mutate(depth = case_when(depth %in% c(48,50)~50,
#                            depth %in% c(467,485,495,500)~500,
#                            TRUE~depth))
# 
# ### depth of max POD map -------------------------------------------------------
# depth_max_detect_baleen <- list()
# 
# species <- c("Megaptera novaeangliae",
#              "Megaptera novaeangliae_poop")
# names(species) <- c("Mnov","Mnov_poop")
# 
# for (i in 1:length(species)){
#   depth_max_detect_baleen[[i]] <- ggplot(westcoast_land) +
#     geom_tile(data = maxPOD_depth_clipped_baleen %>% 
#                 filter(BestTaxon == species[i]), 
#               aes(x = lon_plain, y = lat_plain, fill = depth)) +
#     geom_sf(fill = "grey50", colour = NA) +
#     scale_fill_viridis_c(name = "Depth(m)",
#                          option = "mako",
#                          trans = "reverse",
#                          begin = 0.4, end = 0.9, na.value = "transparent") +
#     ggspatial::geom_spatial_point(data = pos_detect_baleen %>%
#                                     filter(BestTaxon == species[i]),
#                                   aes(x = lon, y = lat,
#                                       color = as.factor(depth)),
#                                   size = 1,
#                                   alpha = 0.8,
#                                   stroke = 1,
#                                   position = position_jitter(width = 0.05,
#                                                              height = 0.05)) +
#     scale_color_manual(values = det_colors) +
#     #scale_color_viridis_d("Detection depth (m)", option = "rocket",
#     #                      direction = -1, begin = 0.3, end = 0.8, guide = "none") +
#     ggtitle(names(species)[i]) +
#     theme_classic() +
#     theme(axis.text = element_blank(),
#           axis.ticks = element_blank(),
#           axis.title = element_blank(),
#           plot.margin = margin(0, 0, 0, 0),
#           legend.position = c(0.54, 0.45),    # <<-- Adjust to place legend over land
#           legend.justification = c("left"),
#           legend.background = element_blank(),
#           legend.box.background = element_blank(),
#           panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(),
#           legend.title = element_text(size = 10),
#           legend.text = element_text(size = 9)) +
#     geom_contour(data = bath_df, aes(x = x, y = y, z = z),
#                  breaks = c(-500, -1000, -2000), color = "grey70", linewidth = 0.3) +
#     guides(color = "none")
# }
# 
# #depth_max_detect[[1]] + depth_max_detect[[2]] + depth_max_detect[[3]] 
# 
# # "legend-only" plot 
# dummy_df_baleen <- data.frame(
#   x = 5, y = 1,
#   depth = factor(unique(pos_detect_baleen$depth), 
#                  levels = sort(unique(pos_detect_baleen$depth))))
# 
# dummy_plot_baleen <- ggplot(dummy_df_baleen, aes(x = x, y = y, color = depth)) +
#   geom_point() +
#   scale_color_viridis_d(
#     name = "Detection depth (m)",
#     option = "rocket",
#     direction = -1,
#     begin = 0.3, end = 0.8) +
#   theme_minimal() +
#   theme(legend.position = "bottom",
#         plot.margin = margin(0, 0, 0, 0),
#         legend.title = element_text(size = 11),
#         legend.text = element_text(size = 10))
# 
# # Extract legend 
# g_baleen <- ggplotGrob(dummy_plot_baleen)
# leg_index_baleen <- which(sapply(g_baleen$grobs, function(x) x$name) == "guide-box")
# shared_legend_baleen <- g_baleen$grobs[[leg_index_baleen]]
# 
# # Combine maps and legend
# 
# combined_plot_baleen <- cowplot::plot_grid(
#   plot_grid(plotlist = depth_max_detect_baleen, nrow = 1),
#   shared_legend_baleen,
#   ncol = 1,
#   rel_heights = c(1, 0.08))
# 
# combined_plot_baleen
# 
# # 95% CI of max POD
# 
# ci95_plot_list_baleen <- list()
# 
# for (i in 1:length(species)){
#   ci95_plot_list_baleen[[i]] <- ggplot(westcoast_land) +
#     geom_tile(data = maxPOD_depth_clipped_baleen %>% 
#                 filter(BestTaxon == species[i]), 
#               aes(x = lon_plain, y = lat_plain, 
#                   fill = depthWidth)) +
#     theme_classic() +
#     #ggtitle(names(species)[i]) +
#     geom_sf(fill = "grey50", colour = NA) +
#     scale_fill_viridis_c(name = "50% CI\ndepth range",
#                          option = "magma",
#                          trans = "reverse",
#                          begin = 0.15, end = 1, na.value = "transparent") +
#     theme(axis.text = element_blank(),
#           axis.ticks = element_blank(),
#           axis.title = element_blank(),
#           plot.margin = margin(0, 0, 0, 0),
#           legend.position = c(0.54, 0.45),    # <<-- Adjust to place legend over land
#           legend.justification = c("left"),
#           legend.background = element_blank(),
#           legend.box.background = element_blank(),
#           legend.title = element_text(size = 10),
#           legend.text = element_text(size = 9)) +
#     geom_contour(data = bath_df, aes(x = x, y = y, z = z),
#                  breaks = c(-500, -1000, -2000), color = "grey70", linewidth = 0.3)
#   
#   if (i == 1){
#     ci95_plot_list_baleen[[i]] <- ci95_plot_list_baleen[[i]] +
#       annotation_north_arrow(location = "bl",
#                              which_north = "true",
#                              style = north_arrow_fancy_orienteering(fill = c("grey30", "white"),
#                                                                     line_col = "grey30"),
#                              height = unit(1, "cm"),
#                              width = unit(1, "cm"))
#   }
# }
# 
# #ci95_plot_list[[1]] + ci95_plot_list[[2]] + ci95_plot_list[[3]]
# 
# depthPODmax_baleen <- plot_grid(
#   plot_grid(plotlist = depth_max_detect_baleen, nrow = 1),
#   shared_legend_baleen,
#   plot_grid(plotlist = ci95_plot_list_baleen, nrow = 1),
#   ncol = 1,
#   rel_heights = c(1, 0.08, 1))
# 
# pdf(height = 15, file = "./Figures/maxPODdepthmap_baleen.pdf")
# depthPODmax_baleen
# dev.off()
# 
# save(depthPODmax_baleen, file = "./Figures/masPODdepthmap_baleen.Rdata")