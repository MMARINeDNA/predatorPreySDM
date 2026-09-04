#### Distribution w environmental variables
#### Summer 2026
#### AVC&MS

#general
library(tidyverse)
library(PNWColors)

#spatial
library(marmap)
library(terra)
library(sf)

#modeling
library(mgcv)
library(MuMIn)

load("./ProcessedData/detect_and_env_ambon.Rdata")
detect_data_merge_ambon <- detect_data_merge_ambon %>% mutate(BestTaxon = as.factor(species)) %>% 
  separate(species, into = c("sp_group",NA), sep = " ", remove = FALSE)

detect_per_species <- detect_data_merge_ambon %>% 
  group_by(species, geo_group) %>% 
  summarize(nDetect = sum(Detected)) %>% 
  st_drop_geometry()

keep_species <- detect_per_species %>%
  filter(nDetect >= 20) %>%
  separate(species, into = c("genus", "species")) %>% 
  mutate(list_name = paste(genus, geo_group, sep = "."))

### Split into species-level datasets ------------------------------------------

detect_data_list <- split(detect_data_merge_ambon, list(detect_data_merge_ambon$sp_group, detect_data_merge_ambon$geo_group))

# filter out species/geo groups with < 20 detections
detect_data_list <- detect_data_list[names(detect_data_list) %in% keep_species$list_name] %>% 
  lapply(st_drop_geometry)

save(detect_data_list, file = "ProcessedData/detect_data_list_ambon.Rdata")

### Get bathymetry for pred grid -----------------------------------------------

bathy <- getNOAA.bathy(lon1 = min(detect_data_merge_ambon$lon), 
                       lon2 = max(detect_data_merge_ambon$lon), 
                       lat1 = min(detect_data_merge_ambon$lat),  
                       lat2 = max(detect_data_merge_ambon$lat),
                       resolution = 1)

bathy_raster <- marmap::as.raster(bathy)
bathy_r <- rast(bathy_raster)

### Q1.0: Depth smoothed over xy with shape and intercept variable by species --

m1.0summary <- data.frame(model = character(),
                          p = numeric(),
                          deviance = numeric(),
                          AIC = numeric(),
                          kcheckp = numeric(),
                          pear.disp = numeric(),
                          nDetect = numeric(),
                          stringsAsFactors = FALSE)
m1.0list <- list()

for (i in 1:length(detect_data_list)){
m1.0_temp <- bam(Detected ~ ti(lat_AAmeters, lon_AAmeters, bs = "tp"),
               family = "binomial",
               method = "fREML",
               data = detect_data_list[[i]],
               discrete = TRUE)

modelSummary <- data.frame(model = names(detect_data_list)[i],
  p = summary(m1.0_temp)$s.table[1, "p-value"],
  deviance = summary(m1.0_temp)$dev.expl,
  AIC = AIC(m1.0_temp),
  kcheckp = k.check(m1.0_temp)[1,"p-value"],
  pear.disp = sum(residuals(m1.0_temp, type = "pearson")^2) / df.residual(m1.0_temp),
  nDetect = sum(detect_data_list[[i]]$Detected),
  stringsAsFactors = FALSE)

m1.0list[[i]] <- m1.0_temp
m1.0summary <- rbind(m1.0summary, modelSummary)
}

### m1.0 predictions ----------------------------------------------------------

m1.0sePreds_list <- list()

for (i in 1:length(detect_data_list)){
m1.0_pred_grid <- expand_grid(lat_AAmeters = seq(min(detect_data_list[[i]]$lat_AAmeters, na.rm = TRUE),
                                            max(detect_data_list[[i]]$lat_AAmeters, na.rm = TRUE),
                                            by = 10000),
                              lon_AAmeters = seq(min(detect_data_list[[i]]$lon_AAmeters, na.rm = TRUE),
                                            max(detect_data_list[[i]]$lon_AAmeters, na.rm = TRUE),
                                            by = 10000))
# response predictions
m1.0_preds <- predict.bam(m1.0list[[i]], m1.0_pred_grid,
                         se.fit = TRUE)

m1.0_sePreds <- data.frame(m1.0_pred_grid,
                           mu   = binomial()$linkinv(m1.0_preds$fit),
                           low  = binomial()$linkinv(m1.0_preds$fit - 1.96 * m1.0_preds$se.fit),
                           high = binomial()$linkinv(m1.0_preds$fit + 1.96 * m1.0_preds$se.fit),
                           low50  = binomial()$linkinv(m1.0_preds$fit - 0.674 * m1.0_preds$se.fit),
                           high50 = binomial()$linkinv(m1.0_preds$fit + 0.674 * m1.0_preds$se.fit))

m1.0sePreds_list[[i]] <- m1.0_sePreds

}

save(m1.0list, m1.0summary, m1.0sePreds_list, file = "ProcessedData/m1.0.Rdata")

### Detection rate smoothed env variables with shape and intercept by species
m1.1_dredge_list <- list()
m1.1_top_models_list <- list()

for (i in 1:length(detect_data_list)){

  print(i)
  print(Sys.time())
  
  tempdat <- detect_data_list[[i]]
  
  # set global model
  global <- gam(Detected ~
                  ti(lat_AAmeters, lon_AAmeters) +
                  s(bathy) +
                  #s(slope) +
                  s(distShore) +
                  #s(BO_O2) +
                  s(BO_meanChl) +
                  s(BO_meanChl_ss) +
                  #s(MS_sst8) +
                  s(MS_sst9) +
                  s(curVel) +
                  #s(MS_sss8) +
                  s(MS_sss9) +
                  s(iceDist) +
                  s(iceCov),
                  #s(iceThick),
                family = binomial(),
                method = "ML",
                data = tempdat,
                na.action = "na.fail")
  
  m1.1 <- dredge(global)

  # Keep the complete dredge object
  m1.1_dredge_list[[i]] <- m1.1
  
  # Keep the delta < 2 models as a table
  m1.1_top_models_list[[i]] <- m1.1 %>%
    filter(delta < 2) %>%
    mutate(species = names(detect_data_list)[i])
}


# model selection

m1.1_allspecies <- bind_rows(m1.1_top_models_list) %>% 
  filter(!(is.na(species)))

selected_model_summaries <- list()
model_selection_summary <- list()

for (i in 1:length(m1.1_top_models_list)) {
  
  tempdat <- detect_data_list[[i]]
  
  top_models <- get.models(m1.1_dredge_list[[i]], subset = delta < 2)
  dev_exp <- sapply(top_models, function(x) summary(x)$dev.expl)
  
  model_summary <- data.frame(model = names(top_models),
                              delta = sapply(top_models, function(x) AICc(x) - min(sapply(top_models, AICc))),
                              m1.1deviance = dev_exp,
                              n_terms = sapply(top_models, function(x) {
                                st <- summary(x)$s.table
                                if (is.null(st)) 0 else nrow(st)
                              }))
  
  max_dev <- max(model_summary$m1.1deviance)
  
  # rules for model selection are:
  # 1. delta AIC < 1.5
  # 2. within 70% of maximum deviance explained
  # 3. minimize number of parameters
  
  final_model <- model_summary %>% 
    filter(delta < 1.5) %>% 
    filter(m1.1deviance >= 0.70 * max_dev) %>% 
    arrange(n_terms, delta) %>%
    slice(1)
  
  selected_model_summaries[[i]] <- data.frame(species = m1.1_top_models_list[[i]]$species[1],
                                                 summary(top_models[[as.character(final_model$model)]])$s.table)
  
  selected_model_summaries[[i]] <- data.frame(
    species = m1.1_top_models_list[[i]]$species[1],
    parameter = attr(
      terms(top_models[[as.character(final_model$model)]]),
      "term.labels"))
  
  model_selection_summary[[i]] <- data.frame(species = m1.1_top_models_list[[i]]$species[1],
                                              model = final_model$model,
                                              delta = final_model$delta,
                                              m1.1deviance = final_model$m1.1deviance,
                                              n_terms = final_model$n_terms)
  
}

m1.1_model_selection <- bind_rows(model_selection_summary)

m1.1_selected_models <- bind_rows(selected_model_summaries) 
 
### combined deviance table

Q1_model_compare <- m1.0summary %>% 
  mutate(m1.0deviance = deviance) %>% 
  dplyr::select(-deviance) %>% 
  left_join(m1.1_model_selection, by = c("model" = "species")) %>% 
  mutate(deltadeviance = m1.1deviance - m1.0deviance)

### Save -----------------------------------------------------------------------

save(m1.1_model_selection, m1.1_selected_models, 
     m1.1_dredge_list, m1.1_top_models_list, 
     Q1_model_compare, file = "ProcessedData/m1.1.Rdata")


### ALL SPECIES MODEL HERE ####################################################             
# m1.0 <-
#   bam(Detected ~ 
#         # main effects of space, depth, taxon
#         # ti(utm.lon, utm.lat,
#         #    d=2,
#         #    k=13,
#         #    bs="tp")+
#         #  ti(depth,
#         #     k=5,
#         #     bs="ts")+
#         ti(BestTaxon,
#            k=6,
#            bs="re")+
#         # interaction between *everything*
#         ti(utm.lon, utm.lat, depth, BestTaxon,
#            d=c(2,1,1),
#            k=c(13, 5, 6),
#            bs=c("tp","ts", "re"))+
#         # space-taxon effect
#         ti(utm.lon, utm.lat, BestTaxon,
#            d=c(2,1),
#            k=c(13,6),
#            bs=c("tp","re"))+
#         # depth-taxon effect
#         ti(depth, BestTaxon,
#            k=c(5,6),
#            bs=c("ts","re")),
#       family = "binomial",
#       method = "fREML",
#       data = detect_data_merge_muri,
#       discrete = TRUE)
# 
# summary(m1.0)
# # Approximate significance of smooth terms:
# #   edf     Ref.df   Chi.sq  p-value    
# # ti(utm.lon,utm.lat)                 -1.366e-13 -8.709e-13    0.000 1.000000    
# # ti(depth)                            1.330e-05  4.000e+00    0.000 0.373754    
# # ti(BestTaxon)                        8.777e-03  5.000e+00    0.015 9.61e-06 ***
# #   ti(BestTaxon,depth,utm.lon,utm.lat)  4.433e+00  2.500e+01  508.743  < 2e-16 ***
# #   ti(BestTaxon,utm.lon,utm.lat)        4.029e+00  6.000e+00 1843.314 0.000802 ***
# #   ti(depth,BestTaxon)                  1.003e+01  2.400e+01 4152.016  < 2e-16 ***
# #   ---
# #   Rank: 405/407
# # R-sq.(adj) =  0.0224   Deviance explained =   11%
# # fREML = 9832.7  Scale est. = 1         n = 9684
# 
# AIC(m1.0)
# # 1704
# 
# #mean squared Pearson residual dispersion parameter
# sum(residuals(m1.0, type = "pearson")^2) / df.residual(m1.0)
# 
# ### m1.0 predictions ----------------------------------------------------------
# 
# m1.0_pred_grid <- expand_grid(depth = seq(from = 0, to = 500, by = 10),
#                               utm.lat = seq(min(detect_data_merge_muri$utm.lat, na.rm = TRUE),
#                                             max(detect_data_merge_muri$utm.lat, na.rm = TRUE),
#                                             by = 5000),
#                               utm.lon = seq(min(detect_data_merge_muri$utm.lon, na.rm = TRUE),
#                                             max(detect_data_merge_muri$utm.lon, na.rm = TRUE),
#                                             by = 5000),
#                               BestTaxon = as.factor(unique(detect_data_merge_muri$BestTaxon)))
# # response predictions
# m1.0preds <- predict.bam(m1.0, m1.0_pred_grid,
#                          se.fit = TRUE)
# 
# m1.0_sePreds <- data.frame(m1.0_pred_grid,
#                            mu   = binomial()$linkinv(m1.0preds$fit),
#                            low  = binomial()$linkinv(m1.0preds$fit - 1.96 * m1.0preds$se.fit),
#                            high = binomial()$linkinv(m1.0preds$fit + 1.96 * m1.0preds$se.fit),
#                            low50  = binomial()$linkinv(m1.0preds$fit - 0.674 * m1.0preds$se.fit),
#                            high50 = binomial()$linkinv(m1.0preds$fit + 0.674 * m1.0preds$se.fit))
# 
# ### Detection rate smoothed over depth and env variables with shape and intercept by species
# m1.1 <-
#   bam(Detected ~ 
#         s(depth, by = BestTaxon, bs= "ts") + 
#         #s(bathy, by = BestTaxon, bs="ts") +
#         s(distShore, by = BestTaxon, bs = "ts") +
#         s(slope, by = BestTaxon, bs="ts") +
#         s(Chla, by = BestTaxon, bs = "ts") +
#         s(SST07, by = BestTaxon, bs="ts") +
#         s(curVel, by = BestTaxon, bs="ts") +
#         s(SSS07, by = BestTaxon, bs="ts"),
#       #s(mld_dr003_7, by = BestTaxon, bs = "ts"),
#       family = "binomial",
#       method = "fREML",
#       data = detect_data_merge_muri,
#       discrete = TRUE)
# 
# summary(m1.1)
# #Deviance explained = 18%
# 
# AIC(m1.1)
# #1569
# 
# ### Save -----------------------------------------------------------------------
# 
# save(m1.0, m1.0preds, m1.0_sePreds,
#      m1.1,
#      file = "./ProcessedData/Q1_models_preds.Rdata")