#### Distribution w environmental variables
#### Summer 2026
#### AVC&MS, m1.0 from zDist m3.0c

#general
library(tidyverse)
library(PNWColors)

#mapping
library(marmap)
library(terra)

#modeling
library(mgcv)

load("./ProcessedData/detect_and_env_ambon.Rdata")
detect_data_merge_ambon <- detect_data_merge_ambon %>% mutate(BestTaxon = as.factor(species)) %>% 
  separate(species, into = c("sp_group",NA), sep = " ", remove = FALSE)

detect_per_species <- detect_data_merge_ambon %>% 
  group_by(species, geo_group) %>% 
  summarize(nDetect = sum(Detected)) %>% 
  st_drop_geometry()

### Split into species-level datasets ------------------------------------------

detect_data_list <- split(detect_data_merge_ambon, list(detect_data_merge_ambon$sp_group, detect_data_merge_ambon$geo_group))

list2env(setNames(detect_data_list,
    paste0("detect_data_", names(detect_data_list))), envir = .GlobalEnv)

rm(detect_data_list)

save(list = ls(pattern = "detect_data_"), file = "ProcessedData/detect_data_subset_ambon.RData")

### Get bathymetry for pred grid -----------------------------------------------

bathy <- getNOAA.bathy(lon1 = min(detect_data_merge_ambon$lon), 
                       lon2 = max(detect_data_merge_ambon$lon), 
                       lat1 = min(detect_data_merge_ambon$lat),  
                       lat2 = max(detect_data_merge_ambon$lat),
                       resolution = 1)

bathy_raster <- marmap::as.raster(bathy)
bathy_r <- rast(bathy_raster)

### Q1.0: Depth smoothed over xy with shape and intercept variable by species --

## Odobenus
m1.0_OrA <- bam(Detected ~ ti(lat, lon, bs = "tp", k = 13),
               family = "binomial",
               method = "fREML",
               data = detect_data_Odobenus.Arctic,
               discrete = TRUE)

summary(m1.0_OrA)
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -2.3317     0.3319  -7.026 2.13e-12 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df Chi.sq p-value  
# ti(lat,lon) 8.404   10.5  18.48  0.0502 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.116   Deviance explained = 21.1%
# fREML = 65.954  Scale est. = 1         n = 237

AIC(m1.0_OrA)
# 141

#mean squared Pearson residual dispersion parameter
sum(residuals(m1.0_OrA, type = "pearson")^2) / df.residual(m1.0_OrA)
#0.65

### m1.0 predictions ----------------------------------------------------------

## Odobenus
m1.0_OrA_pred_grid <- expand_grid(lat = seq(min(detect_data_Odobenus.Arctic$lat, na.rm = TRUE),
                                            max(detect_data_Odobenus.Arctic$lat, na.rm = TRUE),
                                            by = 1),
                              lon = seq(min(detect_data_Odobenus.Arctic$lon, na.rm = TRUE),
                                            max(detect_data_Odobenus.Arctic$lon, na.rm = TRUE),
                                            by = 1))
# response predictions
m1.0_OrApreds <- predict.bam(m1.0_OrA, m1.0_OrA_pred_grid,
                         se.fit = TRUE)

m1.0_OrA_sePreds <- data.frame(m1.0_OrA_pred_grid,
                           mu   = binomial()$linkinv(m1.0_OrApreds$fit),
                           low  = binomial()$linkinv(m1.0_OrApreds$fit - 1.96 * m1.0_OrApreds$se.fit),
                           high = binomial()$linkinv(m1.0_OrApreds$fit + 1.96 * m1.0_OrApreds$se.fit),
                           low50  = binomial()$linkinv(m1.0_OrApreds$fit - 0.674 * m1.0_OrApreds$se.fit),
                           high50 = binomial()$linkinv(m1.0_OrApreds$fit + 0.674 * m1.0_OrApreds$se.fit))

### Detection rate smoothed over depth and env variables with shape and intercept by species
m1.1_OrA <-
  bam(Detected ~ 
        ti(lat, lon, bs = "tp", k = 13) +
        s(bathy, bs="ts") +
        #s(distShore, bs = "ts") +
        s(BO_meanChl, bs = "ts") +
        #s(BO_meanChl_ss, bs = "ts") +
        s(MS_sst8, bs="ts") +
        #s(MS_sst9, bs="ts") +
        #s(curVel, bs="ts") +
        #s(MS_sss8, bs="ts") +
        #s(MS_sss9, bs="ts") +
        #s(iceDist, bs = "ts") +
        s(iceCov, bs = "ts"),
        #s(iceThick, bs = "ts"),
      family = "binomial",
      method = "fREML",
      data = detect_data_Odobenus.Arctic,
      discrete = TRUE)

summary(m1.1_OrA)
#Deviance explained = 21.3%

AIC(m1.1_OrA)
#128

### Save -----------------------------------------------------------------------

save(m1.0_OrA, m1.0_OrA_sePreds, m1.1_OrA, file = "ProcessedData/Q1_OrA.Rdata")

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