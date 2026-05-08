#### 3D Distribution w environmental variables
#### Summer 2026
#### AVC&MS, m1.0 from zDist m3.0c

library(mgcv)
library(tidyverse)
library(PNWColors)
library(marmap)
library(terra)

load("./ProcessedData/detect_data.Rdata")
detect_data <- detect_data %>% mutate(BestTaxon = as.factor(BestTaxon))

### Get bathymetry for pred grid -----------------------------------------------

bathy <- getNOAA.bathy(lon1 = min(detect_data$lon), 
                       lon2 = max(detect_data$lon), 
                       lat1 = min(detect_data$lat),  
                       lat2 = max(detect_data$lat),
                       resolution = 1)

bathy_raster <- as.raster(bathy)
bathy_r <- rast(bathy_raster)

### Q1.0: Depth smoothed over xy with shape and intercept variable by species --

m1.0 <-
  bam(Detected ~ 
        # main effects of space, depth, taxon
        ti(lon, lat,
           d=2,
           k=20,
           bs="tp")+
         ti(depth,
            k=5,
            bs="ts")+
        ti(BestTaxon,
           k=16,
           bs="re")+
        # interaction between *everything*
        ti(lon, lat, depth, BestTaxon,
           d=c(2,1,1),
           k=c(20, 5, 16),
           bs=c("tp","ts", "re"))+
        # space-taxon effect
        ti(lon, lat, BestTaxon,
           d=c(2,1),
           k=c(10,16),
           bs=c("tp","re"))+
        # depth-taxon effect
        ti(depth, BestTaxon,
           k=c(10,16),
           bs=c("ts","re")),
      family = "binomial",
      method = "fREML",
      data = detect_data,
      discrete = TRUE)

summary(m1.0)
# Approximate significance of smooth terms:
#                               edf         Ref.df Chi.sq p-value    
#   ti(lon,lat)                 1.403e+01   16.43  29.86  0.0233 *  
#   ti(depth)                   3.748e-05    4.00   0.00  0.8355    
#   ti(BestTaxon)               1.255e+01   15.00 124.53  <2e-16 ***
#   ti(BestTaxon,depth,lon,lat) 9.049e+01 1216.00 179.41  <2e-16 ***
#   ti(lon,lat,BestTaxon)       4.447e+01  142.00 102.60  <2e-16 ***
#   ti(depth,BestTaxon)         2.834e+01  144.00 100.54  <2e-16 ***
#24.7% deviance explained
#24.6% deviance explained without depth

AIC(m1.0)
# 4985 with all terms
# 4985 with non-significant term (depth) removed

#mean squared Pearson residual dispersion parameter
sum(residuals(m1.0, type = "pearson")^2) / df.residual(m1.0)

### m1.0 predictions ----------------------------------------------------------

m1.0_pred_grid <- expand_grid(depth = seq(from = 0, to = 500, by = 10),
                               lat = seq(min(detect_data$lat, na.rm = TRUE),
                                         max(detect_data$lat, na.rm = TRUE),
                                         by = 0.05),
                               lon = seq(min(detect_data$lon, na.rm = TRUE),
                                         max(detect_data$lon, na.rm = TRUE),
                                         by = 0.05),
                               BestTaxon = as.factor(c("Lagenorhynchus obliquidens",
                                                        "Megaptera novaeangliae",
                                                        "Berardius bairdii")))
# response predictions
m1.0preds <- predict.bam(m1.0, m1.0_pred_grid,
                          se.fit = TRUE)

m1.0_sePreds <- data.frame(m1.0_pred_grid,
                            mu   = binomial()$linkinv(m1.0preds$fit),
                            low  = binomial()$linkinv(m1.0preds$fit - 1.96 * m1.0preds$se.fit),
                            high = binomial()$linkinv(m1.0preds$fit + 1.96 * m1.0preds$se.fit),
                            low50  = binomial()$linkinv(m1.0preds$fit - 0.674 * m1.0preds$se.fit),
                            high50 = binomial()$linkinv(m1.0preds$fit + 0.674 * m1.0preds$se.fit))


### Save -----------------------------------------------------------------------

save(m1.0, m1.0preds, m1.0_sePreds,
     file = "./ProcessedData/Q1_models_preds_0.05degree.Rdata")

