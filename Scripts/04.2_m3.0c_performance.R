### Check m3.0c model performance with AUC
### November 2025

library(mgcv)
library(pROC)
library(tidysdm)

load("./ProcessedData/detect_data.RData")
load("./ProcessedData/detect_data_clean.RData")

#### group data ----------------------------------------------------------------

random_samples <- detect_data %>% 
  mutate(BestTaxon = as.factor(BestTaxon)) %>% 
  ungroup() %>% 
  select(NWFSCsampleID) %>% 
  distinct(NWFSCsampleID) %>% 
  mutate(groupNo = sample(1:10, size = nrow(.), replace = TRUE))


detect_data_AUC <- detect_data %>% 
  filter(BestTaxon %in% c(detect_per_species %>% 
                            filter(nDetect > 10) %>% 
                            pull(BestTaxon))) %>% 
  mutate(BestTaxon = as.factor(BestTaxon)) %>% 
  left_join(random_samples, by = "NWFSCsampleID")

#### AUC and TSS test ----------------------------------------------------------

AUC_list <- list()
TSS_list <- list()

for (i in 1:10){
  
  detect_data_train <- detect_data_AUC %>% 
    filter(!(groupNo == i))
  
  detect_data_test <- detect_data_AUC %>% 
    filter(groupNo == i)
  
  m3.0ctrain <- m3.0c <-
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
        data = detect_data_train,
        discrete = TRUE)
  
  m3.0cpreds <- predict(m3.0ctrain, detect_data_test, type = "response")
  
  roc_object <- roc(detect_data_test$Detected, m3.0cpreds)
  
  AUC_list[[i]] <- auc(roc_object)
  
  threshold <- seq(0.005,max(m3.0cpreds),by=0.001)
  
  TSS_df <- list()
  
  for (j in 1:length(threshold)){
    TSS_df[[j]] <- detect_data_test %>% 
      as.data.frame() %>% 
      bind_cols(as.data.frame(m3.0cpreds)) %>% 
      mutate(Detected = as.factor(Detected)) %>% 
      mutate(m3.0cpreds = case_when(m3.0cpreds > threshold[j]~1,
                                   TRUE~0)) %>% 
      mutate(m3.0cpreds = as.factor(m3.0cpreds)) %>% 
      tss(., Detected, m3.0cpreds)
    
  }
  
  TSS_list[[i]] <- as.data.frame(unlist(TSS_df)) %>% 
    rename("TSS" = 1) %>% 
    filter(!(TSS %in% c("tss", "binary"))) %>% 
    mutate("threshold" = threshold) %>% 
    mutate("test" = i)
}

AUC_df_m3.0c <-do.call(rbind.data.frame, AUC_list) %>% 
  rename("AUC" = 1)

TSS_df_m3.0c <- as.data.frame(do.call(rbind.data.frame, TSS_list)) %>% 
  mutate(TSS = as.numeric(TSS)) %>% 
  group_by(threshold) %>% 
  summarize(mean = mean(TSS), sd = sd(TSS), nSamp = n()) %>% 
  mutate(se = sd/(nSamp^(1/2))) %>% 
  mutate(low95 = mean - 1.96 * se,
         high95 = mean + 1.96 * se)

TSS_plot_m3.0c <- ggplot(TSS_df_m3.0c, aes(x = threshold)) +
  geom_smooth(aes(ymin = low95, ymax = high95, y = mean), 
              stat = "identity") +
  theme_classic() +
  ylab("TSS")

### Do it again with "clean" dataset -------------------------------------------

#### group data ----------------------------------------------------------------

random_samples_clean <- detect_data_clean %>% 
  mutate(BestTaxon = as.factor(BestTaxon)) %>% 
  ungroup() %>% 
  select(NWFSCsampleID) %>% 
  distinct(NWFSCsampleID) %>% 
  mutate(groupNo = sample(1:10, size = nrow(.), replace = TRUE))


detect_data_AUC_clean <- detect_data_clean %>% 
  filter(BestTaxon %in% c(detect_per_species %>% 
                            filter(nDetect > 10) %>% 
                            pull(BestTaxon))) %>% 
  mutate(BestTaxon = as.factor(BestTaxon)) %>% 
  left_join(random_samples_clean, by = "NWFSCsampleID")

#### AUC and TSS test ----------------------------------------------------------

AUC_list_clean <- list()
TSS_list_clean <- list()

for (i in 1:10){
  
  detect_data_train <- detect_data_AUC_clean %>% 
    filter(!(groupNo == i))
  
  detect_data_test <- detect_data_AUC_clean %>% 
    filter(groupNo == i)
  
  m3.0ctrain <- m3.0c <-
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
        data = detect_data_train,
        discrete = TRUE)
  
  m3.0cpreds <- predict(m3.0ctrain, detect_data_test, type = "response")
  
  roc_object <- roc(detect_data_test$Detected, m3.0cpreds)
  
  AUC_list_clean[[i]] <- auc(roc_object)
  
  threshold <- seq(0.005,max(m3.0cpreds),by=0.001)
  
  TSS_df <- list()
  
  for (j in 1:length(threshold)){
    TSS_df[[j]] <- detect_data_test %>% 
      as.data.frame() %>% 
      bind_cols(as.data.frame(m3.0cpreds)) %>% 
      mutate(Detected = as.factor(Detected)) %>% 
      mutate(m3.0cpreds = case_when(m3.0cpreds > threshold[j]~1,
                                    TRUE~0)) %>% 
      mutate(m3.0cpreds = as.factor(m3.0cpreds)) %>% 
      tss(., Detected, m3.0cpreds)
    
  }
  
  TSS_list_clean[[i]] <- as.data.frame(unlist(TSS_df)) %>% 
    rename("TSS" = 1) %>% 
    filter(!(TSS %in% c("tss", "binary"))) %>% 
    mutate("threshold" = threshold) %>% 
    mutate("test" = i)
}

AUC_df_m3.0c_clean <- do.call(rbind.data.frame, AUC_list_clean) %>% 
  rename("AUC" = 1)

TSS_df_m3.0c_clean <- as.data.frame(do.call(rbind.data.frame, TSS_list_clean)) %>% 
  mutate(TSS = as.numeric(TSS)) %>% 
  group_by(threshold) %>% 
  summarize(mean = mean(TSS), sd = sd(TSS), nSamp = n()) %>% 
  mutate(se = sd/(nSamp^(1/2))) %>% 
  mutate(low95 = mean - 1.96 * se,
         high95 = mean + 1.96 * se)

TSS_plot_m3.0c_clean <- ggplot(TSS_df_m3.0c_clean, aes(x = threshold)) +
  geom_smooth(aes(ymin = low95, ymax = high95, y = mean), 
              stat = "identity") +
  theme_classic() +
  ylab("TSS")

save(AUC_df_m3.0c, TSS_df_m3.0c, TSS_plot_m3.0c, 
     AUC_df_m3.0c_clean, TSS_df_m3.0c_clean, TSS_plot_m3.0c_clean,
     file = "./ProcessedData/m3.0cperformance.Rdata")

