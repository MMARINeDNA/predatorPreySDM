#### Q1 AUC and TSS metrics
### Summer 2026
### AVC and MS

#general
library(tidyverse)
library(PNWColors)

#modeling
library(mgcv)
library(pROC)
library(tidysdm)

#load data
#load("ProcessedData/m1.1_ambon.Rdata")
load("ProcessedData/detect_data_list_ambon.Rdata")
#load("ProcessedData/m1.0_ambon.Rdata")

### m1.0 AUC and TSS -----------------------------------------------------------
set.seed(1234)
m1.0_AUC_df <- data.frame()
m1.0_TSS_df <- data.frame()
m1.0_TSS_plot <- list()

for (k in c(1:6)) {
m1.0_AUC_list <- list()
m1.0_TSS_list <- list()
m1.0preds <- list()

detect_data_AUC <- detect_data_list[[k]] %>% 
  mutate(groupNo = sample(1:5, size = nrow(.), replace = TRUE))

for (i in 1:5){
  
  detect_data_train <- detect_data_AUC %>% 
    filter(!(groupNo == i))
  
  detect_data_test <- detect_data_AUC %>% 
    filter(groupNo == i)
  
  m1.0train <- bam(Detected ~ ti(lat_AAmeters, lon_AAmeters, bs = "tp"),
                   family = "binomial",
                   method = "fREML",
                   data = detect_data_train,
                   discrete = TRUE)
  
  m1.0preds[[i]] <- predict(m1.0train, detect_data_test, type = "response")
  
  roc_object <- roc(detect_data_test$Detected, m1.0preds[[i]])
  
  m1.0_AUC_list[[i]] <- auc(roc_object)
 
}
  #TSS
  lower <- max(sapply(m1.0preds, min))
  upper <- min(sapply(m1.0preds, max))
  
  threshold <- seq(lower, upper, by = 0.001)
  
  TSS_df <- list()
  for (i in 1:5) {
    detect_data_test <- detect_data_AUC %>% 
      filter(groupNo == i)
    
  if (length(threshold) > 1) {
  for (j in 1:length(threshold)){
    TSS_df[[length(TSS_df) + 1]] <- (detect_data_test %>% 
      as.data.frame() %>% 
      bind_cols(as.data.frame(m1.0preds[[i]])) %>% 
      mutate(Detected = as.factor(Detected)) %>% 
      mutate(m1.0preds = case_when(m1.0preds[[i]] > threshold[j]~1,
                                   TRUE~0)) %>% 
      mutate(m1.0preds = as.factor(m1.0preds)) %>% 
      tss(., Detected, m1.0preds))$.estimate 
  }} else {
  next
}}
    # test <- detect_data_test %>% 
    #   as.data.frame() %>% 
    #   bind_cols(as.data.frame(m1.0preds)) %>% 
    #   mutate(Detected = as.factor(Detected)) %>% 
    #   mutate(m1.0predD = case_when(m1.0preds > threshold[j]~1,
    #                                TRUE~0)) %>% 
    #   mutate(m1.0predD = as.factor(m1.0predD))
  
  m1.0_TSS_list[[i]] <- as.data.frame(unlist(TSS_df)) %>% 
    rename("TSS" = 1) %>% 
    filter(!(TSS %in% c("tss", "binary"))) %>% 
    mutate("threshold" = rep(threshold, 5)) %>% 
    mutate("test" = i)

AUC_df <- do.call(rbind.data.frame, m1.0_AUC_list) %>% 
  rename("AUC" = 1) %>% 
  mutate(species = names(detect_data_list)[[k]])

m1.0_AUC_df <- rbind(m1.0_AUC_df, AUC_df)

TSS_df <- as.data.frame(do.call(rbind.data.frame, m1.0_TSS_list)) %>% 
  mutate(TSS = as.numeric(TSS)) %>% 
  group_by(threshold) %>% 
  summarize(mean = mean(TSS), sd = sd(TSS), nSamp = n()) %>% 
  mutate(se = sd/(nSamp^(1/2))) %>% 
  mutate(low95 = mean - 1.96 * se,
         high95 = mean + 1.96 * se) %>% 
  mutate(species = names(detect_data_list)[[k]])

m1.0_TSS_df <- rbind(m1.0_TSS_df, TSS_df)

m1.0_TSS_plot[[k]] <- ggplot(TSS_df, aes(x = threshold)) +
  geom_smooth(aes(y = mean), 
              alpha = 0.6, method = "gam") +
  #geom_ribbon(aes(ymin = low95, ymax = high95), alpha = 0.6) +
  #geom_line(aes(y = mean)) +
  theme_minimal() +
  ylab("TSS")

}

rm(m1.0train, m1.0_TSS_list, TSS_df, AUC_df, roc_object, 
   detect_data_train, detect_data_test, detect_data_AUC,
   m1.0_AUC_list, i, j, k, lower, upper, threshold)

save(m1.0_AUC_df, m1.0_TSS_df, m1.0_TSS_plot, file = "ProcessedData/m1.0_performance_ambon.Rdata")
