#### Q1 m1.1 AUC and TSS metrics
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
load("ProcessedData/m1.1.Rdata")
load("ProcessedData/detect_data_list_ambon.Rdata")

### m1.1 AUC and TSS -----------------------------------------------------------
set.seed(1234)
m1.1_AUC_df <- data.frame()
m1.1_TSS_df <- data.frame()
m1.1_TSS_plot <- list()

for (k in c(1:6)) {
  m1.1_AUC_list <- list()
  m1.1_TSS_list <- list()
  m1.1preds <- list()
  
  detect_data_AUC <- detect_data_list[[k]] %>% 
    mutate(groupNo = sample(1:5, size = nrow(.), replace = TRUE))
  
  env_terms <- m1.1_selected_models %>% 
    filter(species == names(detect_data_list)[[k]]) %>% 
    pull(parameter)
  
  m1.1terms <- if (all(c("lat_AAmeters", "lon_AAmeters") %in% env_terms)) {
    
    other_terms <- env_terms[
      !env_terms %in% c("lat_AAmeters", "lon_AAmeters")
    ]
    
    paste(
      c(
        "ti(lat_AAmeters, lon_AAmeters)",
        if (length(other_terms) > 0) paste0("s(", other_terms, ")")
      ),
      collapse = "+"
    )
    
  } else {
    
    paste0("s(", env_terms, ")", collapse = "+")
  }
  
  global_formula <- as.formula(paste("Detected ~", 
                                     m1.1terms))
  
  for (i in 1:5){
    
    detect_data_train <- detect_data_AUC %>% 
      filter(!(groupNo == i))
    
    detect_data_test <- detect_data_AUC %>% 
      filter(groupNo == i)
    
    m1.1train <- bam(global_formula,
                     family = "binomial",
                     method = "fREML",
                     data = detect_data_train,
                     discrete = TRUE)
    
    m1.1preds[[i]] <- predict(m1.1train, detect_data_test, type = "response")
    
    roc_object <- roc(detect_data_test$Detected, m1.1preds[[i]])
    
    m1.1_AUC_list[[i]] <- auc(roc_object)
    
  }
  #TSS
  lower <- max(sapply(m1.1preds, min))
  upper <- min(sapply(m1.1preds, max))
  
  threshold <- seq(lower, upper, by = 0.001)
  
  TSS_df <- list()
  for (i in 1:5) {
    detect_data_test <- detect_data_AUC %>% 
      filter(groupNo == i)
    
    if (length(threshold) > 1) {
      for (j in 1:length(threshold)){
        TSS_df[[length(TSS_df) + 1]] <- (detect_data_test %>% 
                                           as.data.frame() %>% 
                                           bind_cols(as.data.frame(m1.1preds[[i]])) %>% 
                                           mutate(Detected = as.factor(Detected)) %>% 
                                           mutate(m1.1preds = case_when(m1.1preds[[i]] > threshold[j]~1,
                                                                        TRUE~0)) %>% 
                                           mutate(m1.1preds = as.factor(m1.1preds)) %>% 
                                           tss(., Detected, m1.1preds))$.estimate 
      }} else {
        next
      }}

  m1.1_TSS_list[[i]] <- as.data.frame(unlist(TSS_df)) %>% 
    rename("TSS" = 1) %>% 
    filter(!(TSS %in% c("tss", "binary"))) %>% 
    mutate("threshold" = rep(threshold, 5)) %>% 
    mutate("test" = i)
  
  AUC_df <- do.call(rbind.data.frame, m1.1_AUC_list) %>% 
    rename("AUC" = 1) %>% 
    mutate(species = names(detect_data_list)[[k]])
  
  m1.1_AUC_df <- rbind(m1.1_AUC_df, AUC_df)
  
  TSS_df <- as.data.frame(do.call(rbind.data.frame, m1.1_TSS_list)) %>% 
    mutate(TSS = as.numeric(TSS)) %>% 
    group_by(threshold) %>% 
    summarize(mean = mean(TSS), sd = sd(TSS), nSamp = n()) %>% 
    mutate(se = sd/(nSamp^(1/2))) %>% 
    mutate(low95 = mean - 1.96 * se,
           high95 = mean + 1.96 * se) %>% 
    mutate(species = names(detect_data_list)[[k]])
  
  m1.1_TSS_df <- rbind(m1.1_TSS_df, TSS_df)
  
  m1.1_TSS_plot[[k]] <- ggplot(TSS_df, aes(x = threshold)) +
    geom_smooth(aes(y = mean), 
                alpha = 0.6, method = "gam") +
    #geom_ribbon(aes(ymin = low95, ymax = high95), alpha = 0.6) +
    #geom_line(aes(y = mean)) +
    theme_minimal() +
    ylab("TSS")
  
}

rm(m1.1train, m1.1_TSS_list, TSS_df, AUC_df, roc_object, 
   detect_data_train, detect_data_test, detect_data_AUC,
   m1.1_AUC_list, i, j, k, lower, upper, threshold)

save(m1.1_AUC_df, m1.1_TSS_df, m1.1_TSS_plot, file = "ProcessedData/m1.1_performance.Rdata")
