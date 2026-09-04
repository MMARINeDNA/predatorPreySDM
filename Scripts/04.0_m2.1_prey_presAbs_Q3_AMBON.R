#### Distribution w prey presence/absence m2.1
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

load("ProcessedData/detect_data_list_ambon.Rdata")
load("ProcessedData/m1.1.Rdata")

### Wrangle prey data to Genus? Add pres/absence columns? Select target prey species?


### Detection rate smoothed over best env variables + prey presence/absence ----
detect_data_list_pres <- lapply(detect_data_list, function(x) {
  x %>% 
    mutate(codDetect = case_when(`Boreogadus saida` > 0 ~ 1,
                                 TRUE ~ 0)) %>% 
    mutate(lanceDetect = case_when(`Ammodytes hexapterus` > 0 ~ 1,
                                   TRUE ~ 0))
})

m2.1_dredge_list <- list()
m2.1_top_models_list <- list()

for (i in 1:length(detect_data_list_pres)){
  
  print(i)
  print(Sys.time())
  
  tempdat <- detect_data_list_pres[[i]]
  
  
  # set global model
  env_terms <- m1.1_selected_models %>% 
    filter(species == names(detect_data_list_pres)[[i]]) %>% 
    pull(parameter)
  
  # m1.1terms <- if (all(c("lat_AAmeters", "lon_AAmeters") %in% env_terms)) {
  #   paste(c("ti(lat_AAmeters, lon_AAmeters)", 
  #           paste0("s(", env_terms[!env_terms %in% c("lat_AAmeters", "lon_AAmeters")], ")")), 
  #         collapse = "+")
  #   } else {
  #     paste0("s(", env_terms, ")", collapse = "+")
  #   }
  
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
                                     m1.1terms, 
                                     "+ codDetect + lanceDetect"))
  
  global <- gam(global_formula,
                family = binomial(),
                method = "ML",
                data = tempdat,
                na.action = "na.fail")
  
  m2.1 <- dredge(global)
  
  # Keep the complete dredge object
  m2.1_dredge_list[[i]] <- m2.1
  
  # Keep the delta < 2 models as a table
  m2.1_top_models_list[[i]] <- m2.1 %>%
    filter(delta < 2) %>%
    mutate(species = names(detect_data_list_pres)[i])
}

m2.1_allspecies <- bind_rows(m2.1_top_models_list) %>% 
  filter(!(is.na(species)))

selected_model_summaries <- list()
model_selection_summary <- list()

for (i in 1:length(m2.1_top_models_list)) {
  
  tempdat <- detect_data_list_pres[[i]]
  
  top_models <- get.models(m2.1_dredge_list[[i]], subset = TRUE)
  dev_exp <- sapply(top_models, function(x) summary(x)$dev.expl)
  
 
  model_summary <- data.frame(model = names(top_models),
                              formula = sapply(top_models, function(x)
                                paste(deparse(formula(x)), collapse = " ")),
                              delta = m2.1_dredge_list[[i]][names(top_models), "delta"],
                              m2.1deviance = sapply(top_models, function(x) summary(x)$dev.expl),
                              n_terms = sapply(top_models, function(x) length(attr(terms(x), "term.labels"))))
  
  m1_dev <- m1.1_model_selection %>%
    filter(species == names(detect_data_list_pres)[i]) %>%
    pull(m1.1deviance)
  
  # rules for prey model selection are:
  # 1. delta AIC < 2
  # 2. deviance explained >= best env model
  # 3. minimize number of parameters
  
  final_model <- model_summary %>% 
    filter(delta < 2) %>% 
    filter(m2.1deviance >= m1_dev) %>% 
    arrange(n_terms, delta) %>%
    slice(1)
  
  selected_model_summaries[[i]] <- data.frame(
    species = m2.1_top_models_list[[i]]$species[1],
    parameter = attr(
      terms(top_models[[as.character(final_model$model)]]),
      "term.labels"))
  
  model_selection_summary[[i]] <- data.frame(species = m2.1_top_models_list[[i]]$species[1],
                                             model = final_model$model,
                                             delta = final_model$delta,
                                             m2.1deviance = final_model$m2.1deviance,
                                             n_terms = final_model$n_terms)
  
}

m2.1_model_selection <- bind_rows(model_selection_summary)

m2.1_selected_models <- bind_rows(selected_model_summaries) 

### combined deviance table

Q1Q2_deviance_compare <- Q1_model_compare %>% 
  dplyr::select(model,m1.0deviance, m1.1deviance, n_terms) %>% 
  rename("n_terms_m1.1" = n_terms) %>% 
  left_join(m2.1_model_selection %>% 
              dplyr::select(species, m2.1deviance, n_terms) %>% 
              rename("n_terms_m2.1" = n_terms), 
            by = c("model" = "species")) %>% 
  left_join(m2.1_selected_models %>%
      group_by(species) %>%
      summarise(m1.2parameter = paste(parameter, collapse = " + ")),
    by = c("model" = "species")) %>% 
  left_join(m1.1_selected_models %>%
      group_by(species) %>%
      summarise(parameter_m1.1 = paste(parameter, collapse = " + ")),
    by = c("model" = "species")) 
