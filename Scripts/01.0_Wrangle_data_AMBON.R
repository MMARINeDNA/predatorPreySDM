#### Wrangle AMBON marine mammal detection data
#### AVC 2025

#### set up environment --------------------------------------------------

#general
library(tidyverse)
library(ggridges)
library(PNWColors)
#maps
library(ggOceanMaps)
library(scatterpie)
library(mapdata)

taxa <- read.csv("./data/AFSC_MV1/taxonomy_20250609_collapsed.csv", row.names = 1) %>% 
  dplyr::select(kingdom,phylum,class,family,genus,species) %>% 
  rownames_to_column("ASV")

metadata <- read.csv("./data/metadata/FAIRe_noaa-afsc-dbo1.csv", row.names = 1) %>% 
  rownames_to_column("sample") %>% 
  dplyr::select(sample, samp_category, materialSampleID, station_id, 
         collection_depth, eventDate, decimalLatitude, decimalLongitude)

detect_data <- read.csv("./data/AFSC_MV1/ASVtable.csv") %>% 
  separate(1, into = c("sample", "rep")) %>% 
  group_by(sample) %>% 
  summarize(across(c(ASV1:ASV880), ~ceiling(mean(.x, na.rm = TRUE)))) %>% 
  ungroup() %>% 
  mutate(across(everything(), ~replace_na(., 0))) %>% 
  pivot_longer(-sample, names_to = "ASV", values_to = "Detected") %>% 
  mutate(Detected = case_when(Detected > 0~1,
                              TRUE~0)) %>% 
  left_join(taxa, by = "ASV") %>% 
  filter(class == "Mammalia") %>% 
  filter(!(family %in% c("Hominidae", "Bovidae", "Suidae", "Canidae"))) %>% 
  left_join(metadata, "sample") %>% 
  filter(samp_category == "sample") %>% 
  filter(!grepl("not applicable", station_id)) %>% 
  group_by(station_id, species) %>% 
  summarize(Detected = ifelse(sum(Detected) > 0,1,0))
  
  
save(detect_data, metadata, file = "./data products/detect_data.Rdata")
  
### OLD ######################################################################
#### format data ----------------------------------------------------------
#PMEL data
pmel_detect_data <- PMEL_ASVs %>% 
  ungroup() %>% 
  filter(Class == "Mammalia") %>% 
  filter(Cruise_ID_short %in% c("SKQ2021", "NO20", "SKQ23-12S")) %>% 
  filter(!(Family %in% c("Hominidae", "Bovidae", "Suidae", "Camelidae",
                         "Castoridae", "Cricetidae", "Felidae",
                         "Leporidae"))) %>%
  filter(!(Genus %in% c("Canis", "Lagenorhynchus",
                        "Orcaella", "Enhydra",
                        "Delphinus", "Stenella", "NA"))) %>% 
  filter(!(is.na(Species))) %>% 
  mutate(Species = case_when(Species == "Balaena mysticetus" ~ "bowhead whale",
                             Species == "Balaenoptera acutorostrata" ~ "minke whale",
                             Species == "Balaenoptera physalus" ~ "fin whale",
                             Species == "Delphinapterus leucas" ~ "beluga whale",
                             Species == "Erignathus barbatus" ~ "bearded seal",
                             Species == "Eschrichtius robustus" ~ "grey whale",
                             Species == "Lagenorhynchus obliquidens" ~ "Pacific white-sided dolphin",
                             Species == "Megaptera novaeangliae" ~ "humpback whale",
                             Species == "Odobenus rosmarus" ~ "walrus",
                             Species == "Phoca fasciata" ~ "ribbon seal",
                             Species == "Phocoena phocoena" ~ "harbor porpoise",
                             Species == "Pusa hispida" ~ "ringed seal",
                             Species == "Balaenoptera musculus" ~ "blue whale",
                             Species == "Orcinus orca" ~ "killer whale",
                             Species == "Phocoenoides dalli" ~ "Dall's porpoise",
                             Species == "Phoca fasciata" ~ "ribbon seal",
                             Species == "Phoca vitulina" ~ "harbor seal",
                             Species == "Eumetopias jubatus" ~ "Steller sea lion",
                             Species == "Phoca largha" ~ "spotted seal",
                             TRUE ~ Species)) %>% 
  group_by(Sample_Name, Species) %>% 
  arrange(desc(nReads), .by_group = TRUE) %>%
  slice_head() %>% 
  dplyr::select(Species,nReads,Sample_Name:lon) %>% 
  ungroup() %>% 
  rename("latitude" = "lat", "longitude" = "lon", 
         "PMEL_ID" = "Sample_Name",
         "depth" = "Depth_m") %>% 
  mutate(collection_year = year(as.Date(Time))) %>% 
  mutate(collection_month = month(as.Date(Time))) %>% 
  mutate(collection_day = day(as.Date(Time))) %>% 
  dplyr::select(-Negative_control, -Cruise_ID_short, -Cruise_ID_long,
         -Cast_No., -Sample_volume_ml, -Time, -Rosette_position) %>% 
  mutate(Technical_Replicate = case_when(collection_year == 2021 ~ 2,
                                         collection_year == 2023 ~ 4,
                                         TRUE ~ Technical_Replicate)) %>% 
  separate(PMEL_ID, into = c("PMEL_ID", NA), sep = "_") %>% 
  mutate(Lab = "PMEL") %>% 
  mutate(PMEL_ID = gsub("SKQ23", "SKQ2023", PMEL_ID))

#AFSC data
AFSC_data_long <- detect_data %>% 
  rename("bowhead whale" = "Balaena.mysticetus",
         "minke whale" = "Balaenoptera.acutorostrata",
         "fin whale" = "Balaenoptera.physalus",
         "beluga whale" = "Delphinapterus.leucas",
         "bearded seal" = "Erignathus.barbatus",
         "grey whale" = "Eschrichtius.robustus",
         "Pacific white-sided dolphin" = "Lagenorhynchus.obliquidens",
         "humpback whale" = "Megaptera.novaeangliae",
         "walrus" = "Odobenus.rosmarus",
         "harbor/spotted seal" = "Phoca",
         "ribbon seal" = "Phoca.fasciata",
         "harbor porpoise" = "Phocoena.phocoena",
         "ringed seal" = "Pusa.hispida") %>% 
  pivot_longer("bowhead whale":"ringed seal", 
               names_to = "Species", 
               values_to = "nReads") %>% 
  relocate(Species:nReads, .after = "ABL_ID") %>% 
  filter(sample_type == "sample") %>% 
  dplyr::select(ABL_ID:depth) %>% 
  rename("Technical_Replicate" = "replicate",
         "Station" = "location1") %>% 
  dplyr::select(-location2, -location3, -sample_type, -source) %>% 
  mutate(Technical_Replicate = case_when(Technical_Replicate == "A" ~ 1,
                                         Technical_Replicate == "B" ~ 2,
                                         Technical_Replicate == "C" ~ 3,
                                         collection_year == 2021 ~ 1)) %>% 
  mutate(PMEL_ID = toupper(PMEL_ID)) %>% 
  separate(PMEL_ID, into = c("PMEL_ID", NA), sep = "\\.") %>% 
  mutate(PMEL_ID = paste(PMEL_ID, paste0("SKQ", collection_year), sep = ".")) %>% 
  mutate(Lab = "AFSC")
                                         

##Combine data
detect_data_long <- bind_rows(AFSC_data_long, pmel_detect_data)

detections_by_station <- detect_data_long %>% 
  mutate(detect = case_when(nReads > 1 ~ 1,
                            TRUE ~ 0)) %>% 
  relocate(detect, .after = nReads) %>% 
  group_by(Station, depth, Species, collection_year) %>% 
  mutate(nDetect = sum(detect), nObs = n()) %>% 
  relocate(nDetect:nObs, .after = detect)

replicate_detections <- detect_data_long %>% 
  mutate(detect = case_when(nReads > 0 ~ 1,
                            TRUE ~ 0)) %>% 
  relocate(detect, .after = nReads) %>% 
  group_by(Station, depth, collection_year, Species) %>% 
  mutate(nDetect = sum(detect), nObs = n()) %>% 
  relocate(nDetect:nObs, .after = detect) %>% 
  slice_head()

replicate_detections %>% group_by(Station, depth, collection_year) %>% 
  filter(nObs >= 3) %>% slice_head()

biorep_detections <- detect_data_long %>% 
  mutate(detect = case_when(nReads > 0 ~ 1,
                            TRUE ~ 0)) %>% 
  relocate(detect, .after = nReads) %>% 
  group_by(location1, collection_year, Species) %>% 
  mutate(nDetect = sum(detect), nObs = n()) %>% 
  relocate(nDetect:nObs, .after = detect) %>% 
  slice_head() %>% 
  filter(nDetect > 0)

#detection years----------------------------------------------------------------

species_by_yr <- detections_by_station %>% 
  group_by(collection_year, Species) %>% 
  summarise(nDetects = sum(nDetect))

### Density ridgeplot ----------------------------------------------------------

positive_detections <- detections_by_station %>% 
  filter(detect == 1) 

rare_detections <- positive_detections %>% 
  group_by(Species) %>% 
  filter(n() < 3)
  
depth_detection <- ggplot(positive_detections, aes(y = Species, x = depth, 
                                                   fill = Species,
                                                   color = Species)) +
  geom_density_ridges(scale = 0.5, bandwidth = 15,
                      jittered_points = TRUE,
                      point_alpha = 1,
                      point_shape = 21,
                      alpha = 0.6) +
  theme_minimal() + 
  coord_flip() +
  scale_x_reverse() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_fill_manual(values = pnw_palette("Cascades",15, type = "continuous")) +
  scale_color_manual(values = pnw_palette("Cascades",15, type = "continuous")) +
  geom_point(data = rare_detections) +
  theme(legend.position = "none") 

#### Detections by depth model -------------------------------------------------

depth.gam <- gam(detect ~ s(depth),
                 data = detections_by_station,
                 family = "binomial",
                 method = "REML")

summary(depth.gam)

spDepth.gam <- gam(detect ~ s(depth, by = factor(Species)),
                 data = detections_by_station,
                 family = "binomial",
                 method = "REML")

summary(spDepth)

#### Map of sampling locations -------------------------------------------------

detect_stations <- detect_data_long %>% 
  drop_na(longitude) %>% 
  group_by(Station, collection_year) %>% 
  slice_head()

number_of_stations <- detect_data_long %>% 
  drop_na(longitude) %>% 
  group_by(Station, depth, collection_year) %>% 
  slice_head() %>% 
  ungroup() %>% 
  group_by(collection_year) %>% 
  summarize(n())

sampling_map <- basemap(limits=c(-175,-150,60,75), bathymetry = TRUE, rotate = TRUE) +
  ggspatial::geom_spatial_point(data = detect_stations, 
                                aes(x = longitude, y = latitude, 
                                    color = as.factor(collection_year),
                                    shape = as.factor(collection_year)),
                                alpha = 0.4,
                                stroke = 1,
                                size = 3) +
  scale_color_manual(values = c("#0f85a0","#4a9152",
                                "#33454e", "#015b58"))+
  ggspatial::stat_spatial_identity(position = "dodge") +
    theme(legend.position = "none")

### Detections by species -----------------------------------------------

positive_detect_species <- detections_by_station %>% 
  filter(detect == 1) %>% 
  group_by(PMEL_ID, Species) %>% 
  slice_head()

species_detections <- positive_detect_species %>% 
  ungroup() %>% 
  group_by(Species) %>% 
  summarize(n())

positive_detect_species_map <- positive_detect_species %>% 
  group_by(Station, Species) %>% 
  mutate(nDetect = n()) %>% 
  dplyr::select(Station,Species, nDetect, latitude,longitude) %>% 
  group_by(Station) %>% 
  # mutate(totalDetect = sum(nDetect)) %>% 
  # mutate(propDetect = nDetect/totalDetect) %>% 
  # select(-nDetect, -totalDetect) %>% 
  group_by(Station, Species) %>% 
  slice_head() %>% 
  pivot_wider(names_from = Species, values_from = nDetect, values_fill = 0) %>% 
  mutate(total = rowSums(across(3:14))) %>% 
  mutate(total = total/25) %>% 
  mutate(total = case_when(total < 0.1~0.1,
                           TRUE~total)) %>% 
  mutate(longitude = longitude + 360)

ak<-map_data('world2Hires','USA:Alaska')

ggplot()+geom_polygon(data=ak,aes(long,lat,group=group),fill=8,color="black") +
  coord_map(xlim=c(-175, -155), ylim=c(60,72)) +
  geom_scatterpie(data = positive_detect_species_map,
                  aes(x=longitude, y=latitude, group=Station),
                  cols=unique(positive_detections$Species),
                              legend_name = "species",
                              sorted_by_radius = TRUE) + 

  scale_fill_manual(values = c(pnw_palette("Cascades",10, type = "continuous"),
                               pnw_palette("Shuksan",10, type = "continuous"))) +
  theme(legend.position = "none") +
  theme_minimal()


save(detect_data_long, file = "data products/AMBON_total_detections_MM.csv")

save(depth_detection, sampling_map, file = "data products/AMBON_sampling_detections.Rdata")

