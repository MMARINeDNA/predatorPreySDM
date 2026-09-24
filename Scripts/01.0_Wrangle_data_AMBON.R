#### predatorPrey Distribution
#### Wrangle AMBON marine mammal detection data
#### AVC 2026

#### set up environment --------------------------------------------------

#general
library(tidyverse)
library(ggridges)
library(PNWColors)
#maps
library(ggOceanMaps)
library(scatterpie)
library(mapdata)
library(sf)

#Note: As of July 2026, this code gloms all ASVs to species level (or genus if 
#there is no species ID) and keeps marine mammal detections and prey species
#read counts separate among technical replicates

### Get MV1 data ---------------------------------------------------------------
#Data downloaded from Kim Ledger's GitHub repo: dbo-marver-mb
#https://github.com/kjledger-NOAA/dbo_marver_mb/tree/main/outputs

taxa <- read.csv("./data/AMBON/taxonomy_20250609_collapsed.csv", row.names = 1) %>% 
  dplyr::select(kingdom,phylum,class,family,genus,species) %>% 
  rownames_to_column("ASV")

metadata <- read.csv("./data/AMBON/FAIRe_noaa-afsc-dbo1.csv", row.names = 1) %>% 
  rownames_to_column("sample") %>% 
  dplyr::select(sample, samp_category, materialSampleID, station_id, 
         collection_depth, eventDate, decimalLatitude, decimalLongitude)

fish_data <- read.csv("./data/AMBON/ASVtable.csv") %>% 
  separate(1, into = c("sample", "rep")) %>% 
  mutate(across(everything(), ~replace_na(., 0))) %>% 
  pivot_longer(-c(sample, rep), names_to = "ASV", values_to = "nReads") %>%
  left_join(taxa, by = "ASV") %>% 
  #filter(class != "Mammalia") %>% 
  mutate(species = case_when(is.na(genus)~ASV,
                            is.na(species)~paste(genus, " spp."),
                             TRUE~species)) %>% 
  dplyr::select(sample, rep, nReads, species) %>% 
  group_by(sample, rep, species) %>% 
  mutate(totReads = sum(nReads)) %>% 
  slice_head() %>% 
  dplyr::select(-nReads) %>% 
  filter(!grepl("ASV", species)) %>% 
  pivot_wider(names_from = "species", values_from = "totReads")
  
# test <- fish_data %>% left_join(metadata, "sample") %>% 
#   filter(samp_category == "sample") %>% 
#   filter(!grepl("not applicable", station_id)) %>%
#   filter(totReads > 0) %>% 
#   group_by(sample) %>% 
#   summarize(sampleReads = sum(totReads), sampleSpecies = n())

# meantest <- test %>% ungroup() %>% 
#   summarize(mean(sampleReads), mean(sampleSpecies))

detect_data_ambon <- read.csv("./data/AMBON/ASVtable.csv") %>% 
  separate(1, into = c("sample", "rep")) %>% 
  mutate(across(everything(), ~replace_na(., 0))) %>% 
  pivot_longer(-c(sample,rep), names_to = "ASV", values_to = "nReads") %>% 
  left_join(taxa, by = "ASV") %>% 
  filter(class == "Mammalia") %>% 
  filter(!(family %in% c("Hominidae", "Bovidae", "Suidae", "Canidae"))) %>% 
  left_join(metadata, "sample") %>% 
  filter(samp_category == "sample") %>% 
  filter(!grepl("not applicable", station_id)) %>% 
  group_by(sample,rep,species) %>% 
  mutate(totReads = sum(nReads)) %>% 
  mutate(Detected = case_when(totReads > 0~1,
                              TRUE~0)) %>% 
  slice_head() %>% 
  ungroup() %>% 
  dplyr::select(-nReads) %>% 
  #The next lines group marine mammal detections by station, which I'm not sure if we should do.
  #group_by(station_id, species) %>% 
  #mutate(Detected = ifelse(sum(Detected) > 0,1,0)) %>% 
  #ungroup()
  left_join(fish_data, by = c("sample","rep")) %>% 
  mutate(lat = as.numeric(decimalLatitude), lon = as.numeric(decimalLongitude)) %>% 
  mutate(geo_group = case_when(lat > 69~"Arctic",
                               TRUE~"Bering")) %>% 
  st_as_sf(coords = c("lon", "lat"),
    crs = 4326,
    remove = FALSE) %>% 
  st_transform(3338) %>% 
  mutate(lon_AAmeters = st_coordinates(.)[,1],
    lat_AAmeters = st_coordinates(.)[,2]) %>% 
  st_drop_geometry()

save(detect_data_ambon, metadata, file = "./ProcessedData/detect_data_ambon.Rdata")

### Summary of species in dataset ----------------------------------------------
# detection frequency (N samples/site) for each species 
# to see what is common in the dataset, to inform pred-prey interactions

detect_freq_ambon <- detect_data_ambon %>%
  filter(samp_category == "sample") %>%
  pivot_longer(`Acipenser baerii`:`Xiphister  spp.`, names_to = "species_detected", values_to = "reads") %>%
  mutate(genus_detected = word(species_detected, 1),          
         site = paste(decimalLatitude, decimalLongitude, collection_depth, sep = "_"),
         samp_rep = paste(sample, rep, sep = "_")) %>%
  group_by(genus_detected) %>%
  summarise(n_samples = n_distinct(samp_rep[reads > 0]),
            n_sites   = n_distinct(site[reads > 0]),
            total_samples = n_distinct(samp_rep),
            total_sites = n_distinct(site),
            pct_samples = round(100 * n_samples / total_samples, 1),
            pct_sites   = round(100 * n_sites   / total_sites,   1),
            .groups = "drop") %>%
  arrange(desc(n_sites)) %>%
  print(n = Inf)

detect_freq_ambon %>%
  filter(n_samples > 0) %>%
  ggplot(aes(pct_sites)) +
  geom_histogram(binwidth = 5) +
  geom_vline(xintercept = c(5, 10, 20), linetype = "dashed")

save(detect_freq_ambon, file = "./ProcessedData/detect_freq_ambon.Rdata")
  