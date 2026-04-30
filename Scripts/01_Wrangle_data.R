#### predatorPrey Distribution
#### Wrangle data from MURI: 2019 Hake Survey
#### Spring 2026
#### AVC, from Van Cise et al. MURI zDistribution paper


library(tidyverse)

## Get data --------------------------------------------------------------------

metadata <- read.csv("./Data/Hake_2019_metadata.csv")

detect_data_raw <- read.csv("./Data/M3_compiled_taxon_table_wide.csv") %>% 
  pivot_longer(-c(BestTaxon, Class), names_to = "SampleUID", values_to = "nReads") %>% 
  group_by(SampleUID) %>% 
  mutate(totalReads = sum(nReads)) %>% 
  separate(SampleUID, into = c("Sample_name", NA), remove = FALSE, sep = "_") %>% 
  separate(Sample_name, into = c("run", "primer", "pop", "sample", "dilution", "techRep", "seqRep"), remove = FALSE, sep = "\\.") %>% 
  unite(pop:sample, col = "NWFSCsampleID", sep = "-") %>% 
  mutate(techRep = as.numeric(techRep)) %>% 
  mutate(run = gsub("a","",run)) %>% 
  mutate(run = gsub("b","",run)) %>% 
  mutate(run = gsub("c","",run)) %>% 
  mutate(Detected = ifelse(nReads>0, 1, 0)) %>% 
  filter(primer == "MV1") %>% 
  #filter(techRep < 4) %>% 
  filter(!BestTaxon %in% c("Moschus", "Equus caballus")) %>% 
  ungroup()

#test raw data
detect_data_raw %>% filter(Class == "Mammalia") %>% distinct(BestTaxon)
detect_data_raw %>% distinct(NWFSCsampleID) %>% summarize(n())
detect_data_raw %>% group_by(run) %>% distinct(run)

## Reduce sequencing reps ------------------------------------------------------

detect_data_1seq <- detect_data_raw %>% 
  group_by(run, primer, NWFSCsampleID, dilution, techRep, seqRep) %>% 
  mutate(totReads = sum(nReads)) %>% 
  ungroup() %>% 
  group_by(primer, NWFSCsampleID, dilution, techRep) %>% 
  filter(totReads == max(totReads)) %>% 
  ungroup() %>% 
  mutate(seqRep = replace_na(seqRep, "sr1")) %>% 
  filter(!(totReads == 0 & seqRep %in% c("sr2", "sr3"))) %>% 
  select(-totReads)

detect_data_1seq %>% distinct(NWFSCsampleID) %>% summarize(n())
detect_data_1seq %>% filter(Class == "Mammalia") %>% distinct(BestTaxon)
detect_data_1seq %>% group_by(NWFSCsampleID) %>% 
  summarize(nReps = n()/330) %>% arrange(desc(nReps))
test <- detect_data_1seq %>% filter(NWFSCsampleID == "52193-555")

## Reduce dilutions ------------------------------------------------------------

detect_data_1dil <- detect_data_1seq %>% 
  group_by(run, primer, NWFSCsampleID, dilution, techRep) %>% 
  mutate(totReads = sum(nReads)) %>% 
  ungroup() %>% 
  group_by(primer, NWFSCsampleID, techRep) %>% 
  filter(totReads == max(totReads)) %>% 
  mutate(dilution = as.numeric(substr(dilution, 2, nchar(dilution)))) %>% 
  select(-totReads) %>% 
  filter(dilution == min(dilution)) %>% 
  ungroup()

detect_data_1dil <- detect_data_1seq %>% 
  group_by(NWFSCsampleID, dilution) %>% 
  mutate(nReps = max(techRep)) %>% 
  ungroup() %>% 
  group_by(NWFSCsampleID) %>% 
  filter(nReps == max(nReps)) %>%
  mutate(dilution = as.numeric(substr(dilution, 2, nchar(dilution)))) %>%
  filter(dilution == min(dilution)) %>% 
  select(-nReps) %>% 
  ungroup()
 #filter(nReps == 3)

detect_data_1dil %>% distinct(NWFSCsampleID) %>% summarize(n())
detect_data_1dil %>% filter(Class == "Mammalia") %>% 
  distinct(BestTaxon)
detect_data_1dil %>% group_by(NWFSCsampleID) %>% 
  summarize(nReps = n()/330) %>% arrange(nReps)
detect_data_1dil %>% group_by(NWFSCsampleID) %>% 
  summarize(nReps = n()/330) %>% pull(nReps) %>% max()
#should be 3

## Add metadata ----------------------------------------------------------------

detect_data_meta <- detect_data_1dil %>% 
  left_join(metadata, by = c("NWFSCsampleID" = "sampleID")) 

detect_data_meta %>% group_by(station, depth) %>% 
  distinct(techRep) %>% summarize(nReps = n()) %>% arrange(nReps)

## Remove duplicate samples and stations with fewer than 3 replicates ----------

detect_data <- detect_data_meta %>% 
  group_by(station, depth) %>%
  #should we include bioReps?
  mutate(nReps = max(techRep)) %>% 
  filter(nReps == 3) %>% 
  ungroup()
  

detect_data %>% group_by(station) %>% n_groups() #177 station
detect_data %>% group_by(depth, station) %>% n_groups() #527 station/depths

## Remove Delphinidae family
detect_data <- detect_data_meta %>% 
  filter(!(BestTaxon %in% c('Delphinidae')))

## count number of marine mammal detections by species -------------------------

detect_per_species <- detect_data %>% 
  filter(Class == "Mammalia") %>% 
  group_by(BestTaxon) %>% 
  summarize(nDetect = sum(Detected))

