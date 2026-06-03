#### MURI fish ~ marine mammal predator-prey analysis
#### Adapted from 03_0_predatorPrey_correlation_Q2_AMBON.R

library(tidyverse)
library(phyloseq)
library(PNWColors)
library(microViz)
library(ANCOMBC)
library(sf)

### Get data -------------------------------------------------------------------

load("ProcessedData/detect_data_muri.Rdata")

# Identify species (fish/other) columns — everything after the metadata block
species_cols <- names(detect_data_muri %>% dplyr::select(Agonidae:Zoarcidae) %>% 
                        st_drop_geometry())

# Build seqtab: one row per SampleUID (species cols are identical across the 21
# BestTaxon rows for each SampleUID, so take the first), then collapse
# dilution/techRep replicates within NWFSCsampleID via ceiling(mean()),
# matching the AMBON rep-collapse approach.

seqtab <- detect_data_muri %>%
  st_drop_geometry() %>% 
  distinct(SampleUID, NWFSCsampleID, .keep_all = TRUE) %>%
  dplyr::select(NWFSCsampleID, all_of(species_cols)) %>%
  group_by(NWFSCsampleID) %>%
  summarize(across(all_of(species_cols), ~ ceiling(mean(.x, na.rm = TRUE)))) %>%
  mutate(across(everything(), ~ replace_na(., 0))) %>%
  column_to_rownames("NWFSCsampleID") %>% 
  ungroup() 

# Build samdf: one row per NWFSCsampleID (take first rep for metadata)
samdf <- detect_data_muri %>%
  distinct(NWFSCsampleID, .keep_all = TRUE) %>%
  dplyr::select(NWFSCsampleID, sample, station, Niskin, depth, transect, utm.lat, utm.lon,
         lat_deg, lon_deg, water.depth, bathy.bottom.depth, bottom.depth.consensus,
         year, month, day, date, volume, Fluor, Zymo,
         control, drop.sample, field.negative.type, totalReads, nReps) %>%
  column_to_rownames("NWFSCsampleID")

# Build taxon table: BestTaxon is already at species level in this dataset.
# Construct a minimal tax_table from detect_per_species names (all Mammalia).
# For fish we don't have a taxonomy file — use species names as both genus and species.

taxa_fish <- data.frame(ASV = species_cols, label = species_cols) %>% 
  mutate(n_words = str_count(label, "\\S+"),
         rank = case_when(str_detect(label, "idae$") ~ "family",
                          n_words >= 2 ~ "species",
                          TRUE ~ "genus")) %>% 
  mutate(genus = case_when(rank == "species" ~ word(label, 1), 
                           rank == "genus" ~ label,
                           TRUE ~ NA_character_),
         species = case_when(rank == "species" ~ label,
                             TRUE ~ NA_character_),
         family = case_when(rank == "family" ~ label,
                            TRUE ~ NA_character_),
         order = NA_character_,
         class = NA_character_,
         phylum  = NA_character_,
         kingdom = NA_character_) %>% 
  select(ASV, kingdom, phylum, class, order, family, genus, species) %>% 
  column_to_rownames("ASV") %>% 
  as.matrix()

### Make phyloseq object -------------------------------------------------------

ps <- phyloseq(otu_table(seqtab, taxa_are_rows = FALSE),
               sample_data(samdf),
               tax_table(taxa_fish))

### Add mm detects to ps object ------------------------------------------------

# Collapse marine mammal detections: one row per NWFSCsampleID per BestTaxon,
# then take ceiling(mean(Detected)) across dilution/techRep reps
mm_detect_long <- detect_data_muri %>%
  dplyr::select(NWFSCsampleID, BestTaxon, Detected) %>%
  group_by(NWFSCsampleID, BestTaxon) %>%
  summarize(Detected = ceiling(mean(Detected, na.rm = TRUE)), .groups = "drop")

mm_detect <- mm_detect_long %>%
  st_drop_geometry() %>% 
  pivot_wider(names_from = BestTaxon, values_from = Detected, values_fill = 0) %>%
  mutate(total_sp = rowSums(across(-NWFSCsampleID))) %>%
  column_to_rownames("NWFSCsampleID")

detect_by_station <- mm_detect %>% filter(total_sp > 0) %>% count()

# Add marine mammal detections to sample data
samdf_mm <- merge(samdf, mm_detect, by = "row.names") %>%
  column_to_rownames("Row.names")

sample_data(ps) <- sample_data(samdf_mm)

### QAQC and filter dataset ----------------------------------------------------

# All taxa in seqtab are fish (Actinopteri) or other non-mammal — no class filter needed
ps.fish.sp <- ps

# Remove samples with no reads
ps.fish.sp <- prune_samples(sample_sums(ps.fish.sp) > 0, ps.fish.sp)

# glom to genus
ps.gen <- tax_glom(ps.fish.sp, "genus")

# Transform to proportional space
ps.prop <- transform_sample_counts(ps.fish.sp, function(otu) otu / sum(otu))

ps.prop.gen <- transform_sample_counts(ps.gen, function(otu) otu / sum(otu))

### Biodiversity ---------------------------------------------------------------

plot_richness(ps.prop.gen, x = "depth", measures = c("Shannon", "Simpson"))

ord.nmds.bray <- ordinate(ps.prop.gen, method = "NMDS", distance = "bray")
p <- plot_ordination(ps.prop.gen, ord.nmds.bray, color = "depth", title = "Bray NMDS")

p + geom_text(aes(label = sample_names(ps.prop.gen)), size = 2.5,
              nudge_x = 0.02, nudge_y = 0.02) +
  scale_color_viridis_c()

### Community analysis of variance ---------------------------------------------

# Use only species with at least 10 detection
# target_species <- detect_per_species %>%
#   filter(nDetect > 10) %>%
#   pull(BestTaxon) %>%
#   gsub(" ", "\\.", .)

# Use only species in MURI environmental models
target_species <- c("Berardius.bairdii", "Lagenorhynchus.obliquidens",
                    "Megaptera.novaeangliae", "Mirounga.angustirostris",
                    "Phocoena.phocoena","Zalophus.californianus")

sig_results_list <- list()

for (sp in target_species) {
  
  # Skip if column not in sample data (species with 0 detections excluded above)
  if (!sp %in% colnames(sample_data(ps.gen))) next
  
  # Run ANCOMBC2
  out <- ancombc2(data = ps.gen, tax_level = "genus", fix_formula = sp,
                  p_adj_method = "fdr", struc_zero = FALSE, neg_lb = FALSE,
                  pseudo = 0)
  
  res <- out$res
  colnames(res)[13] <- "diff"
  colnames(res)[17] <- "robust"
  
  sig_res <- res[which(res$diff == TRUE), ]
  
  sig_results_list[[sp]] <- sig_res %>%
    pivot_longer(-taxon, names_to = "metric", values_to = "value") %>%
    mutate(predator = sp)
}

# Combine results
sig_results_all <- dplyr::bind_rows(sig_results_list) %>%
  filter(!grepl("Intercept", metric)) %>%
  separate(metric, into = c("metric", NA)) %>%
  pivot_wider(id_cols = c(taxon, predator), names_from = metric, values_from = value) %>%
  mutate(across(c(predator, taxon), ~ gsub("\\.", " ", .))) %>%
  rename(species = taxon)   # taxon == species name directly (no ASV lookup needed)

muriPreydiff <- ggplot(sig_results_all, aes(x = species, y = lfc, color = predator)) +
  geom_jitter(size = 6, alpha = 0.5, height = 0, width = 0.05) +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90))

png("Figures/MURI_ANCOM_prey.png")
muriPreydiff
dev.off()

### Boxplots of potential prey species -----------------------------------------

otu_de <- as(otu_table(ps.prop.gen), "matrix")

# convert column names when using ps.prop.gen
colnames(otu_de) <- as.data.frame(tax_table(ps.prop.gen)) %>% pull(genus)


prey_sp <- c(unique(sig_results_all$species), "Merluccius", "Thunnus")
pred_sp <- unique(sig_results_all$predator)
pred_cols <- gsub(" ", "\\.", pred_sp)

de_prey <- as.data.frame(otu_de) %>%
  rownames_to_column("NWFSCsampleID") %>%
  select(NWFSCsampleID, any_of(prey_sp)) %>%
  left_join(
    mm_detect %>%
      select(any_of(pred_sp)) %>%          # pred_sp with spaces, not pred_cols
      rownames_to_column("NWFSCsampleID"),
    by = "NWFSCsampleID"
  ) %>%
  pivot_longer(cols = any_of(pred_sp),     # same here
               names_to = "predator", values_to = "detected") %>%
  pivot_longer(cols = any_of(prey_sp),
               names_to = "prey", values_to = "pReads") 
  #semi_join(sig_results_all, by = c("prey" = "species", "predator" = "predator"))

muriPreyBox <- ggplot(de_prey, aes(y = pReads, x = prey, fill = as.factor(detected))) +
  geom_boxplot(outliers = FALSE) +
  facet_wrap("predator", scales = "free", ncol = 2) +
  theme_minimal() +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  theme(strip.text = element_text(size = 14))

png("Figures/MURI_prey_boxplot.png", width = 1500, height = 2000)
muriPreyBox
dev.off()

#### Save data -----------------------------------------------------------------

save(muriPreydiff, muriPreyBox, file = "./data products/MURI_prey_plots.Rdata")
save(ps.fish.sp, ps.prop,
     detect_by_station, detect_by_species,
     mm_detect, de_prey,
     sig_results_all,
     file = "./data products/MURI_prey_exploration.Rdata")

### Combine with detection data ------------------------------------------------

de_prey_wide <- de_prey %>% 
  pivot_wider(names_from = prey, values_from = pReads)

load("./ProcessedData/detect_and_env_muri.Rdata")

detect_data_all_muri <- detect_data_merge_muri %>% 
  st_drop_geometry() %>% 
  dplyr::select(-(Agonidae:Zoarcidae)) %>% 
  group_by(NWFSCsampleID, depth, 
           date, year, month, day, water.depth, bathy.bottom.depth,    
           bottom.depth.consensus, utm.lon, utm.lat, transect.dist.km, 
           control, nReps, lat_deg, lon_deg, x, y, bathy, distShore, slope,
           SSS07, SST07, curVel, Chla, mld_dr003_7, BestTaxon) %>% 
  summarize(Detected = ceiling(mean(Detected, na.rm = TRUE)), .groups = "drop") %>% 
  left_join(de_prey_wide, by = c("NWFSCsampleID" = "NWFSCsampleID", 
                                 "BestTaxon" = "predator",
                                 "Detected" = "detected"))
  
save(env_data_muri, env_df_muri, detect_data_target,
     detect_data_merge_muri, detect_data_all_muri,
     de_prey, de_prey_wide,
     file = "ProcessedData/detect_env_prey_muri.Rdata")

