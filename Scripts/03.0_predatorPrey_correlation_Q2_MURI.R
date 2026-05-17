#### MURI fish ~ marine mammal predator-prey analysis
#### Adapted from 03_0_predatorPrey_correlation_Q2_AMBON.R

library(tidyverse)
library(phyloseq)
library(PNWColors)
library(microViz)
library(ANCOMBC)

### Get data -------------------------------------------------------------------

load("ProcessedData/detect_data_muri.Rdata")

# Identify species (fish/other) columns — everything after the metadata block
species_cols <- colnames(detect_data_muri)[(which(colnames(detect_data_muri) == "nReps") + 1):ncol(detect_data_muri)]

# Build seqtab: one row per SampleUID (species cols are identical across the 21
# BestTaxon rows for each SampleUID, so take the first), then collapse
# dilution/techRep replicates within NWFSCsampleID via ceiling(mean()),
# matching the AMBON rep-collapse approach.

seqtab <- detect_data_muri %>%
  distinct(SampleUID, NWFSCsampleID, .keep_all = TRUE) %>%
  select(NWFSCsampleID, all_of(species_cols)) %>%
  group_by(NWFSCsampleID) %>%
  summarize(across(all_of(species_cols), ~ ceiling(mean(.x, na.rm = TRUE)))) %>%
  mutate(across(everything(), ~ replace_na(., 0))) %>%
  column_to_rownames("NWFSCsampleID")

# Build samdf: one row per NWFSCsampleID (take first rep for metadata)
samdf <- detect_data_muri %>%
  distinct(NWFSCsampleID, .keep_all = TRUE) %>%
  select(NWFSCsampleID, sample, station, Niskin, depth, transect, lat, lon,
         water.depth, bathy.bottom.depth, bottom.depth.consensus,
         year, month, day, date, volume, Fluor, Zymo,
         control, drop.sample, field.negative.type, totalReads, nReps) %>%
  column_to_rownames("NWFSCsampleID")

# Build taxa table: BestTaxon is already at species level in this dataset.
# Construct a minimal tax_table from detect_per_species names (all Mammalia).
# For fish we don't have a taxonomy file — use species names as both genus and species.
taxa_fish <- data.frame(
  species = species_cols,
  genus   = word(species_cols, 1),
  family  = NA_character_,
  class   = NA_character_,
  phylum  = NA_character_,
  kingdom = NA_character_,
  row.names = species_cols
) %>%
  as.matrix()

### Make phyloseq object -------------------------------------------------------

ps <- phyloseq(otu_table(seqtab, taxa_are_rows = FALSE),
               sample_data(samdf),
               tax_table(taxa_fish))

### Add mm detects to ps object ------------------------------------------------

# Collapse marine mammal detections: one row per NWFSCsampleID per BestTaxon,
# then take ceiling(mean(Detected)) across dilution/techRep reps
mm_detect_long <- detect_data_muri %>%
  select(NWFSCsampleID, BestTaxon, Detected) %>%
  group_by(NWFSCsampleID, BestTaxon) %>%
  summarize(Detected = ceiling(mean(Detected, na.rm = TRUE)), .groups = "drop")

mm_detect <- mm_detect_long %>%
  pivot_wider(names_from = BestTaxon, values_from = Detected, values_fill = 0) %>%
  mutate(total_sp = rowSums(across(-NWFSCsampleID))) %>%
  column_to_rownames("NWFSCsampleID")

# detect_per_species is already in the Rdata (equivalent to detect_by_species in AMBON)
detect_by_species <- detect_per_species %>%
  rename(species = BestTaxon) %>%
  arrange(-nDetect)

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

# Transform to proportional space
ps.prop <- transform_sample_counts(ps.fish.sp, function(otu) otu / sum(otu))

### Biodiversity ---------------------------------------------------------------

plot_richness(ps.prop, x = "depth", measures = c("Shannon", "Simpson"))

ord.nmds.bray <- ordinate(ps.prop, method = "NMDS", distance = "bray")
p <- plot_ordination(ps.prop, ord.nmds.bray, color = "depth", title = "Bray NMDS")

p + geom_text(aes(label = sample_names(ps.prop)), size = 2.5,
              nudge_x = 0.02, nudge_y = 0.02)

### Community analysis of variance ---------------------------------------------

# Use only species with at least 1 detection
target_species <- detect_by_species %>%
  filter(nDetect > 0) %>%
  pull(species) %>%
  gsub(" ", "\\.", .)

sig_results_list <- list()

for (sp in target_species) {
  
  # Skip if column not in sample data (species with 0 detections excluded above)
  if (!sp %in% colnames(sample_data(ps.fish.sp))) next
  
  # Run ANCOMBC2
  out <- ancombc2(data = ps.fish.sp, tax_level = NULL, fix_formula = sp,
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

otu_de <- as(otu_table(ps.prop), "matrix")
# column names are already species names
colnames(otu_de) <- colnames(seqtab)[colnames(seqtab) %in% colnames(otu_de)]

prey_sp <- unique(sig_results_all$species)
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
               names_to = "prey", values_to = "pReads") %>%
  semi_join(sig_results_all, by = c("prey" = "species", "predator" = "predator"))

muriPreyBox <- ggplot(de_prey, aes(y = pReads, x = prey, fill = as.factor(detected))) +
  geom_boxplot(outliers = FALSE) +
  facet_wrap("predator", scales = "free", ncol = 2) +
  theme_minimal() +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))

png("Figures/MURI_prey_boxplot.png")
muriPreyBox
dev.off()

#### Save data -----------------------------------------------------------------

save(muriPreydiff, muriPreyBox, file = "./data products/MURI_prey_plots.Rdata")
save(ps.fish.sp, ps.prop,
     detect_by_station, detect_by_species,
     mm_detect, de_prey,
     sig_results_all,
     file = "./data products/MURI_prey_exploration.Rdata")