####AMBON MV1 fish data exploration
#### AVC October 2025

library(tidyverse)
library(phyloseq)
library(PNWColors)
library(microViz)
library(ANCOMBC)

### Get data -------------------------------------------------------------------

seqtab <- read.csv("Data/AMBON/ASVtable.csv") %>% 
  separate(1, into = c("sample", "rep")) %>% 
  group_by(sample) %>% 
  summarize(across(c(ASV1:ASV880), ~ceiling(mean(.x, na.rm = TRUE)))) %>% 
  mutate(across(everything(), ~replace_na(., 0))) %>% 
  column_to_rownames("sample")

taxa <- read.csv("Data/AMBON/taxonomy_20250609_collapsed.csv", row.names = 1) %>% 
  select(kingdom,phylum,class,family,genus,species) %>% 
  as.matrix()

samdf <- read.csv("Data/AMBON/FAIRe_noaa-afsc-dbo1.csv", row.names = 1)

### Make phyloseq object -------------------------------------------------------

ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

# remove control samples
ps <- subset_samples(ps, geo_loc_name != "not applicable: control sample")

### Add mm detects to ps object ------------------------------------------------
# get marine mammal detections
ps.mm <- subset_taxa(ps, class == "Mammalia")
ps.mm <- subset_taxa(ps.mm, !(family %in% c("Hominidae", "Suidae", "Bovidae", "Canidae")))
ps.mm.sp <- tax_glom(ps.mm, taxrank = "species")

# convert mm detections to a table
otu_mm <- as(otu_table(ps.mm.sp), "matrix")
tax_table_mm <- as.data.frame(tax_table(ps.mm.sp))
colnames(otu_mm) <- tax_table_mm$species

mm_detect <- as.data.frame((otu_mm > 0) * 1) %>% 
  mutate(total_sp = rowSums(across(everything())))

detect_by_species <- mm_detect %>% 
  select(-total_sp) %>% 
  rownames_to_column("sample") %>% 
  pivot_longer(-sample, names_to = "species", values_to = "detected") %>% 
  group_by(species) %>% 
  summarise(nDetect = sum(detected))

detect_by_station <- mm_detect %>% filter(total_sp > 0) %>% count()

rm(ps.mm, ps.mm.sp, otu_mm, tax_table_mm)

# add marine mammal detections to sample data sheet

samdf_mm <- merge(samdf, mm_detect, by = "row.names") %>% 
  column_to_rownames("Row.names")

sample_data(ps) <- sample_data(samdf_mm)

### QAQC and filter dataset ----------------------------------------------------

# remove non-fish
ps.fish <- subset_taxa(ps, class == "Actinopteri")

# merge to species
ps.fish.sp <- tax_glom(ps.fish, taxrank = "species")

# remove rare species
#ps.fish.sp <- tax_filter(ps.fish.sp, min_prevalence = 5, 
#                    min_sample_abundance = 10)

# remove samples with no reads
ps.fish.sp <- prune_samples(sample_sums(ps.fish.sp) > 0, ps.fish.sp)

# transform to proportional space
ps.prop <- transform_sample_counts(ps.fish.sp, function(otu) otu/sum(otu))

### Biodiversity ---------------------------------------------------------------

plot_richness(ps.prop, x = "geo_loc_name", measures=c("Shannon", "Simpson"))


ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
plot_ordination(ps.prop, ord.nmds.bray, color="geo_loc_name", title="Bray NMDS")

plot_bar(ps.prop, x = "Erignathus.barbatus", fill = "family") +
  theme(legend.position = "none")


### Community analysis of variance----------------------------------------------
# List of target species
target_species <-  gsub(" ", "\\.", detect_by_species$species)
  
sig_results_list <- list()

for (sp in target_species) {
  
  #Run ANCOMBC
  out <- ancombc2(data = ps.fish.sp, tax_level = "species", fix_formula = sp, 
                  p_adj_method = "fdr", struc_zero = FALSE, neg_lb = FALSE, 
                  pseudo = 0)
  
  # Extract results
  res <- out$res
  colnames(res)[13] <- "diff"
  colnames(res)[17] <- "robust"
  
  # Filter significant results
  sig_res <- res[which(res$diff == TRUE),]
  
  # pivot
  sig_results_list[[sp]] <- sig_res %>% 
    pivot_longer(-taxon, names_to = "metric", values_to = "value") %>% 
    mutate(predator = sp) 
  
}

# Combine all results into one data frame
sig_results_all <- dplyr::bind_rows(sig_results_list) %>% 
  filter(!(grepl("Intercept", metric))) %>% 
  separate(metric, into = c("metric", NA)) %>% 
  pivot_wider(id_cols = c(taxon, predator), names_from = metric, values_from = value) %>% 
  mutate(across(c(predator, taxon), ~ gsub("\\.", " ", .))) %>%
  left_join(ps.fish.sp@tax_table %>% as.data.frame() %>% rownames_to_column("ASV") %>% select(ASV,species), by = c("taxon" = "ASV"))

mv1Preydiff <- ggplot(sig_results_all, aes(x = species, y = lfc, color = predator)) +
                  geom_jitter(size = 6, alpha = 0.5, height = 0, width = 0.05) +
                  theme_minimal() +
                  theme(legend.position = "bottom")

png("Figures/MV1_ANCOM_prey.png")
mv1Preydiff
dev.off()

### Boxplots of potential prey species -----------------------------------------
# convert detections to a table
otu_de <- as(otu_table(ps.prop), "matrix")
tax_table_de <- as.data.frame(tax_table(ps.prop))
colnames(otu_de) <- tax_table_de$species

de_prey <- as.data.frame(otu_de) %>% 
  rownames_to_column("sample") %>% 
  select(sample, unique(sig_results_all$species)) %>% 
  left_join(mm_detect %>% rownames_to_column("sample"), by = "sample") %>% 
  pivot_longer(12:24, names_to = "predator", values_to = "detected") %>% 
  select(-total_sp) %>% 
  pivot_longer(2:11, names_to = "prey", values_to = "pReads") %>% 
  semi_join(sig_results_all, by = c("prey" = "species", "predator" = "predator"))

mv1PreyBox <- ggplot(de_prey, aes(y = pReads, x = prey, fill = as.factor(detected))) +
                geom_boxplot(outliers = FALSE) +
                facet_wrap("predator", scales = "free", ncol = 2) +
                theme_minimal() +
                scale_x_discrete(guide = guide_axis(n.dodge = 2)) 

png("Figures/MV1_prey_boxplot.png")
mv1PreyBox
dev.off()

#### Save data -----------------------------------------------------------------

save(mv1Preydiff, mv1PreyBox, file = "./data products/MV1_prey_plots.Rdata")
save(ps.fish.sp, ps.prop, 
     detect_by_station, detect_by_species, 
     mm_detect, de_prey,
     sig_results_all,
     file = "./data products/MV1_prey_exploration.Rdata")
