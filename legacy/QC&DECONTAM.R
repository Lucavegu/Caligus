
######################
# DECONTAMINATION PART
######################

ps <- readRDS(file.path(DIRS$rds, "phyloseq_caligus_microbiome_SILVA.rds"))

suppressPackageStartupMessages({
  library(decontam)
  library(ggplot2)
})
message("decontam version: ", as.character(packageVersion("decontam")))

# Inspect library sizes by group (Sex)
df <- as.data.frame(sample_data(ps))
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize), ]
df$Index <- seq_len(nrow(df))

ggplot(df, aes(x = Index, y = LibrarySize, color = Sex)) +
  geom_point() +
  theme_minimal() +
  labs(title = "Library sizes by sample", x = "Sample index", y = "Reads")

# Identify contaminants using prevalence method
sample_data(ps)$is.neg <- sample_data(ps)$Sex == "Kit"
contamdf.prev05 <- isContaminant(ps, method = "prevalence", neg = "is.neg", threshold = 0.5)
message("Contaminants detected: ")
print(table(contamdf.prev05$contaminant))

# Remove contaminants
ps.noncontam <- prune_taxa(!contamdf.prev05$contaminant, ps)

# Root the tree (required for some phylogenetic metrics)
set.seed(123)
phy_tree(ps.noncontam) <- root(phy_tree(ps.noncontam), sample(taxa_names(ps.noncontam), 1), resolve.root = TRUE)

# Remove kit samples
Ps_kit <- subset_samples(ps.noncontam, Sex == "Kit")
Filtered_ps <- subset_samples(ps.noncontam, Sex != "Kit")

# Remove chloroplast, mitochondria, and eukaryotic sequences
filterPhyla <- c(NA, "Chloroplast", "Mitochondria", "Eukaryota")
Filtered_ps2 <- subset_taxa(Filtered_ps, !Kingdom %in% filterPhyla)
Filtered_ps2 <- subset_taxa(Filtered_ps2, !Phylum %in% filterPhyla)

# Retain taxa present >3 reads in >20% samples
Filtered_ps3 <- filter_taxa(Filtered_ps2, function(x) sum(x > 3) > (0.2 * length(x)), TRUE)

# Summary read counts
message("Total reads (raw): ", sum(sample_sums(ps)))
message("After filtering taxa: ", sum(sample_sums(Filtered_ps2)))
message("After abundance filter: ", sum(sample_sums(Filtered_ps3)))

# Subset to Bacteria only
bacteria_physeq <- subset_taxa(Filtered_ps3, Kingdom == "Bacteria")

# Taxonomic richness summary
tax_table_data <- tax_table(bacteria_physeq)
cat("Unique phyla:", length(unique(tax_table_data[, "Phylum"])), "\n")
cat("Unique classes:", length(unique(tax_table_data[, "Class"])), "\n")
cat("Unique orders:", length(unique(tax_table_data[, "Order"])), "\n")
cat("Unique families:", length(unique(tax_table_data[, "Family"])), "\n")
cat("Unique genera:", length(unique(tax_table_data[, "Genus"])), "\n")

# Rarefaction curve
library(vegan)
mat <- as(t(otu_table(bacteria_physeq)), "matrix")
raremax <- 10000
rarecurve(mat, step = 100, sample = raremax, col = "blue", label = FALSE)

saveRDS(bacteria_physeq, file.path(DIRS$rds, "bacteria_physeq.rds"))

# Split Skin vs Water samples
Ps_water <- subset_samples(bacteria_physeq, Pittag == "CTRL")
Skin_ps <- subset_samples(bacteria_physeq, Pittag != "CTRL")

# Rarefy skin samples
Skin_rare <- rarefy_even_depth(Skin_ps, sample.size = 10000,
                               rngseed = 123, replace = FALSE, trimOTUs = TRUE, verbose = TRUE)
saveRDS(Skin_rare, file.path(DIRS$rds, "Skin_rare_SILVA.rds"))

## Merge Water (unrarefied) + Skin (rarefied)
merged_physeq <- merge_phyloseq(Ps_water, Skin_rare)

# Optionally confirm counts
# sample_sums(merged_physeq)

# Persist QC outputs
saveRDS(ps.noncontam,   file.path(DIRS$rds, "ps_noncontam.rds"))
saveRDS(Filtered_ps2,   file.path(DIRS$rds, "Filtered_ps2.rds"))
saveRDS(Filtered_ps3,   file.path(DIRS$rds, "Filtered_ps3.rds"))
saveRDS(merged_physeq,  file.path(DIRS$rds, "merged_physeq_SILVA.rds"))

message("Decontamination & QC complete. Key artifacts saved in ", DIRS$rds)
