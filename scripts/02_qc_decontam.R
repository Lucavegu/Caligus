#!/usr/bin/env Rscript

###############################################################################
# Skin microbiome & Caligus: QC, decontamination, filtering and rarefaction
# Input : phyloseq_caligus_microbiome_SILVA.rds (from the preprocessing script)
# Output: RDS objects for the later scripts, plus tables and plots in logs/qc/
# Author: Lucas Venegas (adapted & curated; revised 2026-10)
#
# What this version does differently from the original:
#   - sequencing batch (run1 / run2) is added to the sample data and checked
#   - the prevalence filter is applied to skin and to water separately
#   - chloroplast and mitochondria are removed at the correct SILVA ranks
#   - graphics compare the negative controls (and the ASVs flagged as
#     contaminants) with the skin samples
#   - contaminant list, control read counts and per-sample depths are saved
#   - rarefaction is repeated many times to measure how much one random draw
#     matters, and alternative normalizations (proportions, Hellinger, CLR,
#     robust CLR, optional CSS) are produced and compared
###############################################################################

## --------------------------- User Parameters ---------------------------------
PROJECT_DIR <- getwd()
RDS_DIR     <- file.path(PROJECT_DIR)                 # where the phyloseq RDS files live
QC_DIR      <- file.path(PROJECT_DIR, "logs", "qc")   # tables and plots from this script
if (!dir.exists(QC_DIR)) dir.create(QC_DIR, recursive = TRUE)

# How controls are labelled in the metadata
CTRL_PITTAG <- "CTRL"   # Pittag value of water controls
KIT_SEX     <- "Kit"    # Sex value of extraction-kit negative controls

# Sequencing run. If the phyloseq object already has a seq_batch column (written
# by the preprocessing script) it is used; otherwise the samples listed here
# are "run2" and every other sample is "run1".
BATCH2_SAMPLES <- c("SG21_0185", "SG21_0186", "SG21_0190",
                    "SG21_0369", "SG21_0371", "SG21_0374")

# decontam (prevalence method)
DECONTAM_THRESHOLD <- 0.5

# TRUE removes Chloroplast (Order) and Mitochondria (Family); FALSE reproduces
# the original filter, which did not remove them
FIX_ORGANELLE_FILTER <- TRUE

# Prevalence filter, applied to each component separately:
# keep taxa with more than MIN_READS reads in more than MIN_PREV of the samples
SKIN_MIN_READS  <- 3
SKIN_MIN_PREV   <- 0.2
WATER_MIN_READS <- 3
WATER_MIN_PREV  <- 0.2    # with 5 water samples this means at least 2 samples

# Rarefaction of skin samples
RAREFY_DEPTH <- 10000
RAREFY_SEED  <- 123

# Normalization alternatives
N_RAREFY_ITER   <- 100    # repeated rarefactions
CLR_PSEUDOCOUNT <- 0.5    # added to every count before the CLR (zeros cannot be logged)
N_PERM          <- 999    # permutations for the PERMANOVA tests
PERM_SEED       <- 123

# Tree rooting: "random_asv" reproduces the original analysis; "midpoint" does
# not depend on a random draw, but it changes Faith's PD values
ROOT_METHOD <- "random_asv"
ROOT_SEED   <- 123

## ------------------------------ Libraries ------------------------------------
suppressPackageStartupMessages({
  library(phyloseq)
  library(decontam)
  library(ggplot2)
  library(vegan)
  library(ape)
  library(phangorn)
})
message("decontam version: ", as.character(packageVersion("decontam")))

## ------------------------------- Helpers -------------------------------------
otu_mat <- function(x) {            # counts with ASVs in rows
  m <- as(otu_table(x), "matrix")
  if (!taxa_are_rows(x)) m <- t(m)
  m
}
tax_mat <- function(x) as(tax_table(x), "matrix")
write_tsv <- function(x, name) {
  write.table(x, file.path(QC_DIR, name), sep = "\t", quote = FALSE, row.names = FALSE)
}
stage_row <- function(stage, x) {
  data.frame(stage = stage, n_samples = nsamples(x), n_ASVs = ntaxa(x), total_reads = sum(sample_sums(x)))
}
to_num <- function(v) suppressWarnings(as.numeric(as.character(v)))

## -------------------------------- Input --------------------------------------
ps <- readRDS(file.path(RDS_DIR, "phyloseq_caligus_microbiome_SILVA.rds"))

sd0 <- as(sample_data(ps), "data.frame")
is_kit   <- sd0$Sex %in% KIT_SEX
is_water <- sd0$Pittag %in% CTRL_PITTAG & !is_kit
is_skin  <- !is_kit & !is_water
sample_type <- ifelse(is_kit, "kit", ifelse(is_water, "water", "skin"))
names(sample_type) <- sample_names(ps)

if ("seq_batch" %in% colnames(sd0)) {
  seq_batch <- as.character(sd0$seq_batch)
  message("Sequencing run taken from the seq_batch column of the phyloseq object.")
} else {
  not_found <- setdiff(BATCH2_SAMPLES, sample_names(ps))
  if (length(not_found) > 0) warning("BATCH2_SAMPLES not in the data: ", paste(not_found, collapse = ", "))
  seq_batch <- ifelse(sample_names(ps) %in% BATCH2_SAMPLES, "run2", "run1")
}
names(seq_batch) <- sample_names(ps)

sample_data(ps)$sample_type <- sample_type
sample_data(ps)$seq_batch   <- seq_batch
sample_data(ps)$is.neg      <- is_kit

message("Samples by type and sequencing run:")
print(table(sample_type, seq_batch))
if (sum(is_kit) == 0) stop("No kit negative controls found (Sex == '", KIT_SEX, "').")
runs_without_neg <- setdiff(unique(seq_batch), unique(seq_batch[is_kit]))
if (length(runs_without_neg) > 0) {
  message("Note: no negative control in ", paste(runs_without_neg, collapse = ", "),
          "; contaminants are identified from the controls of the other run and removed from all samples.")
}

stages <- stage_row("1_input", ps)

## --------------------------- Library sizes -----------------------------------
df <- data.frame(sample = sample_names(ps), sample_type = sample_type, seq_batch = seq_batch,
                 LibrarySize = sample_sums(ps))
df <- df[order(df$LibrarySize), ]
df$Index <- seq_len(nrow(df))

p_lib <- ggplot(df, aes(x = Index, y = LibrarySize, color = sample_type, shape = seq_batch)) +
  geom_point(size = 2) +
  theme_minimal() +
  labs(title = "Library sizes by sample", x = "Sample index", y = "Reads",
       color = "Sample type", shape = "Sequencing run")
ggsave(file.path(QC_DIR, "library_sizes.pdf"), p_lib, width = 7, height = 5)

## ----------------------------- Decontamination -------------------------------
# Prevalence method: ASVs more prevalent in kit negative controls than in samples
contamdf <- isContaminant(ps, method = "prevalence", neg = "is.neg", threshold = DECONTAM_THRESHOLD)
contamdf <- contamdf[taxa_names(ps), , drop = FALSE]
message("Contaminants detected: ")
print(table(contamdf$contaminant))

# Counts and relative abundances before any removal
m0   <- otu_mat(ps)[taxa_names(ps), , drop = FALSE]
rel0 <- sweep(m0, 2, pmax(colSums(m0), 1), "/")
tax0 <- tax_mat(ps)[taxa_names(ps), , drop = FALSE]

# One row per ASV: decontam result, presence in controls and skin, taxonomy
asv_info <- data.frame(
  ASV            = taxa_names(ps),
  decontam_p     = contamdf$p,
  contaminant    = contamdf$contaminant,
  n_kits_present = rowSums(m0[, is_kit, drop = FALSE] > 0),
  reads_in_kits  = rowSums(m0[, is_kit, drop = FALSE]),
  mean_rel_kits  = rowMeans(rel0[, is_kit, drop = FALSE]),
  n_skin_present = rowSums(m0[, is_skin, drop = FALSE] > 0),
  skin_prev_pct  = round(100 * rowMeans(m0[, is_skin, drop = FALSE] > 0), 1),
  reads_in_skin  = rowSums(m0[, is_skin, drop = FALSE]),
  mean_rel_skin  = rowMeans(rel0[, is_skin, drop = FALSE]),
  reads_in_water = rowSums(m0[, is_water, drop = FALSE]),
  tax0,
  row.names = NULL, check.names = FALSE
)

# Read counts of every negative control
write_tsv(data.frame(sample = sample_names(ps)[is_kit],
                     seq_batch = seq_batch[is_kit],
                     reads  = sample_sums(ps)[is_kit],
                     n_ASVs = colSums(m0[, is_kit, drop = FALSE] > 0)),
          "negative_controls_summary.tsv")

# Remove contaminants
keep_taxa    <- rownames(contamdf)[!contamdf$contaminant]
ps.noncontam <- prune_taxa(keep_taxa, ps)
stages <- rbind(stages, stage_row("2_after_decontam", ps.noncontam))

## ------------------------------ Root the tree --------------------------------
# A rooted tree is required for some phylogenetic metrics (Faith's PD)
if (ROOT_METHOD == "midpoint") {
  phy_tree(ps.noncontam) <- phangorn::midpoint(phy_tree(ps.noncontam))
  root_note <- "midpoint"
} else {
  set.seed(ROOT_SEED)
  outgroup <- sample(taxa_names(ps.noncontam), 1)
  phy_tree(ps.noncontam) <- ape::root(phy_tree(ps.noncontam), outgroup, resolve.root = TRUE)
  root_note <- paste0("random ASV (", outgroup, ", seed ", ROOT_SEED, ")")
}
message("Tree rooted: ", root_note)

## ---------------------------- Remove kit samples -----------------------------
st <- sample_data(ps.noncontam)$sample_type
Filtered_ps <- prune_samples(sample_names(ps.noncontam)[st != "kit"], ps.noncontam)

## ------------------------------ Taxonomy filter ------------------------------
# Remove ASVs with no Kingdom or Phylum, eukaryotes, chloroplasts and mitochondria
tt <- tax_mat(Filtered_ps)
filterPhyla <- c(NA, "Chloroplast", "Mitochondria", "Eukaryota")
drop_taxa <- tt[, "Kingdom"] %in% filterPhyla | tt[, "Phylum"] %in% filterPhyla
if (FIX_ORGANELLE_FILTER) {
  drop_taxa <- drop_taxa | tt[, "Order"] %in% "Chloroplast" | tt[, "Family"] %in% "Mitochondria"
}
Filtered_ps2 <- prune_taxa(rownames(tt)[!drop_taxa], Filtered_ps)

# Bacteria only (this also removes Archaea)
tt2 <- tax_mat(Filtered_ps2)
bacteria_all <- prune_taxa(rownames(tt2)[tt2[, "Kingdom"] %in% "Bacteria"], Filtered_ps2)
stages <- rbind(stages, stage_row("3_no_kits_taxonomy_filter_bacteria", bacteria_all))

## ------------------ Prevalence filter: skin and water separately -------------
stb         <- sample_data(bacteria_all)$sample_type
skin_names  <- sample_names(bacteria_all)[stb == "skin"]
water_names <- sample_names(bacteria_all)[stb == "water"]
mb <- otu_mat(bacteria_all)

pass_skin  <- rowSums(mb[, skin_names, drop = FALSE] > SKIN_MIN_READS) > (SKIN_MIN_PREV * length(skin_names))
pass_water <- if (length(water_names) > 0) {
  rowSums(mb[, water_names, drop = FALSE] > WATER_MIN_READS) > (WATER_MIN_PREV * length(water_names))
} else {
  setNames(rep(FALSE, nrow(mb)), rownames(mb))
}
# The original filter (skin and water pooled), kept only for comparison
pass_pooled <- rowSums(mb > SKIN_MIN_READS) > (SKIN_MIN_PREV * ncol(mb))

message("Taxa kept - skin filter: ", sum(pass_skin), "; water filter: ", sum(pass_water),
        "; kept by both: ", sum(pass_skin & pass_water),
        "; original pooled filter: ", sum(pass_pooled))
message("Skin taxa gained vs pooled filter: ", sum(pass_skin & !pass_pooled),
        "; lost: ", sum(pass_pooled & !pass_skin))

# Skin : skin samples x taxa that pass the skin filter  (input for diversity and ML)
# Water: water samples x taxa that pass the water filter (describes the water community)
Skin_ps  <- prune_taxa(rownames(mb)[pass_skin],  prune_samples(skin_names,  bacteria_all))
Ps_water <- prune_taxa(rownames(mb)[pass_water], prune_samples(water_names, bacteria_all))

# Skin + water together, with every taxon kept by either filter and the original
# counts. Later scripts take the water samples from this object.
bacteria_physeq <- prune_taxa(rownames(mb)[pass_skin | pass_water], bacteria_all)

stages <- rbind(stages, stage_row("4a_skin_prevalence_filter", Skin_ps),
                        stage_row("4b_water_prevalence_filter", Ps_water))

# Which filter keeps which taxon
any_pass <- pass_skin | pass_water | pass_pooled
tb <- tax_mat(bacteria_all)
write_tsv(data.frame(ASV = rownames(mb)[any_pass],
                     kept_skin_filter   = pass_skin[any_pass],
                     kept_water_filter  = pass_water[any_pass],
                     kept_pooled_filter = pass_pooled[any_pass],
                     skin_prev_pct  = round(100 * rowMeans(mb[any_pass, skin_names, drop = FALSE] > 0), 1),
                     water_n_present = rowSums(mb[any_pass, water_names, drop = FALSE] > 0),
                     skin_reads  = rowSums(mb[any_pass, skin_names, drop = FALSE]),
                     water_reads = rowSums(mb[any_pass, water_names, drop = FALSE]),
                     tb[any_pass, , drop = FALSE], check.names = FALSE),
          "prevalence_filter_by_component.tsv")

# Taxonomic richness summary (skin)
tax_table_data <- tax_mat(Skin_ps)
for (rk in c("Phylum", "Class", "Order", "Family", "Genus")) {
  cat("Skin, unique ", tolower(rk), ": ", length(unique(tax_table_data[, rk])), "\n", sep = "")
}
message("Water sample depths after the water filter: ", paste(sample_sums(Ps_water), collapse = ", "))

## ---------------------------- Rarefaction curve ------------------------------
pdf(file.path(QC_DIR, "rarefaction_curves_skin.pdf"), width = 7, height = 5)
rarecurve(t(otu_mat(Skin_ps)), step = 100, sample = RAREFY_DEPTH, col = "blue", label = FALSE)
dev.off()

## --------------------------- Rarefy skin samples -----------------------------
skin_depth <- sample_sums(Skin_ps)

# How many skin samples would be kept at other depths (to justify RAREFY_DEPTH)
alt_depths <- sort(unique(c(1000, 2000, 5000, 7500, 10000, 15000, 20000, RAREFY_DEPTH)))
keep_tab <- data.frame(depth = alt_depths,
                       skin_samples_kept = sapply(alt_depths, function(d) sum(skin_depth >= d)),
                       skin_samples_lost = sapply(alt_depths, function(d) sum(skin_depth < d)))
keep_tab$samples_kept_pct <- round(100 * keep_tab$skin_samples_kept / length(skin_depth), 1)
keep_tab$reads_kept_pct   <- round(100 * keep_tab$skin_samples_kept * keep_tab$depth / sum(skin_depth), 1)
write_tsv(keep_tab, "samples_kept_by_rarefaction_depth.tsv")

p_depth <- ggplot(data.frame(depth = skin_depth, seq_batch = seq_batch[names(skin_depth)]),
                  aes(x = depth, fill = seq_batch)) +
  geom_histogram(bins = 30) +
  geom_vline(xintercept = RAREFY_DEPTH, linetype = "dashed") +
  theme_minimal() +
  labs(title = "Skin sample depth after filtering", x = "Reads", y = "Samples", fill = "Sequencing run")
ggsave(file.path(QC_DIR, "skin_depth_distribution.pdf"), p_depth, width = 7, height = 5)

Skin_rare <- rarefy_even_depth(Skin_ps, sample.size = RAREFY_DEPTH,
                               rngseed = RAREFY_SEED, replace = FALSE, trimOTUs = TRUE, verbose = TRUE)
stages <- rbind(stages, stage_row("5_skin_rarefied", Skin_rare))

lost <- setdiff(sample_names(Skin_ps), sample_names(Skin_rare))
message("Skin samples removed by rarefaction (below ", RAREFY_DEPTH, " reads): ",
        length(lost), if (length(lost)) paste0(" (", paste(lost, collapse = ", "), ")") else "")

## ----------------- Water (unrarefied) + Skin (rarefied) ----------------------
# Built from the common object so that tree, taxonomy and sequences stay aligned:
# water keeps its original counts, skin counts are replaced by the rarefied ones.
merged_taxa   <- union(taxa_names(Skin_rare), rownames(mb)[pass_water])
merged_physeq <- prune_taxa(merged_taxa, prune_samples(c(water_names, sample_names(Skin_rare)), bacteria_all))
mm <- otu_mat(merged_physeq)
mm[, sample_names(Skin_rare)] <- 0
mm[taxa_names(Skin_rare), sample_names(Skin_rare)] <- otu_mat(Skin_rare)[taxa_names(Skin_rare), sample_names(Skin_rare)]
otu_table(merged_physeq) <- otu_table(mm, taxa_are_rows = TRUE)

## ------------------------- Tables for the reviewers --------------------------
# Per-sample read depth at each stage
depth_tab <- data.frame(sample = sample_names(ps), sample_type = sample_type, seq_batch = seq_batch,
                        raw = sample_sums(ps),
                        after_decontam = sample_sums(ps.noncontam)[sample_names(ps)],
                        row.names = sample_names(ps))
depth_tab$contaminant_read_pct <- round(100 * (depth_tab$raw - depth_tab$after_decontam) / pmax(1, depth_tab$raw), 2)
depth_tab$after_taxonomy_filter   <- NA_real_
depth_tab$after_prevalence_filter <- NA_real_
depth_tab[sample_names(bacteria_all), "after_taxonomy_filter"] <- sample_sums(bacteria_all)
depth_tab[sample_names(Skin_ps),  "after_prevalence_filter"]   <- sample_sums(Skin_ps)
depth_tab[sample_names(Ps_water), "after_prevalence_filter"]   <- sample_sums(Ps_water)
depth_tab$kept_after_rarefaction <- ifelse(depth_tab$sample_type == "skin",
                                           depth_tab$sample %in% sample_names(Skin_rare), NA)
write_tsv(depth_tab, "per_sample_depth_by_stage.tsv")

## ---------------- Rarefaction depth: evidence for the reviewers --------------
# Richness at the rarefaction depth compared with richness at full depth, and
# the slope of the rarefaction curve at that depth, for every sample that
# reaches RAREFY_DEPTH. Done before and after the prevalence filter.
skin_unf <- t(otu_mat(prune_samples(skin_names, bacteria_all)))
skin_unf <- skin_unf[, colSums(skin_unf) > 0, drop = FALSE]
skin_flt <- t(otu_mat(Skin_ps))

rich_tab <- function(mat, label) {
  ok  <- rowSums(mat) >= RAREFY_DEPTH
  sub <- mat[ok, , drop = FALSE]
  obs <- rowSums(sub > 0)
  ex  <- as.numeric(vegan::rarefy(sub, sample = RAREFY_DEPTH))
  slp <- as.numeric(vegan::rareslope(sub, sample = RAREFY_DEPTH))
  data.frame(dataset = label, sample = rownames(sub), depth = rowSums(sub),
             richness_full_depth = obs,
             richness_at_rarefaction_depth = round(ex, 1),
             pct_richness_captured = round(100 * ex / obs, 1),
             new_ASVs_per_1000_extra_reads = round(1000 * slp, 2))
}
rich <- rbind(rich_tab(skin_unf, "Before prevalence filter"),
              rich_tab(skin_flt, "After prevalence filter"))
write_tsv(rich, "rarefaction_richness_captured.tsv")
rich_summary <- do.call(rbind, lapply(split(rich, rich$dataset), function(d) data.frame(
  dataset = d$dataset[1], n_samples = nrow(d),
  median_pct_richness_captured = median(d$pct_richness_captured),
  min_pct_richness_captured = min(d$pct_richness_captured),
  median_new_ASVs_per_1000_extra_reads = median(d$new_ASVs_per_1000_extra_reads))))
write_tsv(rich_summary, "rarefaction_richness_captured_summary.tsv")
print(rich_summary, row.names = FALSE)

# The samples that rarefaction removes, with their depth at every stage
lost_tab <- depth_tab[depth_tab$sample %in% lost, , drop = FALSE]
lost_tab$Total_caligus <- to_num(sd0[lost_tab$sample, "Total_caligus"])
write_tsv(lost_tab, "samples_removed_by_rarefaction.tsv")

keep_long <- rbind(data.frame(depth = keep_tab$depth, measure = "Skin samples kept (%)", pct = keep_tab$samples_kept_pct),
                   data.frame(depth = keep_tab$depth, measure = "Skin reads kept (%)",   pct = keep_tab$reads_kept_pct))
p_keep <- ggplot(keep_long, aes(x = depth, y = pct, color = measure)) +
  geom_line() + geom_point() +
  geom_vline(xintercept = RAREFY_DEPTH, linetype = "dashed") +
  theme_minimal() +
  labs(title = "Samples and reads kept at each candidate rarefaction depth",
       x = "Rarefaction depth (reads per sample)", y = "%", color = NULL)
p_rich <- ggplot(rich, aes(x = dataset, y = pct_richness_captured)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.2) +
  theme_minimal() +
  labs(title = paste0("ASV richness captured at ", RAREFY_DEPTH, " reads, relative to full depth"),
       x = NULL, y = "% of the sample's full-depth richness")

pdf(file.path(QC_DIR, "rarefaction_justification.pdf"), width = 8, height = 6)
print(p_depth)
rarecurve(skin_unf, step = 500, sample = RAREFY_DEPTH, col = "grey30", label = FALSE,
          main = "Rarefaction curves, skin samples, before prevalence filter", xlab = "Reads", ylab = "ASVs")
rarecurve(skin_flt, step = 500, sample = RAREFY_DEPTH, col = "grey30", label = FALSE,
          main = "Rarefaction curves, skin samples, after prevalence filter", xlab = "Reads", ylab = "ASVs")
print(p_keep)
print(p_rich)
dev.off()

# ASVs flagged as contaminants, and every ASV seen in a negative control
asv_info$in_final_skin  <- asv_info$ASV %in% taxa_names(Skin_ps)
asv_info$in_final_water <- asv_info$ASV %in% taxa_names(Ps_water)
write_tsv(asv_info[asv_info$contaminant, ][order(-asv_info$reads_in_skin[asv_info$contaminant]), ],
          "contaminant_ASVs_removed.tsv")
neg_asvs <- asv_info[asv_info$reads_in_kits > 0, ]
write_tsv(neg_asvs[order(-neg_asvs$reads_in_kits), ], "ASVs_detected_in_negative_controls.tsv")
message("ASVs seen in negative controls that remain in the final skin dataset: ", sum(neg_asvs$in_final_skin))

###############################################################################
# GRAPHICS 1: negative controls and contaminant ASVs compared with skin samples
# (all computed on the data before any removal)
###############################################################################
kit_asv   <- asv_info$reads_in_kits > 0
skin_asv  <- asv_info$reads_in_skin > 0
is_contam <- asv_info$contaminant
rho <- suppressWarnings(cor(asv_info$mean_rel_skin[kit_asv], asv_info$mean_rel_kits[kit_asv], method = "spearman"))

# Numbers behind the plots
overlap <- data.frame(
  measure = c("Kit negative controls (n)",
              "ASVs detected in kit controls",
              "...of which also detected in skin",
              "% of kit-control reads in ASVs also found in skin",
              "% of skin reads in ASVs detected in kit controls",
              "ASVs flagged as contaminants",
              "% of skin reads in flagged ASVs",
              "Spearman rho, mean abundance in kits vs skin (ASVs detected in kits)"),
  value = c(sum(is_kit),
            sum(kit_asv),
            sum(kit_asv & skin_asv),
            round(100 * sum(m0[kit_asv & skin_asv, is_kit]) / max(1, sum(m0[, is_kit])), 2),
            round(100 * sum(m0[kit_asv, is_skin]) / max(1, sum(m0[, is_skin])), 2),
            sum(is_contam),
            round(100 * sum(m0[is_contam, is_skin]) / max(1, sum(m0[, is_skin])), 2),
            round(rho, 3))
)
write_tsv(overlap, "kit_vs_skin_overlap_summary.tsv")
print(overlap)

status <- ifelse(asv_info$contaminant, "Flagged contaminant (removed)",
                 ifelse(asv_info$in_final_skin, "Kept, in final skin dataset", "Kept, removed by later filters"))
asv_info$status <- factor(status, levels = c("Flagged contaminant (removed)", "Kept, in final skin dataset",
                                             "Kept, removed by later filters"))
status_cols <- c("Flagged contaminant (removed)" = "#D55E00", "Kept, in final skin dataset" = "#0072B2",
                 "Kept, removed by later filters" = "grey70")

# (a) decontam score distribution
p_a <- ggplot(asv_info[!is.na(asv_info$decontam_p), ], aes(x = decontam_p)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = DECONTAM_THRESHOLD, linetype = "dashed") +
  theme_minimal() +
  labs(title = "decontam scores (prevalence method)",
       subtitle = paste0("ASVs left of the dashed line (score < ", DECONTAM_THRESHOLD, ") are flagged as contaminants"),
       x = "decontam score", y = "ASVs")

# (b) prevalence in skin against presence in the negative controls
p_b <- ggplot(asv_info, aes(x = factor(n_kits_present), y = skin_prev_pct, color = status)) +
  geom_jitter(width = 0.25, height = 0, alpha = 0.6, size = 1.3) +
  scale_color_manual(values = status_cols) +
  theme_minimal() +
  labs(title = "Prevalence in skin vs presence in negative controls",
       x = "Number of kit controls in which the ASV was detected",
       y = "Skin samples with the ASV (%)", color = NULL)

# (c) abundance in the negative controls against abundance in skin
eps <- 1e-6
p_c <- ggplot(asv_info[kit_asv, ], aes(x = log10(mean_rel_skin + eps), y = log10(mean_rel_kits + eps), color = status)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.7, size = 1.6) +
  scale_color_manual(values = status_cols) +
  theme_minimal() +
  labs(title = "ASVs detected in negative controls: abundance in controls vs skin",
       subtitle = paste0("Spearman rho = ", round(rho, 2), "; points on the dashed line are equally abundant in both"),
       x = "log10 mean relative abundance in skin", y = "log10 mean relative abundance in kit controls", color = NULL)

# (d) ordination of all samples (Bray-Curtis on relative abundance)
d0 <- vegdist(t(rel0), method = "bray")
pc <- cmdscale(d0, k = 2, eig = TRUE)
ve <- round(100 * pc$eig[1:2] / sum(pc$eig[pc$eig > 0]), 1)
pc_df <- data.frame(sample = sample_names(ps), PC1 = pc$points[, 1], PC2 = pc$points[, 2],
                    sample_type = sample_type, seq_batch = seq_batch)
p_d <- ggplot(pc_df, aes(x = PC1, y = PC2, color = sample_type, shape = seq_batch)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_text(data = pc_df[pc_df$sample_type != "skin", ], aes(label = sample), size = 3, vjust = -1, show.legend = FALSE) +
  theme_minimal() +
  labs(title = "PCoA of all samples before decontamination (Bray-Curtis)",
       x = paste0("PCoA1 (", ve[1], "%)"), y = paste0("PCoA2 (", ve[2], "%)"),
       color = "Sample type", shape = "Sequencing run")

# (e) how far each skin sample is from the controls, from water and from other skin samples
D <- as.matrix(d0)
dist_df <- rbind(
  data.frame(comparison = "Skin vs kit controls", bray = rowMeans(D[is_skin, is_kit, drop = FALSE])),
  data.frame(comparison = "Skin vs other skin",   bray = rowSums(D[is_skin, is_skin, drop = FALSE]) / (sum(is_skin) - 1))
)
if (any(is_water)) {
  dist_df <- rbind(dist_df, data.frame(comparison = "Skin vs water", bray = rowMeans(D[is_skin, is_water, drop = FALSE])))
}
p_e <- ggplot(dist_df, aes(x = comparison, y = bray)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
  theme_minimal() +
  labs(title = "Mean Bray-Curtis dissimilarity of each skin sample",
       subtitle = "If skin is as close to the kit controls as to other skin, the controls mostly contain sample-derived reads",
       x = NULL, y = "Bray-Curtis dissimilarity (0 = identical, 1 = nothing shared)")

# (f) genus composition of each kit control next to water and skin
genus <- tax0[, "Genus"]
genus[is.na(genus)] <- "Unclassified"
gen_rel <- rowsum(rel0, group = genus)
grp <- ifelse(is_kit, sample_names(ps),
              ifelse(is_water, "Water (mean)", paste0("Skin ", seq_batch, " (mean)")))
grp_levels <- unique(c(sample_names(ps)[is_kit], "Water (mean)", sort(unique(grp[is_skin]))))
grp_levels <- grp_levels[grp_levels %in% grp]
grp_mean <- sapply(grp_levels, function(g) rowMeans(gen_rel[, grp == g, drop = FALSE]))
top_gen  <- names(sort(rowMeans(grp_mean), decreasing = TRUE))
top_gen  <- head(setdiff(top_gen, "Unclassified"), 12)
comp <- rbind(grp_mean[top_gen, , drop = FALSE], `Other / unclassified` = 1 - colSums(grp_mean[top_gen, , drop = FALSE]))
comp_df <- data.frame(group = factor(rep(colnames(comp), each = nrow(comp)), levels = grp_levels),
                      genus = factor(rep(rownames(comp), times = ncol(comp)), levels = rownames(comp)),
                      rel   = as.vector(comp))
p_f <- ggplot(comp_df, aes(x = group, y = 100 * rel, fill = genus)) +
  geom_col() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(title = "Genus composition: kit controls, water and skin (before decontamination)",
       x = NULL, y = "Relative abundance (%)", fill = "Genus")

# (g) share of each skin sample's reads in flagged ASVs and in ASVs seen in the controls
share_df <- rbind(
  data.frame(sample = sample_names(ps)[is_skin], seq_batch = seq_batch[is_skin],
             measure = "Reads in flagged contaminant ASVs (%)",
             pct = 100 * colSums(m0[is_contam, is_skin, drop = FALSE]) / pmax(1, colSums(m0[, is_skin, drop = FALSE]))),
  data.frame(sample = sample_names(ps)[is_skin], seq_batch = seq_batch[is_skin],
             measure = "Reads in ASVs detected in kit controls (%)",
             pct = 100 * colSums(m0[kit_asv, is_skin, drop = FALSE]) / pmax(1, colSums(m0[, is_skin, drop = FALSE])))
)
p_g <- ggplot(share_df, aes(x = seq_batch, y = pct)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.2) +
  facet_wrap(~measure, scales = "free_y") +
  theme_minimal() +
  labs(title = "Share of each skin sample's reads linked to the negative controls",
       x = "Sequencing run", y = "% of the sample's reads")

pdf(file.path(QC_DIR, "contaminants_vs_skin.pdf"), width = 9, height = 6)
for (p in list(p_a, p_b, p_c, p_d, p_e, p_f, p_g)) print(p)
dev.off()

###############################################################################
# GRAPHICS 2: sequencing-batch checks on the final skin dataset
###############################################################################
sdr <- as(sample_data(Skin_rare), "data.frame")
mr  <- otu_mat(Skin_rare)
batch_lines <- c("Skin samples per sequencing run (rarefied dataset):", capture.output(print(table(sdr$seq_batch))))

if (length(unique(sdr$seq_batch)) > 1) {
  # (a) ordination coloured by run
  dr  <- vegdist(t(mr), method = "bray")
  pcr <- cmdscale(dr, k = 2, eig = TRUE)
  ver <- round(100 * pcr$eig[1:2] / sum(pcr$eig[pcr$eig > 0]), 1)
  pcr_df <- data.frame(PC1 = pcr$points[, 1], PC2 = pcr$points[, 2], seq_batch = sdr$seq_batch)
  p_b1 <- ggplot(pcr_df, aes(x = PC1, y = PC2, color = seq_batch)) +
    geom_point(size = 2.5, alpha = 0.8) +
    theme_minimal() +
    labs(title = "PCoA of rarefied skin samples by sequencing run (Bray-Curtis)",
         x = paste0("PCoA1 (", ver[1], "%)"), y = paste0("PCoA2 (", ver[2], "%)"), color = "Sequencing run")

  # (b) prevalence of each final taxon in run 1 against run 2
  r2 <- sdr$seq_batch == "run2"
  prev_df <- data.frame(ASV = rownames(mr),
                        prev_run1_pct = round(100 * rowMeans(mr[, !r2, drop = FALSE] > 0), 1),
                        prev_run2_pct = round(100 * rowMeans(mr[, r2, drop = FALSE] > 0), 1),
                        mean_reads_run1 = round(rowMeans(mr[, !r2, drop = FALSE]), 1),
                        mean_reads_run2 = round(rowMeans(mr[, r2, drop = FALSE]), 1),
                        tax_mat(Skin_rare)[rownames(mr), , drop = FALSE], check.names = FALSE)
  prev_df <- prev_df[order(-abs(prev_df$prev_run1_pct - prev_df$prev_run2_pct)), ]
  write_tsv(prev_df, "taxon_prevalence_by_run.tsv")
  p_b2 <- ggplot(prev_df, aes(x = prev_run1_pct, y = prev_run2_pct)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_point(alpha = 0.6) +
    theme_minimal() +
    labs(title = "Prevalence of each final skin taxon in run 1 vs run 2",
         subtitle = "Points far from the dashed line are taxa whose detection differs between runs",
         x = "Run 1 skin samples with the taxon (%)", y = "Run 2 skin samples with the taxon (%)")

  # (c) depth before rarefaction by run
  p_b3 <- ggplot(data.frame(depth = skin_depth, seq_batch = seq_batch[names(skin_depth)]),
                 aes(x = seq_batch, y = depth)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.5) +
    geom_hline(yintercept = RAREFY_DEPTH, linetype = "dashed") +
    theme_minimal() +
    labs(title = "Skin sample depth after filtering, by sequencing run", x = "Sequencing run", y = "Reads")

  pdf(file.path(QC_DIR, "batch_checks.pdf"), width = 8, height = 6)
  for (p in list(p_b1, p_b2, p_b3)) print(p)
  dev.off()

  # Does community composition differ between runs?
  set.seed(123)
  perm <- adonis2(dr ~ seq_batch, data = sdr, permutations = 999)
  disp <- anova(betadisper(dr, sdr$seq_batch))
  batch_lines <- c(batch_lines, "", "PERMANOVA (Bray-Curtis ~ sequencing run, 999 permutations):",
                   capture.output(print(perm)), "",
                   "Homogeneity of dispersion between runs (betadisper):", capture.output(print(disp)))
}

# Lice data of the run-2 fish, to check them against the burden classes later
b2 <- depth_tab[depth_tab$seq_batch == "run2" & depth_tab$sample_type == "skin",
                c("sample", "after_prevalence_filter", "kept_after_rarefaction")]
b2$Total_caligus <- to_num(sd0[b2$sample, "Total_caligus"])
if ("Final_Length" %in% colnames(sd0)) {
  lice_cm <- to_num(sd0$Total_caligus) / to_num(sd0$Final_Length)
  names(lice_cm) <- rownames(sd0)
  b2$lice_per_cm <- round(lice_cm[b2$sample], 3)
  b2$lice_per_cm_percentile <- round(100 * ecdf(lice_cm[is_skin])(lice_cm[b2$sample]), 0)
}
write_tsv(b2, "run2_samples.tsv")
writeLines(c(batch_lines, "", "Run-2 skin samples:", capture.output(print(b2, row.names = FALSE))),
           file.path(QC_DIR, "batch_checks.txt"))
message(paste(batch_lines, collapse = "\n"))

###############################################################################
# RAREFACTION PROCEDURE AND NORMALIZATION ALTERNATIVES (skin dataset)
# Every alternative starts from the same filtered, unrarefied skin counts.
# They are compared on technical grounds only: samples kept, how much of the
# community variation is explained by sequencing depth and by sequencing run,
# and agreement with the single rarefaction used so far.
###############################################################################
counts <- t(otu_mat(Skin_ps))                       # samples x taxa
counts <- counts[, colSums(counts) > 0, drop = FALSE]
lib    <- rowSums(counts)                           # depth before normalization

## ---- Repeated rarefaction ----------------------------------------------------
repeated_rarefy <- function(mat, depth, n_iter, seed) {
  sub <- mat[rowSums(mat) >= depth, , drop = FALSE]
  set.seed(seed)
  acc_tab  <- matrix(0, nrow(sub), ncol(sub), dimnames = dimnames(sub))
  acc_dist <- matrix(0, nrow(sub), nrow(sub), dimnames = list(rownames(sub), rownames(sub)))
  rich <- matrix(NA_real_, nrow(sub), n_iter, dimnames = list(rownames(sub), NULL))
  shan <- rich
  for (k in seq_len(n_iter)) {
    r <- vegan::rrarefy(sub, depth)
    acc_tab   <- acc_tab + r
    acc_dist  <- acc_dist + as.matrix(vegdist(r, method = "bray"))
    rich[, k] <- rowSums(r > 0)
    shan[, k] <- vegan::diversity(r, index = "shannon")
  }
  list(mean  = acc_tab / n_iter,
       dist  = as.dist(acc_dist / n_iter),
       alpha = data.frame(sample = rownames(sub),
                          richness_mean = rowMeans(rich), richness_sd = apply(rich, 1, sd),
                          shannon_mean  = rowMeans(shan), shannon_sd  = apply(shan, 1, sd)))
}

message("Repeated rarefaction: ", N_RAREFY_ITER, " draws at ", RAREFY_DEPTH, " reads, and at the minimum skin depth (", min(lib), ")")
rr_main <- repeated_rarefy(counts, RAREFY_DEPTH, N_RAREFY_ITER, RAREFY_SEED)   # same samples as the single draw
rr_all  <- repeated_rarefy(counts, min(lib),     N_RAREFY_ITER, RAREFY_SEED)   # keeps every skin sample

# How much does one random draw matter?
single      <- t(otu_mat(Skin_rare))
single_dist <- vegdist(single, method = "bray")
stab <- rr_main$alpha
stab$richness_single_draw <- rowSums(single > 0)[stab$sample]
stab$shannon_single_draw  <- vegan::diversity(single, index = "shannon")[stab$sample]
write_tsv(stab, "rarefaction_stability_per_sample.tsv")

sm_common <- intersect(attr(single_dist, "Labels"), attr(rr_main$dist, "Labels"))
d_single  <- as.matrix(single_dist)[sm_common, sm_common]
d_mean    <- as.matrix(rr_main$dist)[sm_common, sm_common]
stab_lines <- c(
  paste0("Repeated rarefaction: ", N_RAREFY_ITER, " draws at ", RAREFY_DEPTH, " reads, ", nrow(stab), " samples"),
  paste0("Richness: median SD across draws = ", round(median(stab$richness_sd), 2),
         " ASVs (median richness ", round(median(stab$richness_mean), 1), ")"),
  paste0("Shannon: median SD across draws = ", signif(median(stab$shannon_sd), 3),
         " (median Shannon ", round(median(stab$shannon_mean), 3), ")"),
  paste0("Largest Shannon SD in any sample = ", signif(max(stab$shannon_sd), 3)),
  paste0("Bray-Curtis, single draw vs mean of all draws: Spearman rho = ",
         round(cor(d_single[lower.tri(d_single)], d_mean[lower.tri(d_mean)], method = "spearman"), 4))
)
writeLines(stab_lines, file.path(QC_DIR, "rarefaction_stability_summary.txt"))
message(paste(stab_lines, collapse = "\n"))

## ---- Normalizations without rarefaction (all skin samples kept) --------------
tss  <- counts / lib                                         # proportions
hell <- sqrt(tss)                                            # Hellinger
lc   <- log(counts + CLR_PSEUDOCOUNT)
clr  <- lc - rowMeans(lc)                                    # centred log-ratio
lr   <- log(counts); lr[counts == 0] <- NA
rclr <- lr - rowMeans(lr, na.rm = TRUE); rclr[is.na(rclr)] <- 0   # robust CLR: zeros are left out, no pseudocount

norm <- list(
  rarefied_single     = list(X = single,       dist = single_dist,                   distance = "Bray-Curtis"),
  rarefied_mean       = list(X = rr_main$mean, dist = rr_main$dist,                  distance = "Bray-Curtis (mean of draws)"),
  rarefied_mean_all   = list(X = rr_all$mean,  dist = rr_all$dist,                   distance = "Bray-Curtis (mean of draws)"),
  proportions         = list(X = tss,          dist = vegdist(tss, method = "bray"), distance = "Bray-Curtis"),
  hellinger           = list(X = hell,         dist = dist(hell),                    distance = "Euclidean (Hellinger)"),
  clr                 = list(X = clr,          dist = dist(clr),                     distance = "Euclidean (Aitchison)"),
  rclr                = list(X = rclr,         dist = dist(rclr),                    distance = "Euclidean (robust Aitchison)")
)

# Optional: cumulative sum scaling, only if metagenomeSeq is installed
if (requireNamespace("metagenomeSeq", quietly = TRUE)) {
  css <- tryCatch({
    mr <- metagenomeSeq::newMRexperiment(t(counts))
    pq <- tryCatch(metagenomeSeq::cumNormStatFast(mr), error = function(e) metagenomeSeq::cumNormStat(mr))
    mr <- metagenomeSeq::cumNorm(mr, p = pq)
    t(metagenomeSeq::MRcounts(mr, norm = TRUE, log = TRUE))
  }, error = function(e) { warning("CSS normalization failed: ", conditionMessage(e)); NULL })
  if (!is.null(css)) norm$css <- list(X = css, dist = vegdist(css, method = "bray"), distance = "Bray-Curtis")
} else {
  message("metagenomeSeq not installed: CSS normalization skipped.")
}

## ---- Compare the alternatives ------------------------------------------------
eval_norm <- function(name, obj, ref_dist) {
  d  <- obj$dist
  sm <- attr(d, "Labels")
  meta <- data.frame(log_depth = log10(lib[sm]), seq_batch = factor(seq_batch[sm]), row.names = sm)
  has_batch <- nlevels(meta$seq_batch) > 1

  # Marginal effects: each term tested after accounting for the other
  set.seed(PERM_SEED)
  fm   <- if (has_batch) d ~ log_depth + seq_batch else d ~ log_depth
  perm <- adonis2(fm, data = meta, permutations = N_PERM, by = "margin")
  disp_p <- if (has_batch) anova(betadisper(d, meta$seq_batch))[1, "Pr(>F)"] else NA

  # Agreement of the between-sample distances with the single rarefaction
  common <- intersect(sm, attr(ref_dist, "Labels"))
  a <- as.matrix(d)[common, common]
  b <- as.matrix(ref_dist)[common, common]

  pc <- cmdscale(d, k = 2, eig = TRUE)
  list(
    row = data.frame(
      method = name, distance = obj$distance, n_samples = length(sm), n_taxa = ncol(obj$X),
      depth_R2 = round(perm["log_depth", "R2"], 4), depth_p = perm["log_depth", "Pr(>F)"],
      run_R2   = if (has_batch) round(perm["seq_batch", "R2"], 4) else NA,
      run_p    = if (has_batch) perm["seq_batch", "Pr(>F)"] else NA,
      run_dispersion_p = round(disp_p, 4),
      agreement_with_single_rarefaction = round(cor(a[lower.tri(a)], b[lower.tri(b)], method = "spearman"), 3)),
    pcoa = data.frame(method = name, sample = sm, PC1 = pc$points[, 1], PC2 = pc$points[, 2],
                      log_depth = meta$log_depth, seq_batch = meta$seq_batch)
  )
}

res <- lapply(names(norm), function(nm) eval_norm(nm, norm[[nm]], single_dist))
norm_tab <- do.call(rbind, lapply(res, `[[`, "row"))
pc_all   <- do.call(rbind, lapply(res, `[[`, "pcoa"))
pc_all$method <- factor(pc_all$method, levels = names(norm))
write_tsv(norm_tab, "normalization_comparison.tsv")
message("Normalization comparison (lower depth_R2 and run_R2 are better):")
print(norm_tab, row.names = FALSE)

r2_long <- rbind(data.frame(method = norm_tab$method, term = "Sequencing depth", R2 = norm_tab$depth_R2),
                 data.frame(method = norm_tab$method, term = "Sequencing run",   R2 = norm_tab$run_R2))
r2_long$method <- factor(r2_long$method, levels = names(norm))
p_n1 <- ggplot(r2_long, aes(x = method, y = 100 * R2, fill = term)) +
  geom_col(position = "dodge") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(title = "Community variation explained by technical factors, by normalization",
       subtitle = "Marginal PERMANOVA R2; lower is better",
       x = NULL, y = "Variation explained (%)", fill = NULL)
p_n2 <- ggplot(pc_all, aes(x = PC1, y = PC2, color = seq_batch)) +
  geom_point(size = 1.6, alpha = 0.8) +
  facet_wrap(~method, scales = "free") +
  theme_minimal() +
  labs(title = "PCoA of skin samples under each normalization, by sequencing run", color = "Sequencing run")
p_n3 <- ggplot(pc_all, aes(x = PC1, y = PC2, color = log_depth)) +
  geom_point(size = 1.6, alpha = 0.8) +
  facet_wrap(~method, scales = "free") +
  scale_color_viridis_c() +
  theme_minimal() +
  labs(title = "PCoA of skin samples under each normalization, by sequencing depth", color = "log10 reads")
p_n4 <- ggplot(stab, aes(x = shannon_single_draw, y = shannon_mean)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_errorbar(aes(ymin = shannon_mean - shannon_sd, ymax = shannon_mean + shannon_sd), width = 0, alpha = 0.5) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(title = "Shannon diversity: single rarefaction vs repeated rarefaction",
       subtitle = paste0("Points: mean of ", N_RAREFY_ITER, " draws; bars: +/- 1 SD across draws"),
       x = "Single draw (seed used so far)", y = "Mean of repeated draws")

pdf(file.path(QC_DIR, "normalization_comparison.pdf"), width = 10, height = 7)
for (p in list(p_n1, p_n2, p_n3, p_n4)) print(p)
dev.off()

## ---- Save each alternative as a phyloseq object ------------------------------
# Same samples, taxonomy, tree and sequences as the skin object; only the
# abundance table differs. Skin_rare_SILVA.rds remains the single rarefaction.
make_ps <- function(X) {
  out <- prune_taxa(colnames(X), prune_samples(rownames(X), Skin_ps))
  otu_table(out) <- otu_table(t(X), taxa_are_rows = TRUE)
  out
}
for (nm in setdiff(names(norm), "rarefied_single")) {
  saveRDS(make_ps(norm[[nm]]$X), file.path(RDS_DIR, paste0("Skin_", nm, "_SILVA.rds")))
}

## --------------------------- Persist QC outputs ------------------------------
write_tsv(stages, "stage_summary.tsv")
print(stages)

saveRDS(ps.noncontam,    file.path(RDS_DIR, "ps_noncontam.rds"))
saveRDS(Filtered_ps2,    file.path(RDS_DIR, "Filtered_ps2.rds"))
saveRDS(bacteria_all,    file.path(RDS_DIR, "bacteria_unfiltered_physeq.rds"))  # before the prevalence filter
saveRDS(bacteria_physeq, file.path(RDS_DIR, "bacteria_physeq.rds"))             # skin + water, taxa kept by either filter
saveRDS(Skin_ps,         file.path(RDS_DIR, "Skin_ps_SILVA.rds"))               # skin, skin filter, not rarefied
saveRDS(Ps_water,        file.path(RDS_DIR, "Water_ps_SILVA.rds"))              # water, water filter
saveRDS(Skin_rare,       file.path(RDS_DIR, "Skin_rare_SILVA.rds"))             # skin, rarefied
saveRDS(merged_physeq,   file.path(RDS_DIR, "merged_physeq_SILVA.rds"))         # water + rarefied skin

writeLines(c(
  paste("decontam: prevalence method, threshold", DECONTAM_THRESHOLD, "; negatives =", sum(is_kit), "kit controls"),
  paste("contaminant ASVs removed:", sum(contamdf$contaminant)),
  paste("sequencing run 2 samples:", paste(BATCH2_SAMPLES, collapse = ", ")),
  paste("tree rooting:", root_note),
  paste("organelle fix applied:", FIX_ORGANELLE_FILTER),
  paste("skin prevalence filter: >", SKIN_MIN_READS, "reads in >", SKIN_MIN_PREV, "of", length(skin_names), "skin samples;", sum(pass_skin), "taxa kept"),
  paste("water prevalence filter: >", WATER_MIN_READS, "reads in >", WATER_MIN_PREV, "of", length(water_names), "water samples;", sum(pass_water), "taxa kept"),
  paste("rarefaction: depth", RAREFY_DEPTH, ", seed", RAREFY_SEED, ", skin samples removed:", length(lost)),
  paste("repeated rarefaction:", N_RAREFY_ITER, "draws; CLR pseudocount:", CLR_PSEUDOCOUNT),
  paste("normalizations saved:", paste(setdiff(names(norm), "rarefied_single"), collapse = ", ")),
  capture.output(sessionInfo())
), file.path(QC_DIR, "qc_run_record.txt"))

message("Decontamination & QC complete. RDS files in ", RDS_DIR, "; tables and plots in ", QC_DIR)
