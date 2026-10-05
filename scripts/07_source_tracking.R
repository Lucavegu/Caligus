#!/usr/bin/env Rscript

###############################################################################
# Skin microbiome & Caligus: where do the skin sequences come from?
# Needs the outputs of the preprocessing and QC scripts:
#   phyloseq_caligus_microbiome_SILVA.rds   (before decontamination, with kits)
#   bacteria_unfiltered_physeq.rds          (after decontamination, no kits)
#
# The script has three parts:
#   1. Direct overlap: share of each skin sample's reads in ASVs that were also
#      detected in water and in the kit controls (no model, R only)
#   2. Writes the input files and the commands for SourceTracker2
#   3. Summarises the SourceTracker2 results, once they exist
# Run it, run the printed commands in a terminal, then run it again for part 3.
#
# Why it replaces the original SourceTracker step: the original used data from
# which the kit-like (and therefore water-like) ASVs had already been removed,
# and from which the prevalence filter had removed most water taxa. Skin was
# then compared with a water profile of a few hundred reads, so a result of
# ">99% Unknown" was close to inevitable.
###############################################################################

## --------------------------- User Parameters ---------------------------------
PROJECT_DIR <- getwd()
RDS_DIR     <- file.path(PROJECT_DIR)
ST_DIR      <- file.path(PROJECT_DIR, "sourcetracker")      # inputs, outputs, tables, plots

# How controls are labelled in the metadata
CTRL_PITTAG <- "CTRL"   # Pittag value of water samples
KIT_SEX     <- "Kit"    # Sex value of extraction-kit negative controls

# SourceTracker2
# The executable inside the conda environment (not the environment folder)
SOURCETRACKER_BIN <- "/Users/lucavegu/miniconda3/envs/st2/bin/sourcetracker2"
RUN_FROM_R   <- FALSE   # TRUE runs the SourceTracker2 commands from this script (can take a while)
SOURCE_DEPTH <- 5000    # reads drawn from each pooled source (lowered automatically if a source has fewer)
SINK_DEPTH   <- 5000    # reads drawn from each sink sample (sinks with fewer reads are left out and listed)
JOBS         <- 4      # parallel workers; set to 1 if the parallel start-up fails (slower, but needs no helper)

# "Detected in water" for the direct overlap, strict version:
# more than WATER_MIN_READS reads in at least WATER_MIN_SAMPLES water samples
WATER_MIN_READS   <- 3
WATER_MIN_SAMPLES <- 2

## ------------------------------ Libraries ------------------------------------
suppressPackageStartupMessages({
  library(phyloseq)
  library(biomformat)
  library(ggplot2)
})

## ------------------------------- Helpers -------------------------------------
otu_mat <- function(x) {            # counts with ASVs in rows
  m <- as(otu_table(x), "matrix")
  if (!taxa_are_rows(x)) m <- t(m)
  m
}
tax_mat <- function(x) as(tax_table(x), "matrix")
for (d in c(ST_DIR, file.path(ST_DIR, "inputs"), file.path(ST_DIR, "outputs"))) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}
write_tsv <- function(x, name) {
  write.table(x, file.path(ST_DIR, name), sep = "\t", quote = FALSE, row.names = FALSE)
}
pct_summary <- function(v) {
  c(median = median(v), q25 = unname(quantile(v, 0.25)), q75 = unname(quantile(v, 0.75)),
    min = min(v), max = max(v))
}

## -------------------------------- Input --------------------------------------
ps_raw   <- readRDS(file.path(RDS_DIR, "phyloseq_caligus_microbiome_SILVA.rds"))
ps_after <- readRDS(file.path(RDS_DIR, "bacteria_unfiltered_physeq.rds"))

sample_types <- function(ps) {
  sd <- as(sample_data(ps), "data.frame")
  is_kit   <- sd$Sex %in% KIT_SEX
  is_water <- sd$Pittag %in% CTRL_PITTAG & !is_kit
  st <- ifelse(is_kit, "Kit", ifelse(is_water, "Water", "Skin"))
  names(st) <- sample_names(ps)
  st
}

# Before decontamination: keep bacteria, remove chloroplasts and mitochondria.
# No contaminant removal and no prevalence filter.
tt <- tax_mat(ps_raw)
keep <- tt[, "Kingdom"] %in% "Bacteria" & !is.na(tt[, "Phylum"]) &
        !(tt[, "Order"] %in% "Chloroplast") & !(tt[, "Family"] %in% "Mitochondria")
m_before  <- otu_mat(ps_raw)[keep, , drop = FALSE]
st_before <- sample_types(ps_raw)

# After decontamination (as produced by the QC script, before its prevalence filter)
m_after  <- otu_mat(ps_after)
st_after <- sample_types(ps_after)

message("Before decontamination: ", nrow(m_before), " ASVs; samples: ",
        paste(names(table(st_before)), table(st_before), collapse = ", "))
message("After decontamination:  ", nrow(m_after), " ASVs; samples: ",
        paste(names(table(st_after)), table(st_after), collapse = ", "))

###############################################################################
# PART 1. Direct overlap between skin, water and kit controls
###############################################################################
skin  <- names(st_before)[st_before == "Skin"]
water <- names(st_before)[st_before == "Water"]
kit   <- names(st_before)[st_before == "Kit"]

in_water_any    <- rowSums(m_before[, water, drop = FALSE] > 0) > 0
in_water_strict <- rowSums(m_before[, water, drop = FALSE] > WATER_MIN_READS) >= WATER_MIN_SAMPLES
in_kit          <- rowSums(m_before[, kit, drop = FALSE] > 0) > 0
in_skin         <- rowSums(m_before[, skin, drop = FALSE] > 0) > 0

# Direction of the overlap. A shared ASV does not say who gave it to whom: fish
# shed skin bacteria into the water, and abundant skin ASVs can leak into other
# samples of the same sequencing run at low counts. ASVs that are relatively
# more abundant in water than in skin are the ones that could plausibly have
# come from the water.
rel_before     <- sweep(m_before, 2, pmax(1, colSums(m_before)), "/")
mean_rel_water <- rowMeans(rel_before[, water, drop = FALSE])
mean_rel_skin  <- rowMeans(rel_before[, skin, drop = FALSE])
mean_rel_kit   <- rowMeans(rel_before[, kit, drop = FALSE])
water_enriched   <- in_water_strict & mean_rel_water > mean_rel_skin
water_enriched10 <- in_water_strict & mean_rel_water > 10 * mean_rel_skin
skin_enriched    <- in_water_strict & !water_enriched

shared_tab <- data.frame(ASV = rownames(m_before),
                         mean_pct_water = round(100 * mean_rel_water, 4),
                         mean_pct_skin  = round(100 * mean_rel_skin, 4),
                         mean_pct_kit   = round(100 * mean_rel_kit, 4),
                         water_to_skin_ratio = signif(mean_rel_water / pmax(mean_rel_skin, 1e-9), 3),
                         in_kit = in_kit,
                         enriched_in = ifelse(water_enriched, "water", "skin"),
                         tt[keep, , drop = FALSE], check.names = FALSE)[in_water_strict, ]
write_tsv(shared_tab[order(-shared_tab$mean_pct_skin), ], "water_ASVs_abundance_in_water_skin_kit.tsv")

share <- function(asv_set, samples) {
  100 * colSums(m_before[asv_set, samples, drop = FALSE]) / pmax(1, colSums(m_before[, samples, drop = FALSE]))
}
overlap_skin <- data.frame(
  sample = skin,
  reads  = colSums(m_before[, skin, drop = FALSE]),
  pct_in_water_ASVs_any        = round(share(in_water_any, skin), 2),
  pct_in_water_ASVs_strict     = round(share(in_water_strict, skin), 2),
  pct_in_kit_ASVs              = round(share(in_kit, skin), 2),
  pct_in_water_only_ASVs       = round(share(in_water_strict & !in_kit, skin), 2),
  pct_in_water_or_kit_ASVs     = round(share(in_water_strict | in_kit, skin), 2),
  pct_in_skin_enriched_water_ASVs    = round(share(skin_enriched, skin), 2),
  pct_in_water_enriched_ASVs         = round(share(water_enriched, skin), 2),
  pct_in_water_enriched10_ASVs       = round(share(water_enriched10, skin), 2),
  pct_in_water_enriched_not_kit_ASVs = round(share(water_enriched & !in_kit, skin), 2)
)
write_tsv(overlap_skin, "direct_overlap_per_skin_sample.tsv")

measures <- c(pct_in_water_ASVs_any    = "ASVs detected in any water sample",
              pct_in_water_ASVs_strict = paste0("ASVs detected in water (>", WATER_MIN_READS, " reads in >=", WATER_MIN_SAMPLES, " samples)"),
              pct_in_kit_ASVs          = "ASVs detected in a kit control",
              pct_in_water_only_ASVs   = "ASVs detected in water (strict) but in no kit control",
              pct_in_water_or_kit_ASVs = "ASVs detected in water (strict) or a kit control",
              pct_in_skin_enriched_water_ASVs    = "Water (strict) ASVs that are more abundant in skin than in water",
              pct_in_water_enriched_ASVs         = "Water (strict) ASVs that are more abundant in water than in skin",
              pct_in_water_enriched10_ASVs       = "Water (strict) ASVs at least 10x more abundant in water than in skin",
              pct_in_water_enriched_not_kit_ASVs = "Water-enriched ASVs found in no kit control")
overlap_summary <- do.call(rbind, lapply(names(measures), function(k) {
  data.frame(share_of_skin_reads_in = measures[[k]], t(round(pct_summary(overlap_skin[[k]]), 2)))
}))
write_tsv(overlap_summary, "direct_overlap_summary.tsv")
message("Share of each skin sample's reads (%), before decontamination:")
print(overlap_summary, row.names = FALSE)

asv_sets <- data.frame(
  set = c("ASVs in skin", "ASVs in water (any)", "ASVs in water (strict)", "ASVs in kit controls",
          "Water (strict) ASVs more abundant in water than in skin",
          "Water (strict) ASVs at least 10x more abundant in water",
          "Water (strict) ASVs also in skin", "Water (strict) ASVs also in kits",
          "% of water reads in ASVs also found in skin", "% of water reads in ASVs also found in kits",
          "% of kit reads in ASVs also found in water (any)"),
  value = c(sum(in_skin), sum(in_water_any), sum(in_water_strict), sum(in_kit),
            sum(water_enriched), sum(water_enriched10),
            sum(in_water_strict & in_skin), sum(in_water_strict & in_kit),
            round(100 * sum(m_before[in_skin, water]) / sum(m_before[, water]), 1),
            round(100 * sum(m_before[in_kit, water]) / sum(m_before[, water]), 1),
            round(100 * sum(m_before[in_water_any, kit]) / sum(m_before[, kit]), 1))
)
write_tsv(asv_sets, "direct_overlap_asv_sets.tsv")
print(asv_sets, row.names = FALSE)

ov_long <- do.call(rbind, lapply(names(measures)[-1], function(k) {
  data.frame(measure = measures[[k]], pct = overlap_skin[[k]])
}))
ov_long$measure <- factor(ov_long$measure, levels = unname(measures[-1]))
p_ov <- ggplot(ov_long, aes(x = measure, y = pct)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Share of each skin sample's reads in ASVs shared with water or kit controls",
       subtitle = "Before decontamination; one point per fish",
       x = NULL, y = "% of the skin sample's reads")
ggsave(file.path(ST_DIR, "direct_overlap.pdf"), p_ov, width = 11, height = 6)

###############################################################################
# PART 2. SourceTracker2 inputs
###############################################################################
# Scenario A: skin as sink; Water and Kit as sources; before decontamination
# Scenario B: skin as sink; Water as the only source; after decontamination
# Scenario C: water samples as sinks; Kit as the only source; before decontamination
#             (how much of each water sample looks like the kit controls)
build_scenario <- function(name, m, env, role, auto_sink_depth = FALSE) {
  m <- m[, names(env), drop = FALSE]
  m <- m[rowSums(m) > 0, , drop = FALSE]
  totals  <- colSums(m)
  is_sink <- role == "sink"

  sink_depth <- if (auto_sink_depth) min(SINK_DEPTH, min(totals[is_sink])) else SINK_DEPTH
  dropped    <- names(totals)[is_sink & totals < sink_depth]
  keep_s     <- !(names(env) %in% dropped)
  m <- m[, keep_s, drop = FALSE]; env <- env[keep_s]; role <- role[keep_s]; totals <- totals[keep_s]

  pooled    <- tapply(totals[role == "source"], env[role == "source"], sum)
  src_depth <- min(SOURCE_DEPTH, min(pooled))
  if (src_depth < SOURCE_DEPTH) message("  ", name, ": source depth lowered to ", src_depth, " (smallest pooled source)")

  in_dir  <- file.path(ST_DIR, "inputs", name)
  out_dir <- file.path(ST_DIR, "outputs", name)
  if (!dir.exists(in_dir)) dir.create(in_dir, recursive = TRUE)
  biom_f <- file.path(in_dir, "table.biom")
  map_f  <- file.path(in_dir, "map.txt")
  write_biom(make_biom(data = m), biom_file = biom_f)
  map <- data.frame(`#SampleID` = colnames(m), SourceSink = role, Env = env, check.names = FALSE)
  write.table(map, map_f, sep = "\t", quote = FALSE, row.names = FALSE)

  cmd <- paste(SOURCETRACKER_BIN, "gibbs",
               "-i", shQuote(biom_f), "-m", shQuote(map_f), "-o", shQuote(out_dir),
               "--source_rarefaction_depth", src_depth,
               "--sink_rarefaction_depth", sink_depth,
               "--jobs", JOBS)
  message("  ", name, ": ", sum(role == "sink"), " sinks, sources: ",
          paste(names(pooled), collapse = " + "),
          if (length(dropped)) paste0("; sinks left out (below ", sink_depth, " reads): ", paste(dropped, collapse = ", ")) else "")
  list(name = name, out_dir = out_dir, cmd = cmd, dropped = dropped,
       source_depth = src_depth, sink_depth = sink_depth, sinks = colnames(m)[role == "sink"])
}

message("Writing SourceTracker2 inputs:")
envA  <- st_before
roleA <- ifelse(envA == "Skin", "sink", "source")
scA <- build_scenario("A_before_decontam_water_and_kit", m_before, envA, roleA)

envB  <- st_after[st_after %in% c("Skin", "Water")]
roleB <- ifelse(envB == "Skin", "sink", "source")
scB <- build_scenario("B_after_decontam_water_only", m_after, envB, roleB)

envC  <- st_before[st_before %in% c("Water", "Kit")]
roleC <- ifelse(envC == "Water", "sink", "source")
scC <- build_scenario("C_water_as_sink_kit_as_source", m_before, envC, roleC, auto_sink_depth = TRUE)

scenarios <- list(scA, scB, scC)

# Shell script with the three commands. SourceTracker2 refuses to write into an
# existing folder, so an earlier result is renamed first, not deleted.
# With --jobs above 1, SourceTracker2 starts a helper program (ipcluster) that
# lives in the same conda environment, so that folder is put on the PATH.
sh <- c("#!/bin/bash", "set -e",
        paste0("export PATH=", shQuote(dirname(path.expand(SOURCETRACKER_BIN))), ":\"$PATH\""))
for (sc in scenarios) {
  sh <- c(sh, "",
          paste0("if [ -d ", shQuote(sc$out_dir), " ]; then mv ", shQuote(sc$out_dir), " ",
                 shQuote(paste0(sc$out_dir, "_old_")), "$(date +%Y%m%d_%H%M%S); fi"),
          sc$cmd)
}
sh_file <- file.path(ST_DIR, "run_sourcetracker.sh")
writeLines(sh, sh_file)

if (dir.exists(SOURCETRACKER_BIN) || !file.exists(SOURCETRACKER_BIN)) {
  warning("SOURCETRACKER_BIN is not an executable file: ", SOURCETRACKER_BIN,
          "\nIt should end in /bin/sourcetracker2 inside the conda environment.", immediate. = TRUE)
} else if (RUN_FROM_R) {
  message("Running SourceTracker2 (output is saved in sourcetracker/run_log.txt)...")
  status <- system2("bash", shQuote(sh_file), stdout = file.path(ST_DIR, "run_log.txt"), stderr = file.path(ST_DIR, "run_log.txt"))
  if (status != 0) warning("SourceTracker2 stopped with an error; see ", file.path(ST_DIR, "run_log.txt"), immediate. = TRUE)
}

###############################################################################
# PART 3. Summarise SourceTracker2 results (only for scenarios that have run)
###############################################################################
read_mix <- function(f) {
  x <- as.matrix(read.delim(f, row.names = 1, check.names = FALSE))
  # Some versions write sources in rows and sinks in columns
  if ("Unknown" %in% rownames(x) && !("Unknown" %in% colnames(x))) x <- t(x)
  x
}

# Mean with a bootstrap 95% confidence interval. Source proportions pile up
# near 0 or 100, so a normal-theory interval (mean +/- 1.96 SE) is not suitable.
boot_ci <- function(v, n_boot = 2000) {
  if (length(v) < 2) return(c(NA_real_, NA_real_))
  b <- replicate(n_boot, mean(sample(v, replace = TRUE)))
  unname(quantile(b, c(0.025, 0.975)))
}

done <- FALSE
sd_raw <- as(sample_data(ps_raw), "data.frame")
all_summ   <- list()
plots      <- list()
test_lines <- c("Comparison of estimated sources within sinks.",
                "The same sinks are measured for every source, so paired tests are used",
                "(Friedman, then paired Wilcoxon with Holm correction), not Kruskal-Wallis/Dunn.", "")
set.seed(123)
for (sc in scenarios) {
  f <- file.path(sc$out_dir, "mixing_proportions.txt")
  if (!file.exists(f)) next
  done <- TRUE
  mix <- 100 * read_mix(f)
  src_cols <- colnames(mix)                       # the sources SourceTracker estimated
  if (all(c("Water", "Kit") %in% colnames(mix))) {
    mix <- cbind(mix, `Water + Kit` = mix[, "Water"] + mix[, "Kit"])
  }

  # Uncertainty of each estimate across the Gibbs draws, if SourceTracker wrote it
  f_sd <- file.path(sc$out_dir, "mixing_proportions_stds.txt")
  gibbs_sd <- if (file.exists(f_sd)) 100 * read_mix(f_sd) else NULL

  per_sample <- data.frame(sample = rownames(mix), round(mix, 3), check.names = FALSE)
  if ("seq_batch" %in% colnames(sd_raw)) per_sample$seq_batch <- sd_raw[per_sample$sample, "seq_batch"]
  write_tsv(per_sample, paste0("sourcetracker_", sc$name, "_per_sample.tsv"))

  summ <- do.call(rbind, lapply(colnames(mix), function(src) {
    v  <- mix[, src]
    ci <- boot_ci(v)
    data.frame(scenario = sc$name, source = src, n_sinks = length(v),
               mean = round(mean(v), 3), sd = round(sd(v), 3),
               mean_ci95_low = round(ci[1], 3), mean_ci95_high = round(ci[2], 3),
               t(round(pct_summary(v), 3)),
               n_sinks_above_1pct = sum(v > 1),
               median_gibbs_sd = if (!is.null(gibbs_sd) && src %in% colnames(gibbs_sd)) round(median(gibbs_sd[, src]), 3) else NA)
  }))
  all_summ[[sc$name]] <- summ

  long <- do.call(rbind, lapply(colnames(mix), function(src) data.frame(source = src, pct = mix[, src])))
  long$source <- factor(long$source, levels = colnames(mix))
  sub_txt <- paste0(nrow(mix), " sinks; source depth ", sc$source_depth, ", sink depth ", sc$sink_depth)
  plots[[length(plots) + 1]] <- ggplot(long, aes(x = source, y = pct)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
    theme_minimal() +
    labs(title = paste0("SourceTracker2, scenario ", sc$name), subtitle = sub_txt,
         x = "Estimated source", y = "% of the sink sample")
  plots[[length(plots) + 1]] <- ggplot(long, aes(x = pct, fill = source)) +
    geom_density(alpha = 0.5) +
    facet_wrap(~source, scales = "free") +
    theme_bw(base_size = 12) +
    theme(legend.position = "none") +
    labs(title = paste0("Density of source contributions, scenario ", sc$name), subtitle = sub_txt,
         x = "% of the sink sample", y = "Density")

  # Paired comparison of the estimated sources
  tl <- tryCatch(suppressWarnings({
    if (length(src_cols) >= 3 && nrow(mix) >= 3) {
      fr <- friedman.test(mix[, src_cols])
      pw <- pairwise.wilcox.test(as.vector(mix[, src_cols]), rep(src_cols, each = nrow(mix)),
                                 paired = TRUE, p.adjust.method = "holm")
      c(capture.output(print(fr)), capture.output(print(pw)))
    } else if (length(src_cols) == 2 && nrow(mix) >= 3) {
      capture.output(print(wilcox.test(mix[, src_cols[1]], mix[, src_cols[2]], paired = TRUE)))
    } else {
      "Too few sinks or sources for a test."
    }
  }), error = function(e) paste("Test not run:", conditionMessage(e)))
  test_lines <- c(test_lines, paste0("== Scenario ", sc$name, " (", nrow(mix), " sinks) =="), tl, "")
}

if (done) {
  summ_tab <- do.call(rbind, all_summ)
  write_tsv(summ_tab, "sourcetracker_summary.tsv")
  message("SourceTracker2 results (% of each sink sample):")
  print(summ_tab, row.names = FALSE)
  writeLines(test_lines, file.path(ST_DIR, "sourcetracker_tests.txt"))
  pdf(file.path(ST_DIR, "sourcetracker_results.pdf"), width = 8, height = 6)
  for (p in plots) print(p)
  dev.off()
}
missing_sc <- vapply(scenarios, function(sc) !file.exists(file.path(sc$out_dir, "mixing_proportions.txt")), logical(1))
if (any(missing_sc)) {
  message("\nSourceTracker2 results not found for: ",
          paste(vapply(scenarios[missing_sc], function(sc) sc$name, character(1)), collapse = ", "))
  message("Run this in a terminal where sourcetracker2 is available, then run this script again:\n  bash ",
          shQuote(sh_file))
}

writeLines(capture.output(sessionInfo()), file.path(ST_DIR, "sessionInfo.txt"))
message("Source tracking: tables and plots in ", ST_DIR)
