#!/usr/bin/env Rscript

###############################################################################
# Skin microbiome & Caligus: Preprocessing pipeline for amplicon sequencing
# based on the 16S rRNA gene (DADA2 -> taxonomy -> phyloseq)
# Collaborative PhD Project between University of Chile, Laval University, and
# Benchmark Genetics
# Author: Lucas Venegas (adapted & curated)
# Email: lucas.venegas@ug.uchile.cl
# Date: 2025-08-09 (revised 2026-10)
# R version: >=4.2; Packages: dada2 (>=1.26), ShortRead, Biostrings, phyloseq,
#            DECIPHER, phangorn, dplyr, ggplot2
###############################################################################

## --------------------------- Reproducibility --------------------------------
set.seed(777)
options(stringsAsFactors = FALSE)

## --------------------------- User Parameters ---------------------------------
# Root project directory (change if needed). All outputs created relative to this.
PROJECT_DIR <- getwd()                    # e.g., "/home/lucas/Caligus_microbiome/First_sequencing"
RAW_DIR     <- file.path(PROJECT_DIR)     # directory containing raw FASTQ(.gz) files
OUTPUT_DIR  <- file.path(PROJECT_DIR)     # write outputs here by default

# Number of threads for multithreaded steps (tune to your machine)
THREADS <- parallel::detectCores(logical = TRUE)

# Primer sequences (V4: 515F/806R)
FWD_PRIMER <- "GTGYCAGCMGCCGCGGTAA"       # 515F (Parada)
REV_PRIMER <- "GGACTACNVGGGTWTCTAAT"       # 806R (Apprill)

# Cutadapt executable (the program itself, not the conda environment folder)
CUTADAPT_BIN <- "/Users/lucavegu/miniconda3/envs/cutadapt/bin/cutadapt"

# SILVA reference files (v138.2)
SILVA_TRAIN <- "/Users/lucavegu/Desktop/Databases/silva_nr99_v138.2_toGenus_trainset.fa.gz"
SILVA_SPEC  <- "/Users/lucavegu/Desktop/Databases/silva_v138.2_assignSpecies.fa.gz"

# Metadata file (tab-delimited; rownames = sample IDs)
METADATA_TSV <- file.path(PROJECT_DIR, "Metadata_microbiota_caligus.txt")

# How controls are labelled in the metadata (same labels as QC&DECONTAM.R)
CTRL_PITTAG <- "CTRL"   # Pittag value of water controls
KIT_SEX     <- "Kit"    # Sex value of extraction-kit negative controls

# Samples sequenced in the second run; every other sample is "run1".
# DADA2 error models are learned, and samples denoised, separately for each run.
BATCH2_SAMPLES <- c("SG21_0185", "SG21_0186", "SG21_0190",
                    "SG21_0369", "SG21_0371", "SG21_0374")

# Output subfolders
DIRS <- list(
  filtN    = file.path(OUTPUT_DIR, "filtN"),
  cutadapt = file.path(OUTPUT_DIR, "cutadapt"),
  filtered = file.path(OUTPUT_DIR, "cutadapt", "filtered"),
  fasta    = file.path(OUTPUT_DIR),
  tables   = file.path(OUTPUT_DIR),
  rds      = file.path(OUTPUT_DIR),
  logs     = file.path(OUTPUT_DIR, "logs"),
  cutlogs  = file.path(OUTPUT_DIR, "logs", "cutadapt")
)

# Create output directories
invisible(lapply(DIRS, function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE)))

## ------------------------------ Libraries ------------------------------------
suppressPackageStartupMessages({
  library(dada2)
  library(ShortRead)
  library(Biostrings)
  library(phyloseq)
  library(DECIPHER)
  library(phangorn)
  library(dplyr)
  library(ggplot2)
})

message("dada2 version: ", as.character(packageVersion("dada2")))
message("ShortRead version: ", as.character(packageVersion("ShortRead")))
message("Biostrings version: ", as.character(packageVersion("Biostrings")))

## ----------------------------- Input FASTQs ----------------------------------
# Expected file name pattern: SAMPLENAME_R1_001.fastq(.gz) and SAMPLENAME_R2_001.fastq(.gz)
fnFs <- sort(list.files(RAW_DIR, pattern = "_R1_001\\.fastq(\\.gz)?$", full.names = TRUE))
fnRs <- sort(list.files(RAW_DIR, pattern = "_R2_001\\.fastq(\\.gz)?$", full.names = TRUE))
stopifnot(length(fnFs) == length(fnRs), length(fnFs) > 0)

# --- Sample names are everything before "_V" ---
get_names <- function(f) sapply(strsplit(basename(f), "_[V]"), `[`, 1)
sample.names <- get_names(fnFs)
# Sanity: forward/reverse derive the same names, and names are unique
stopifnot(identical(sample.names, get_names(fnRs)), !anyDuplicated(sample.names))
all_samples <- sample.names

## ------------------------------- Metadata ------------------------------------
# Loaded early so that problems are found before the long steps
if (!file.exists(METADATA_TSV)) stop("Metadata file not found: ", METADATA_TSV)
metadata <- read.delim(METADATA_TSV, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

needed_cols <- c("Pittag", "Sex", "Total_caligus", "Final_Weight")
if (!all(needed_cols %in% colnames(metadata))) {
  stop("Metadata is missing column(s): ",
       paste(setdiff(needed_cols, colnames(metadata)), collapse = ", "))
}

# Controls are identified by their labels, not by row position
flag_controls <- function(md) md$Pittag %in% CTRL_PITTAG | md$Sex %in% KIT_SEX
message("Controls found in metadata: ", sum(flag_controls(metadata)),
        " (kit: ", sum(metadata$Sex %in% KIT_SEX), ")")

# Compare FASTQ sample names with the metadata now, before the long steps
missing_meta <- setdiff(all_samples, rownames(metadata))
missing_seq  <- setdiff(rownames(metadata), all_samples)
message(length(all_samples), " FASTQ pairs; ", nrow(metadata), " metadata rows; ",
        length(missing_seq), " metadata rows without FASTQ files")
if (length(missing_meta) > 0) {
  stop("FASTQ samples with no metadata row: ", paste(missing_meta, collapse = ", "))
}
if (length(missing_seq) > 0) {
  writeLines(missing_seq, file.path(DIRS$logs, "metadata_rows_without_fastq.txt"))
  warning(length(missing_seq), " metadata rows have no FASTQ files and will be left out; ",
          "listed in logs/metadata_rows_without_fastq.txt", immediate. = TRUE)
}

## ----------------------------- Primer Checks ---------------------------------
allOrients <- function(primer) {
  dna <- DNAString(primer)
  orients <- c(Forward = dna,
               Complement = Biostrings::complement(dna),
               Reverse = Biostrings::reverse(dna),
               RevComp = Biostrings::reverseComplement(dna))
  sapply(orients, toString)
}
FWD.orients <- allOrients(FWD_PRIMER)
REV.orients <- allOrients(REV_PRIMER)

primerHits <- function(primer, fn) {
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  sum(nhits > 0)
}

## ----------Pre-filter (remove reads with ambiguous bases (Ns)-----------------
fnFs.filtN <- file.path(DIRS$filtN, basename(fnFs))
fnRs.filtN <- file.path(DIRS$filtN, basename(fnRs))

outN <- filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN,
                      maxN = 0, multithread = THREADS > 1)
rownames(outN) <- all_samples

# Samples with no reads left produce no output file; carry on without them
okN <- outN[, 2] > 0 & file.exists(fnFs.filtN) & file.exists(fnRs.filtN)

## ------------------------------- Cutadapt ------------------------------------
# Verify cutadapt is available
cutadapt_version <- tryCatch(
  suppressWarnings(system2(CUTADAPT_BIN, args = "--version", stdout = TRUE, stderr = TRUE)),
  error = function(e) NULL)
if (dir.exists(CUTADAPT_BIN) || !file.exists(CUTADAPT_BIN) || is.null(cutadapt_version) ||
    !is.null(attr(cutadapt_version, "status"))) {
  stop("CUTADAPT_BIN must be the cutadapt executable itself. Current value: ", CUTADAPT_BIN)
}
message("cutadapt version: ", cutadapt_version[1])

path.cut <- DIRS$cutadapt
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD_PRIMER)
REV.RC <- dada2:::rc(REV_PRIMER)
R1.flags <- paste("-g", FWD_PRIMER, "-a", REV.RC)  # trim FWD and rev-comp of REV from R1
R2.flags <- paste("-G", REV_PRIMER, "-A", FWD.RC)   # trim REV and rev-comp of FWD from R2

# Cutadapt output for each sample goes to logs/cutadapt/; the run stops if cutadapt fails
for (i in which(okN)) {
  log_i <- file.path(DIRS$cutlogs, paste0(all_samples[i], ".cutadapt.log"))
  status <- system2(CUTADAPT_BIN,
                    args = c(R1.flags, R2.flags, "-n", 2,
                             "-o", shQuote(fnFs.cut[i]), "-p", shQuote(fnRs.cut[i]),
                             shQuote(fnFs.filtN[i]), shQuote(fnRs.filtN[i])),
                    stdout = log_i, stderr = log_i)
  if (status != 0) stop("Cutadapt failed for sample ", all_samples[i], "; see ", log_i)
}

# Use the files written in this run (not whatever is in the folder)
cutFs <- fnFs.cut[okN]
cutRs <- fnRs.cut[okN]
names_cut <- all_samples[okN]
stopifnot(all(file.exists(cutFs)), all(file.exists(cutRs)))

## Optional sanity check on an arbitrary sample (change index if desired)
idx <- min(50L, length(cutFs))
chk <- rbind(
  FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = cutFs[[idx]]),
  FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = cutRs[[idx]]),
  REV.ForwardReads = sapply(REV.orients, primerHits, fn = cutFs[[idx]]),
  REV.ReverseReads = sapply(REV.orients, primerHits, fn = cutRs[[idx]])
)
write.table(chk, file.path(DIRS$logs, "cutadapt_primer_check.tsv"), sep = "\t", quote = FALSE)

## ----------------------------- Quality profiles ------------------------------
# For the supplementary material
# Raw reads are used, from up to 12 samples with at least 1,000 reads
tryCatch({
  qidx <- head(which(outN[, 1] >= 1000), 12)
  ggsave(file.path(DIRS$logs, "quality_profile_forward.pdf"),
         plotQualityProfile(fnFs[qidx], aggregate = TRUE), width = 7, height = 5)
  ggsave(file.path(DIRS$logs, "quality_profile_reverse.pdf"),
         plotQualityProfile(fnRs[qidx], aggregate = TRUE), width = 7, height = 5)
}, error = function(e) warning("Quality profile plots not saved: ", conditionMessage(e)))

## ----------------------------- Quality filtering -----------------------------
filtFs <- file.path(DIRS$filtered, basename(cutFs))
filtRs <- file.path(DIRS$filtered, basename(cutRs))

# Reasonable defaults for V4; adjust if needed.
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs,
                     maxN = 0, maxEE = c(4, 4), truncQ = 2,
                     minLen = 50, rm.phix = TRUE, compress = TRUE,
                     multithread = THREADS > 1)
rownames(out) <- names_cut
write.table(out, file.path(DIRS$logs, "filterAndTrim_counts.tsv"), sep = "\t", quote = FALSE, col.names = NA)

# Drop samples that lost all reads, and record them
okF <- out[, 2] > 0 & file.exists(filtFs) & file.exists(filtRs)
filtFs <- filtFs[okF]
filtRs <- filtRs[okF]
sample.names <- names_cut[okF]

dropped <- setdiff(all_samples, sample.names)
writeLines(dropped, file.path(DIRS$logs, "samples_dropped_no_reads.txt"))
if (length(dropped) > 0) {
  warning("Samples with no reads after filtering (dropped): ", paste(dropped, collapse = ", "))
}

## ------------- Error models, denoising and merging: one run at a time --------
# Error rates differ between sequencing runs, so each run gets its own error
# model. The merged reads of all runs are then combined into one sequence table.
not_found <- setdiff(BATCH2_SAMPLES, all_samples)
if (length(not_found) > 0) warning("BATCH2_SAMPLES not among the FASTQ files: ", paste(not_found, collapse = ", "))
seq_run <- ifelse(sample.names %in% BATCH2_SAMPLES, "run2", "run1")
names(seq_run) <- sample.names
message("Samples per sequencing run:")
print(table(seq_run))

getN <- function(x) sum(getUniques(x))
mergers  <- list()
n_dada_f <- c()
n_dada_r <- c()

for (run in sort(unique(seq_run))) {
  idx <- which(seq_run == run)
  if (length(idx) < 2) stop("Sequencing run '", run, "' has fewer than 2 samples.")
  message("=== ", run, ": ", length(idx), " samples ===")

  errF <- learnErrors(filtFs[idx], multithread = THREADS > 1)
  errR <- learnErrors(filtRs[idx], multithread = THREADS > 1)
  saveRDS(errF, file.path(DIRS$rds, paste0("errF_", run, ".rds")))
  saveRDS(errR, file.path(DIRS$rds, paste0("errR_", run, ".rds")))
  tryCatch({
    ggsave(file.path(DIRS$logs, paste0("error_model_forward_", run, ".pdf")), plotErrors(errF, nominalQ = TRUE), width = 8, height = 8)
    ggsave(file.path(DIRS$logs, paste0("error_model_reverse_", run, ".pdf")), plotErrors(errR, nominalQ = TRUE), width = 8, height = 8)
  }, error = function(e) warning("Error-model plots not saved: ", conditionMessage(e)))

  # One sample at a time, so only one sample's reads are in memory.
  # With pool = FALSE this gives the same result as denoising all samples together.
  for (i in idx) {
    sam <- sample.names[i]
    message("  ", sam)
    derepF <- derepFastq(filtFs[i])
    derepR <- derepFastq(filtRs[i])
    ddF <- dada(derepF, err = errF, multithread = THREADS, pool = FALSE, verbose = FALSE)
    ddR <- dada(derepR, err = errR, multithread = THREADS, pool = FALSE, verbose = FALSE)
    mergers[[sam]] <- mergePairs(ddF, derepF, ddR, derepR, minOverlap = 12)
    n_dada_f[sam]  <- getN(ddF)
    n_dada_r[sam]  <- getN(ddR)
    rm(derepF, derepR, ddF, ddR)
  }
  rm(errF, errR)
  invisible(gc())
}
mergers <- mergers[sample.names]

## --------------------------- Sequence table & QC -----------------------------
seqtab <- makeSequenceTable(mergers)
writeLines(paste("Sequence table dimensions:", paste(dim(seqtab), collapse = " x ")))

# Length distribution: number of ASVs and number of reads at each length
seq_len_bp <- nchar(getSequences(seqtab))
len_dist <- data.frame(
  length_bp = as.integer(names(table(seq_len_bp))),
  n_ASVs    = as.integer(table(seq_len_bp)),
  n_reads   = as.numeric(tapply(colSums(seqtab), seq_len_bp, sum))
)
write.table(len_dist, file.path(DIRS$logs, "read_length_distribution.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                    multithread = THREADS > 1, verbose = TRUE)

# Fraction of reads kept after chimera removal
nonchim_fraction <- sum(seqtab.nochim) / sum(seqtab)
writeLines(sprintf("Non-chimeric fraction retained: %.4f", nonchim_fraction))

# Checkpoint: the slow DADA2 steps are done. If R stops later, this file can be
# reloaded with readRDS() instead of repeating them.
saveRDS(seqtab.nochim, file.path(DIRS$rds, "seqtab_nochim.rds"))
rm(seqtab)
invisible(gc())

## -------------------------- Tracking summary table ---------------------------
# One row per input sample, from raw reads to final counts.
# Samples dropped along the way keep NA in the later columns.
getN <- function(x) sum(getUniques(x))
summary_tab <- data.frame(
  sample         = all_samples,
  raw            = outN[, 1],
  no_N           = outN[, 2],
  primer_trimmed = NA_real_,
  filtered       = NA_real_,
  dada_f         = NA_real_,
  dada_r         = NA_real_,
  merged         = NA_real_,
  nonchim        = NA_real_,
  row.names      = all_samples
)
summary_tab[names_cut, "primer_trimmed"] <- out[, 1]
summary_tab[names_cut, "filtered"]       <- out[, 2]
summary_tab[sample.names, "dada_f"]      <- n_dada_f[sample.names]
summary_tab[sample.names, "dada_r"]      <- n_dada_r[sample.names]
summary_tab[sample.names, "merged"]      <- sapply(mergers, getN)
summary_tab[rownames(seqtab.nochim), "nonchim"] <- rowSums(seqtab.nochim)
summary_tab$final_perc_of_raw <- round(summary_tab$nonchim / summary_tab$raw * 100, 1)
summary_tab$seq_run    <- ifelse(all_samples %in% BATCH2_SAMPLES, "run2", "run1")
summary_tab$is_control <- flag_controls(metadata[match(all_samples, rownames(metadata)), , drop = FALSE])

track_tsv <- file.path(DIRS$tables, "track_summary.tsv")
write.table(summary_tab, track_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

## --------------------------------- Taxonomy ----------------------------------
# Free memory: these large objects are no longer needed
rm(mergers)
invisible(gc())

if (!file.exists(path.expand(SILVA_TRAIN)) || !file.exists(path.expand(SILVA_SPEC))) {
  stop("SILVA reference files not found. Update SILVA_TRAIN and SILVA_SPEC paths.")
}

set.seed(777)  # taxonomy assignment uses random bootstraps
taxa <- assignTaxonomy(seqtab.nochim, path.expand(SILVA_TRAIN), multithread = THREADS > 1)
# Allow multiple species hits for higher recall at species level
# n = sequences matched at a time; a smaller n uses less memory and gives the same result
saveRDS(taxa, file.path(DIRS$rds, "taxa_genus.rds"))   # checkpoint before the memory-heavy step
invisible(gc())
taxa <- addSpecies(taxa, path.expand(SILVA_SPEC), allowMultiple = TRUE, n = 200)

## ---------------------- Export ASV sequences & count table --------------------
asv_seqs <- colnames(seqtab.nochim)
asv_ids  <- paste0("ASV_", seq_len(ncol(seqtab.nochim)))

# FASTA of ASV sequences
asv_fasta <- c(rbind(paste0(">", asv_ids), asv_seqs))
write(asv_fasta, file.path(DIRS$fasta, "Caligus_microbiome_SILVA.fasta"))

# Count table (ASVs x samples)
asv_tab <- t(seqtab.nochim)
rownames(asv_tab) <- asv_ids
write.table(asv_tab, file.path(DIRS$tables, "ASVs-counts_caligus_SILVA.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

# Taxonomy table (ASVs x ranks)
rownames(taxa) <- asv_ids
write.table(taxa, file.path(DIRS$tables, "ASVs-taxonomy_caligus_SILVA.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

## ---------------------------- Metadata processing ----------------------------
# Every sequenced sample must have metadata. Metadata rows without reads were
# reported at the start of the run and are left out here.
seq_samples  <- colnames(asv_tab)
missing_meta <- setdiff(seq_samples, rownames(metadata))
if (length(missing_meta) > 0) {
  stop("Sequenced samples with no metadata row: ", paste(missing_meta, collapse = ", "))
}
metadata <- metadata[seq_samples, , drop = FALSE]

# Compute lice metrics helper
estimate_lice_metrics <- function(LC, BW) {
  stopifnot(length(LC) == length(BW))
  LogLC <- log(LC + 1)
  liceD <- LC / (BW^(2/3))
  LogLD <- log((LC + 1) / (BW^(2/3)))
  data.frame(LC = LC, BW = BW, LogLC = LogLC, liceD = liceD, LogLD = LogLD)
}

is_ctrl <- flag_controls(metadata)

# Numeric phenotype columns; controls have no fish phenotype, so they become NA
to_num <- function(v) suppressWarnings(as.numeric(as.character(v)))
metadata$Total_caligus <- to_num(metadata$Total_caligus)
metadata$Final_Weight  <- to_num(metadata$Final_Weight)

bad_fish <- !is_ctrl & (is.na(metadata$Total_caligus) | is.na(metadata$Final_Weight))
if (any(bad_fish)) {
  warning("Fish with missing or non-numeric Total_caligus / Final_Weight: ",
          paste(rownames(metadata)[bad_fish], collapse = ", "))
}

# Lice metrics for fish only; NA (not 0) for water and kit controls
metadata$LogLC    <- NA_real_
metadata$liceD    <- NA_real_
metadata$LogLiceD <- NA_real_
res <- estimate_lice_metrics(metadata$Total_caligus[!is_ctrl], metadata$Final_Weight[!is_ctrl])
metadata$LogLC[!is_ctrl]    <- res$LogLC
metadata$liceD[!is_ctrl]    <- res$liceD
metadata$LogLiceD[!is_ctrl] <- res$LogLD

metadata$Sex <- as.factor(metadata$Sex)
metadata$seq_batch <- ifelse(rownames(metadata) %in% BATCH2_SAMPLES, "run2", "run1")

## --------------------------------- Phyloseq ----------------------------------
count_tab_phy <- otu_table(asv_tab, taxa_are_rows = TRUE)
tax_tab_phy   <- tax_table(as.matrix(taxa))

ps <- phyloseq(count_tab_phy, sample_data(metadata), tax_tab_phy)
stopifnot(nsamples(ps) == length(seq_samples))

# Add ASV sequences into refseq
DNA <- Biostrings::DNAStringSet(asv_seqs)
names(DNA) <- asv_ids
ps <- merge_phyloseq(ps, DNA)

## --------------------------- Phylogenetic tree -------------------------------
# Multiple sequence alignment (DECIPHER)
alignment <- AlignSeqs(DNA, anchor = NA, verbose = TRUE)

# Convert to phangorn format and compute distance matrix
phang.align <- phyDat(as(alignment, "matrix"), type = "DNA")
dm <- dist.ml(phang.align)

# Neighbor-Joining tree under JC69
nj_JC69 <- NJ(dm)
phy_tree(ps) <- nj_JC69

# Export ASV sequences as FASTA from refseq (for external tools)
Biostrings::writeXStringSet(refseq(ps), filepath = file.path(DIRS$fasta, "asv_caligus_microbiome_SILVA.fna"),
                            append = FALSE, compress = FALSE, compression_level = NA, format = "fasta")

# Save final phyloseq object
saveRDS(ps, file.path(DIRS$rds, "phyloseq_caligus_microbiome_SILVA.rds"))

# Session info
writeLines(capture.output(sessionInfo()), file.path(DIRS$logs, "sessionInfo.txt"))

message("Pipeline complete. Artifacts written to:")
message(" - Track summary: ", track_tsv)
message(" - ASV FASTA: ", file.path(DIRS$fasta, "Caligus_microbiome_SILVA.fasta"))
message(" - Count table: ", file.path(DIRS$tables, "ASVs-counts_caligus_SILVA.tsv"))
message(" - Taxonomy table: ", file.path(DIRS$tables, "ASVs-taxonomy_caligus_SILVA.tsv"))
message(" - Phyloseq RDS: ", file.path(DIRS$rds, "phyloseq_caligus_microbiome_SILVA.rds"))
message(" - Logs in: ", DIRS$logs)
