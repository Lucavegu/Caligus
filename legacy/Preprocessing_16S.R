#!/usr/bin/env Rscript

###############################################################################
# Skin microbiome & Caligus: Preprocessing pipeline for amplicon sequencing
# based on the 16S rRNA gene (DADA2 → taxonomy → phyloseq)
# Collaborative PhD Project between University of Chile, Laval University, and
# Benchmark Genetics
# Author: Lucas Venegas (adapted & curated)
# Email: lucas.venegas@ug.uchile.cl
# Date: 2025-08-09
# R version: >=4.2; Packages: dada2 (>=1.26), ShortRead, Biostrings, phyloseq,
#            DECIPHER, phangorn, dplyr
###############################################################################

## --------------------------- Reproducibility --------------------------------
set.seed(777)
options(stringsAsFactors = FALSE)

# Print session info at the end (also saved to outputs/logs)
.on.exit({
  si <- capture.output(sessionInfo())
  writeLines(si, file.path(OUTPUT_DIR, "logs", "sessionInfo.txt"))
})

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

# Cutadapt binary location (ensure it is installed and on your system)
CUTADAPT_BIN <- "/usr/bin/cutadapt"       # change if needed

# SILVA reference files (v138.2) — update paths if different
SILVA_TRAIN <- "~/Databases/silva_nr99_v138.2_toGenus_trainset.fa.gz"
SILVA_SPEC  <- "~/Databases/silva_v138.2_assignSpecies.fa.gz"

# Metadata file (tab-delimited; rownames = sample IDs)
METADATA_TSV <- file.path(PROJECT_DIR, "Metadata_microbiota_caligus.txt")

# Output subfolders
DIRS <- list(
  filtN    = file.path(OUTPUT_DIR, "filtN"),
  cutadapt = file.path(OUTPUT_DIR, "cutadapt"),
  filtered = file.path(OUTPUT_DIR, "cutadapt", "filtered"),
  fasta    = file.path(OUTPUT_DIR),
  tables   = file.path(OUTPUT_DIR),
  rds      = file.path(OUTPUT_DIR),
  logs     = file.path(OUTPUT_DIR, "logs")
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
})

message("dada2 version:", as.character(packageVersion("dada2")))
message("ShortRead version:", as.character(packageVersion("ShortRead")))
message("Biostrings version:", as.character(packageVersion("Biostrings")))

## ----------------------------- Input FASTQs ----------------------------------
# Expected file name pattern: SAMPLENAME_R1_001.fastq(.gz) and SAMPLENAME_R2_001.fastq(.gz)
fnFs <- sort(list.files(RAW_DIR, pattern = "_R1_001\\.fastq(\\.gz)?$", full.names = TRUE))
fnRs <- sort(list.files(RAW_DIR, pattern = "_R2_001\\.fastq(\\.gz)?$", full.names = TRUE))
stopifnot(length(fnFs) == length(fnRs), length(fnFs) > 0)

# --- Sample names are everything before "_V" ---
sample.names <- sapply(strsplit(basename(fnFs), "_[V]"), `[`, 1)
# Sanity: forward/reverse derive the same names
sample.names.R2 <- sapply(strsplit(basename(fnRs), "_[V]"), `[`, 1)
stopifnot(identical(sample.names, sample.names.R2))

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

## ------------------------------- Cutadapt ------------------------------------
# Verify cutadapt is available
ca_ok <- tryCatch({
  suppressWarnings(system2(CUTADAPT_BIN, args = "--version", stdout = TRUE, stderr = TRUE))
  TRUE
}, error = function(e) FALSE)
if (!ca_ok) stop("Cutadapt not found at ", CUTADAPT_BIN)

path.cut <- DIRS$cutadapt
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD_PRIMER)
REV.RC <- dada2:::rc(REV_PRIMER)
R1.flags <- paste("-g", FWD_PRIMER, "-a", REV.RC)  # trim FWD and rev-comp of REV from R1
R2.flags <- paste("-G", REV_PRIMER, "-A", FWD.RC)   # trim REV and rev-comp of FWD from R2

for (i in seq_along(fnFs)) {
  system2(CUTADAPT_BIN,
          args = c(R1.flags, R2.flags, "-n", 2,
                   "-o", fnFs.cut[i], "-p", fnRs.cut[i], fnFs.filtN[i], fnRs.filtN[i]))
}

## Optional sanity check on an arbitrary sample (change index if desired)
if (length(fnFs.cut) >= 1) {
  idx <- min(50L, length(fnFs.cut))
  chk <- rbind(
    FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[idx]]),
    FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[idx]]),
    REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[idx]]),
    REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[idx]])
  )
  write.table(chk, file.path(DIRS$logs, "cutadapt_primer_check.tsv"), sep = "\t", quote = FALSE)
}

## ----------------------------- Quality filtering -----------------------------
cutFs <- sort(list.files(path.cut, pattern = "R1_001\\.fastq(\\.gz)?$", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "R2_001\\.fastq(\\.gz)?$", full.names = TRUE))

# Sanity: after cutadapt, names (pre-"_V") still match
sample.names.cut.R1 <- sapply(strsplit(basename(cutFs), "_[V]"), `[`, 1)
sample.names.cut.R2 <- sapply(strsplit(basename(cutRs), "_[V]"), `[`, 1)
stopifnot(identical(sample.names, sample.names.cut.R1),
          identical(sample.names.cut.R1, sample.names.cut.R2))

filtFs <- file.path(DIRS$filtered, basename(cutFs))
filtRs <- file.path(DIRS$filtered, basename(cutRs))

# Reasonable defaults for V4; adjust if needed.
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs,
                     maxN = 0, maxEE = c(4, 4), truncQ = 2,
                     minLen = 50, rm.phix = TRUE, compress = TRUE,
                     multithread = THREADS > 1)
write.table(out, file.path(DIRS$logs, "filterAndTrim_counts.tsv"), sep = "\t", quote = FALSE)

## ----------------------------- Learn error models ----------------------------
errF <- learnErrors(filtFs, multithread = THREADS > 1)
#103895930 total bases in 448513 reads from 13 samples will be used for learning the error rates.
errR <- learnErrors(filtRs, multithread = THREADS > 1)
#103495304 total bases in 448513 reads from 13 samples will be used for learning the error rates.

# Optionally save error models
saveRDS(errF, file.path(DIRS$rds, "errF.rds"))
saveRDS(errR, file.path(DIRS$rds, "errR.rds"))

## --------------------------------- DADA --------------------------------------
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err = errF, multithread = THREADS, pool = FALSE)
dadaRs <- dada(derepRs, err = errR, multithread = THREADS, pool = FALSE)

## ------------------------------- Merge pairs ---------------------------------
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE, minOverlap = 12)

## --------------------------- Sequence table & QC -----------------------------
seqtab <- makeSequenceTable(mergers)
writeLines(paste("Sequence table dimensions:", paste(dim(seqtab), collapse = " x ")))
#Sequence table dimensions: 120 x 5909

len_dist <- table(nchar(getSequences(seqtab)))
write.table(as.data.frame(len_dist), file.path(DIRS$logs, "read_length_distribution.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                    multithread = THREADS > 1, verbose = TRUE)

# Chimera rate summary
chimera_rate <- sum(seqtab.nochim) / sum(seqtab)
writeLines(sprintf("Non-chimeric fraction retained: %.4f", chimera_rate))
#Non-chimeric fraction retained: 0.9960

## -------------------------- Tracking summary table ---------------------------
getN <- function(x) sum(getUniques(x))
summary_tab <- data.frame(
  sample = sample.names,
  dada2_input = out[, 1],
  filtered = out[, 2],
  dada_f = sapply(dadaFs, getN),
  dada_r = sapply(dadaRs, getN),
  merged = sapply(mergers, getN),
  nonchim = rowSums(seqtab.nochim),
  final_perc_reads_retained = round(rowSums(seqtab.nochim) / out[, 1] * 100, 1),
  row.names = sample.names
)

track_tsv <- file.path(DIRS$tables, "track_summary.tsv")
write.table(summary_tab, track_tsv, sep = "\t", quote = FALSE, col.names = NA)

## --------------------------------- Taxonomy ----------------------------------
if (!file.exists(path.expand(SILVA_TRAIN)) || !file.exists(path.expand(SILVA_SPEC))) {
  stop("SILVA reference files not found. Update SILVA_TRAIN and SILVA_SPEC paths.")
}

taxa <- assignTaxonomy(seqtab.nochim, path.expand(SILVA_TRAIN), multithread = THREADS > 1)
# Allow multiple species hits for higher recall at species level
taxa <- addSpecies(taxa, path.expand(SILVA_SPEC), allowMultiple = TRUE)

## ---------------------- Export ASV sequences & count table --------------------
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- paste0(">ASV_", seq_len(ncol(seqtab.nochim)))

# FASTA of ASV sequences
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, file.path(DIRS$fasta, "Caligus_microbiome_SILVA.fasta"))

# Count table (ASVs x samples)
asv_tab <- t(seqtab.nochim)
colnames(asv_tab) <- sample.names
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, file.path(DIRS$tables, "ASVs-counts_caligus_SILVA.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

# Taxonomy table (ASVs x ranks)
rownames(taxa) <- rownames(asv_tab)
write.table(taxa, file.path(DIRS$tables, "ASVs-taxonomy_caligus_SILVA.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

## --------------------------------- Phyloseq ----------------------------------
# Load metadata
metadata <- read.delim(METADATA_TSV, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

# Compute lice metrics helper
estimate_lice_metrics <- function(LC, BW) {
  stopifnot(length(LC) == length(BW))
  LogLC <- log(LC + 1)
  liceD <- LC / (BW^(2/3))
  LogLD <- log((LC + 1) / (BW^(2/3)))
  data.frame(LC = LC, BW = BW, LogLC = LogLC, liceD = liceD, LogLD = LogLD)
}

# Split and preprocess metadata (adapt indices to your design)
metadata_ctrl    <- metadata[1:7, , drop = FALSE]
metadata_caligus <- metadata[8:nrow(metadata), , drop = FALSE]

metadata_caligus$Total_caligus <- as.numeric(metadata_caligus$Total_caligus)
metadata_caligus$Final_Weight  <- as.numeric(metadata_caligus$Final_Weight)

res <- estimate_lice_metrics(metadata_caligus$Total_caligus, metadata_caligus$Final_Weight)
metadata_caligus$LogLC    <- res$LogLC
metadata_caligus$liceD    <- res$liceD
metadata_caligus$LogLiceD <- res$LogLD

# Clean-up temporary columns (if present)
metadata_caligus$LogLD <- NULL
metadata_ctrl$LogLD    <- NULL

# Zeroes for controls
metadata_ctrl$LogLC    <- 0
metadata_ctrl$liceD    <- 0
metadata_ctrl$LogLiceD <- 0

# Merge back & factor handling
metadata <- rbind(metadata_ctrl, metadata_caligus)
if ("Sex" %in% colnames(metadata)) metadata$Sex <- as.factor(metadata$Sex)

# Read tables for phyloseq
count_tab <- read.table(file.path(DIRS$tables, "ASVs-counts_caligus_SILVA.tsv"), header = TRUE, row.names = 1, check.names = FALSE, sep = "\t")
tax_tab   <- read.table(file.path(DIRS$tables, "ASVs-taxonomy_caligus_SILVA.tsv"), header = TRUE, row.names = 1, check.names = FALSE, sep = "\t")

# Build phyloseq object
count_tab_phy <- otu_table(as.matrix(count_tab), taxa_are_rows = TRUE)
# Ensure column order aligns with metadata row order
common_samples <- intersect(colnames(count_tab_phy), rownames(metadata))
count_tab_phy  <- count_tab_phy[, common_samples]
metadata       <- metadata[common_samples, , drop = FALSE]

# Additional sanity check: ensure names pre-"_V" match metadata rownames if they contain suffixes
stopifnot(all(colnames(count_tab_phy) %in% rownames(metadata)))

tax_tab_phy <- tax_table(as.matrix(tax_tab[rownames(count_tab_phy), , drop = FALSE]))

ps <- phyloseq(count_tab_phy, sample_data(metadata), tax_tab_phy)

# Add ASV sequences into refseq
DNA <- Biostrings::DNAStringSet(asv_seqs)
names(DNA) <- taxa_names(ps)
ps <- merge_phyloseq(ps, DNA)
taxa_names(ps) <- paste0("ASV_", seq(ntaxa(ps)))

saveRDS(ps, file.path(DIRS$rds, "phyloseq_caligus_microbiome_SILVA.rds"))

## --------------------------- Phylogenetic tree -------------------------------
# Multiple sequence alignment (DECIPHER)
alignment <- AlignSeqs(DNAStringSet(DNA), anchor = NA, verbose = TRUE)

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

message("Pipeline complete. Artifacts written to:")
message(" - Track summary: ", track_tsv)
message(" - ASV FASTA: ", file.path(DIRS$fasta, "Caligus_microbiome_SILVA.fasta"))
message(" - Count table: ", file.path(DIRS$tables, "ASVs-counts_caligus_SILVA.tsv"))
message(" - Taxonomy table: ", file.path(DIRS$tables, "ASVs-taxonomy_caligus_SILVA.tsv"))
message(" - Phyloseq RDS: ", file.path(DIRS$rds, "phyloseq_caligus_microbiome_SILVA.rds"))
message(" - Logs in: ", DIRS$logs)
