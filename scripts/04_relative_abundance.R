#!/usr/bin/env Rscript

###############################################################################
# Skin microbiome & Caligus: relative abundance of water and skin communities
# (Figure 1: A = water samples, B = skin samples)
# Replaces the "ABUNDANCE PLOT" block of Script_CMS.R. Needs the outputs of
# 02_qc_decontam.R.
#
# Differences from the original block:
#   - runs on its own and writes everything to the folder "abundance"
#   - uses the new objects: rarefied skin table and the separately filtered
#     water table
#   - colours are tied to taxon names, so a taxon has the same colour in both
#     panels whatever taxa are present
#   - ASVs without a name at the plotted rank are kept as "Unclassified"
#     (tax_glom dropped them, so the percentages did not refer to all reads)
#   - taxa outside the most abundant ones are shown as "Other", so every bar
#     sums to 1
#   - samples are ordered by the most abundant taxon of the skin panel
###############################################################################

## --------------------------- User Parameters ---------------------------------
PROJECT_DIR <- getwd()
RDS_DIR     <- file.path(PROJECT_DIR)
OUT_DIR     <- file.path(PROJECT_DIR, "abundance")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

SKIN_RDS  <- "Skin_rare_SILVA.rds"     # rarefied skin samples
WATER_RDS <- "Water_ps_SILVA.rds"      # water samples, water prevalence filter

TOP_N_PHYLUM <- 10     # most abundant phyla shown per panel
TOP_N_GENUS  <- 10     # most abundant genera shown per panel

# Colours of the original figure, by phylum name. Phyla not listed here get a
# colour from EXTRA_COLOURS.
PHYLUM_COLOURS <- c(Actinomycetota    = "#8F6F00",
                    Bacillota         = "#01843a",
                    Bacteroidota      = "#fed976",
                    Bdellovibrionota  = "#FFA73B",
                    Patescibacteria   = "#D7301F",
                    Pseudomonadota    = "#225EA8",
                    Verrucomicrobiota = "#35978F")
EXTRA_COLOURS <- c("#9467bd", "#e377c2", "#17becf", "#8c564b", "#bcbd22", "#1f77b4", "#ff7f0e", "#2ca02c",
                   "#d62728", "#7f7f7f", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5", "#c49c94",
                   "#f7b6d2", "#dbdb8d", "#9edae5", "#393b79")
GENUS_COLOURS <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#bcbd22",
                   "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5", "#c49c94", "#f7b6d2",
                   "#dbdb8d", "#9edae5", "#393b79", "#637939")

## ------------------------------ Libraries ------------------------------------
suppressPackageStartupMessages({
  library(phyloseq)
  library(ggplot2)
  library(patchwork)
})

## ------------------------------- Helpers -------------------------------------
otu_mat <- function(x) {            # counts with ASVs in rows
  m <- as(otu_table(x), "matrix")
  if (!taxa_are_rows(x)) m <- t(m)
  m
}
write_tsv <- function(x, name) {
  write.table(x, file.path(OUT_DIR, name), sep = "\t", quote = FALSE, row.names = FALSE)
}

# Relative abundance summed at one taxonomic rank (taxa in rows, samples in columns)
rel_by_rank <- function(ps, rank) {
  m   <- otu_mat(ps)
  rel <- sweep(m, 2, colSums(m), "/")
  lab <- as(tax_table(ps), "matrix")[rownames(m), rank]
  lab[is.na(lab) | lab == ""] <- "Unclassified"
  rowsum(rel, group = lab)
}

# Mean abundance and prevalence of every taxon
abundance_table <- function(rel, community, rank) {
  tab <- data.frame(community = community, rank = rank, taxon = rownames(rel),
                    mean_rel_abundance_pct = 100 * rowMeans(rel),
                    sd_pct     = 100 * apply(rel, 1, sd),
                    median_pct = 100 * apply(rel, 1, median),
                    min_pct    = 100 * apply(rel, 1, min),
                    max_pct    = 100 * apply(rel, 1, max),
                    prevalence_pct = 100 * rowMeans(rel > 0),
                    row.names = NULL)
  tab <- tab[order(-tab$mean_rel_abundance_pct), ]
  tab[-(1:3)] <- lapply(tab[-(1:3)], round, 4)
  tab
}

top_taxa <- function(rel, n) {
  m <- sort(rowMeans(rel), decreasing = TRUE)
  head(setdiff(names(m), "Unclassified"), n)
}

# Long table for one panel: taxa outside 'keep' are pooled as "Other"
panel_data <- function(rel, keep, order_by) {
  shown <- intersect(keep, rownames(rel))
  out   <- rel[shown, , drop = FALSE]
  if ("Unclassified" %in% rownames(rel)) out <- rbind(out, Unclassified = rel["Unclassified", ])
  other <- 1 - colSums(out)
  if (any(other > 1e-9)) out <- rbind(out, Other = pmax(other, 0))
  ord <- if (order_by %in% rownames(rel)) order(-rel[order_by, ]) else seq_len(ncol(rel))
  data.frame(Sample = factor(rep(colnames(out), each = nrow(out)), levels = colnames(rel)[ord]),
             Taxon = rep(rownames(out), times = ncol(out)),
             Abundance = as.vector(out))
}

bar_panel <- function(df, colours, legend_title) {
  df$Taxon <- factor(df$Taxon, levels = names(colours))
  ggplot(df, aes(x = Sample, y = Abundance, fill = Taxon)) +
    geom_col(position = "stack", show.legend = TRUE) +
    scale_fill_manual(values = colours, drop = FALSE, name = legend_title) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
    labs(x = "Sample", y = "Relative Abundance") +
    theme_minimal(base_size = 13) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(), legend.position = "right")
}

## -------------------------------- Input --------------------------------------
ps_skin  <- readRDS(file.path(RDS_DIR, SKIN_RDS))
ps_water <- readRDS(file.path(RDS_DIR, WATER_RDS))
message("Skin: ", nsamples(ps_skin), " samples, ", ntaxa(ps_skin), " ASVs; water: ",
        nsamples(ps_water), " samples, ", ntaxa(ps_water), " ASVs")

###############################################################################
# One figure and one set of tables per rank
###############################################################################
make_figure <- function(rank, top_n, base_colours, file_stub) {
  rel_skin  <- rel_by_rank(ps_skin, rank)
  rel_water <- rel_by_rank(ps_water, rank)

  # Tables: all taxa, each community
  tab <- rbind(abundance_table(rel_skin, "skin", rank), abundance_table(rel_water, "water", rank))
  write_tsv(tab, paste0(tolower(rank), "_abundance.tsv"))
  message("\n", rank, ", skin (mean % of reads, top ", min(top_n, nrow(rel_skin)), "):")
  print(head(tab[tab$community == "skin", c("taxon", "mean_rel_abundance_pct", "sd_pct", "prevalence_pct")], top_n), row.names = FALSE)
  message(rank, ", water (mean % of reads, top ", min(top_n, nrow(rel_water)), "):")
  print(head(tab[tab$community == "water", c("taxon", "mean_rel_abundance_pct", "sd_pct", "prevalence_pct")], top_n), row.names = FALSE)

  # Taxa shown: the most abundant of each community; one legend for both panels
  keep_skin  <- top_taxa(rel_skin, top_n)
  keep_water <- top_taxa(rel_water, top_n)
  shown <- sort(union(keep_skin, keep_water))
  colours <- base_colours[intersect(names(base_colours), shown)]
  unnamed <- setdiff(shown, names(colours))
  pool    <- if (is.null(names(base_colours))) base_colours else EXTRA_COLOURS
  if (length(unnamed) > length(pool)) stop("Not enough colours for ", length(unnamed), " taxa; lower the top-N setting.")
  colours <- c(colours, setNames(pool[seq_along(unnamed)], unnamed))
  colours <- c(colours[sort(names(colours))], Unclassified = "grey60", Other = "grey85")

  # Samples ordered by the most abundant taxon of the skin community
  lead <- names(which.max(rowMeans(rel_skin[setdiff(rownames(rel_skin), "Unclassified"), , drop = FALSE])))
  d_water <- panel_data(rel_water, keep_water, lead)
  d_skin  <- panel_data(rel_skin,  keep_skin,  lead)
  used <- names(colours) %in% c(as.character(d_water$Taxon), as.character(d_skin$Taxon))
  colours <- colours[used]

  p_water <- bar_panel(d_water, colours, rank)
  p_skin  <- bar_panel(d_skin,  colours, rank)
  fig <- p_water + p_skin +
    plot_layout(widths = c(1, 2), guides = "collect") +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 16),
          legend.key.size = unit(if (length(colours) > 12) 0.4 else 0.6, "cm"))

  ggsave(file.path(OUT_DIR, paste0(file_stub, ".png")), fig, width = 10, height = 5, dpi = 600, bg = "white")
  ggsave(file.path(OUT_DIR, paste0(file_stub, ".pdf")), fig, width = 10, height = 5)
  message("Samples ordered by ", lead, "; figure: ", file.path(OUT_DIR, paste0(file_stub, ".png")))
  invisible(fig)
}

fig_phylum <- make_figure("Phylum", TOP_N_PHYLUM, PHYLUM_COLOURS, "relative_abundance_phylum")
fig_genus  <- make_figure("Genus",  TOP_N_GENUS,  GENUS_COLOURS,  "relative_abundance_genus")
print(fig_phylum)

writeLines(c(paste("skin object:", SKIN_RDS, "-", nsamples(ps_skin), "samples,", ntaxa(ps_skin), "ASVs"),
             paste("water object:", WATER_RDS, "-", nsamples(ps_water), "samples,", ntaxa(ps_water), "ASVs"),
             capture.output(sessionInfo())), file.path(OUT_DIR, "abundance_run_record.txt"))
message("\nRelative abundance: tables and figures in ", OUT_DIR)
