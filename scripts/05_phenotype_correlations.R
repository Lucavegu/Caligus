#!/usr/bin/env Rscript

###############################################################################
# Skin microbiome & Caligus: description of the fish and correlations between
# body size and sea lice measures
# Replaces Correlations.R. Needs Skin_ps_SILVA.rds and Skin_rare_SILVA.rds from
# 02_qc_decontam.R.
#
# What it adds to the original script:
#   - runs on its own (the original needed Skin_ps, DIRS and vars_to_summarize
#     from other sessions, and used ggcorrplot without loading it)
#   - data checks: lice counts that do not add up, fish from another tank,
#     implausible weight changes, missing values
#   - lice per cm (the trait used in the prediction analysis) is included next
#     to the weight-scaled lice density, so the size dependence of each lice
#     measure can be compared                                    [R1-7, R2-11]
#   - Pearson and Spearman correlations with confidence intervals and
#     Benjamini-Hochberg adjusted p-values                        [R1-15]
#   - results for the fish analysed (rarefied dataset) and for all sequenced fish
###############################################################################

## --------------------------- User Parameters ---------------------------------
PROJECT_DIR <- getwd()
RDS_DIR     <- file.path(PROJECT_DIR)
OUT_DIR     <- file.path(PROJECT_DIR, "phenotype")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# A final weight outside this range of the initial weight is listed for checking
WEIGHT_RATIO_RANGE <- c(0.7, 1.5)

## ------------------------------ Libraries ------------------------------------
suppressPackageStartupMessages({
  library(phyloseq)
  library(ggplot2)
})

## ------------------------------- Helpers -------------------------------------
write_tsv <- function(x, name) {
  write.table(x, file.path(OUT_DIR, name), sep = "\t", quote = FALSE, row.names = FALSE)
}
to_num <- function(v) suppressWarnings(as.numeric(as.character(v)))

## -------------------------------- Input --------------------------------------
ps_all  <- readRDS(file.path(RDS_DIR, "Skin_ps_SILVA.rds"))       # all sequenced skin fish
ps_rare <- readRDS(file.path(RDS_DIR, "Skin_rare_SILVA.rds"))     # fish analysed
sd0 <- as(sample_data(ps_all), "data.frame")

num_cols <- intersect(c("Initial_Weight", "Final_Weight", "Final_Length", "Fin_caligus", "Body_caligus", "Total_caligus"),
                      colnames(sd0))
d <- data.frame(sample = rownames(sd0), lapply(sd0[, num_cols, drop = FALSE], to_num), row.names = rownames(sd0))
for (v in intersect(c("Sex", "Tank", "Fish_family", "seq_batch"), colnames(sd0))) d[[v]] <- as.character(sd0[[v]])
d$analysed <- d$sample %in% sample_names(ps_rare)

# Lice measures. Lice per cm is the trait used to define burden classes;
# lice density scales the count by weight^(2/3), a proxy for body surface.
d$Lice_per_cm      <- d$Total_caligus / d$Final_Length
d$Lice_density     <- d$Total_caligus / d$Final_Weight^(2 / 3)
d$Log_lice_count   <- log(d$Total_caligus + 1)
d$Log_lice_density <- log((d$Total_caligus + 1) / d$Final_Weight^(2 / 3))

###############################################################################
# 1. DATA CHECKS
###############################################################################
checks <- NULL
add_check <- function(what, rows, detail) {
  if (nrow(rows) == 0) return(invisible())
  checks <<- rbind(checks, data.frame(check = what, sample = rows$sample, detail = detail, analysed = rows$analysed))
}
# Fin and body counts must add up to the total. Where they do not, the split is
# treated as missing (the total is kept).
if (all(c("Fin_caligus", "Body_caligus") %in% colnames(d))) {
  bad <- which(d$Fin_caligus + d$Body_caligus != d$Total_caligus)
  add_check("fin + body lice differ from total", d[bad, ],
            paste0("fin ", d$Fin_caligus[bad], ", body ", d$Body_caligus[bad], ", total ", d$Total_caligus[bad]))
  d$Fin_caligus[bad] <- NA
  d$Body_caligus[bad] <- NA
}
if ("Tank" %in% colnames(d)) {
  main_tank <- names(which.max(table(d$Tank)))
  other <- which(d$Tank != main_tank)
  add_check(paste("tank other than", main_tank), d[other, ], paste("tank", d$Tank[other]))
}
if (all(c("Initial_Weight", "Final_Weight") %in% colnames(d))) {
  ratio <- d$Final_Weight / d$Initial_Weight
  odd <- which(ratio < WEIGHT_RATIO_RANGE[1] | ratio > WEIGHT_RATIO_RANGE[2])
  add_check("large weight change between initial and final weighing", d[odd, ],
            paste0(d$Initial_Weight[odd], " g to ", d$Final_Weight[odd], " g (x", round(ratio[odd], 2), ")"))
  message("Fish lighter at the final than at the initial weighing: ", sum(ratio < 1, na.rm = TRUE), " of ", sum(!is.na(ratio)))
}
miss <- which(!complete.cases(d[, c("Final_Weight", "Final_Length", "Total_caligus")]))
add_check("missing weight, length or lice count", d[miss, ], "missing value")

if (is.null(checks)) {
  message("Data checks: nothing to report.")
} else {
  write_tsv(checks, "data_checks.tsv")
  message("Data checks (phenotype/data_checks.tsv):")
  print(checks, row.names = FALSE)
}

###############################################################################
# 2. DESCRIPTION AND CORRELATIONS, for the fish analysed and for all fish
###############################################################################
VARS <- intersect(c("Initial_Weight", "Final_Weight", "Final_Length", "Fin_caligus", "Body_caligus", "Total_caligus",
                    "Lice_per_cm", "Lice_density", "Log_lice_count", "Log_lice_density"), colnames(d))
SIZE <- intersect(c("Final_Length", "Final_Weight", "Initial_Weight"), VARS)
LICE <- intersect(c("Total_caligus", "Lice_per_cm", "Lice_density"), VARS)
SETS <- list(analysed = d[d$analysed, ], all_sequenced = d)

cor_pair <- function(x, y) {
  ok <- complete.cases(x, y)
  pe <- cor.test(x[ok], y[ok], method = "pearson")
  sp <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))
  data.frame(n = sum(ok), pearson_r = unname(pe$estimate), pearson_ci_low = pe$conf.int[1], pearson_ci_high = pe$conf.int[2],
             pearson_p = pe$p.value, spearman_rho = unname(sp$estimate), spearman_p = sp$p.value)
}

for (sn in names(SETS)) {
  x <- SETS[[sn]]
  message("\n=== ", sn, ": ", nrow(x), " fish ===")

  # Descriptive statistics
  desc <- do.call(rbind, lapply(VARS, function(v) {
    z <- x[[v]]
    data.frame(variable = v, n = sum(!is.na(z)), mean = mean(z, na.rm = TRUE), sd = sd(z, na.rm = TRUE),
               median = median(z, na.rm = TRUE), min = min(z, na.rm = TRUE), max = max(z, na.rm = TRUE),
               cv_pct = 100 * sd(z, na.rm = TRUE) / mean(z, na.rm = TRUE))
  }))
  desc[-1] <- lapply(desc[-1], signif, 4)
  write_tsv(desc, paste0("descriptive_statistics_", sn, ".tsv"))
  print(desc, row.names = FALSE)
  for (v in intersect(c("Sex", "Tank", "seq_batch"), colnames(x))) {
    message(v, ": ", paste(names(table(x[[v]])), table(x[[v]]), collapse = ", "))
  }
  if ("Fish_family" %in% colnames(x)) {
    ft <- table(table(x$Fish_family))
    message("Families: ", length(unique(x$Fish_family)), " (fish per family: ",
            paste(paste0(names(ft), " fish x ", ft), collapse = ", "), ")")
  }

  # All pairwise correlations, with p-values adjusted over the pairs
  pairs <- t(combn(VARS, 2))
  ct <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) {
    cbind(data.frame(var_1 = pairs[i, 1], var_2 = pairs[i, 2]), cor_pair(x[[pairs[i, 1]]], x[[pairs[i, 2]]]))
  }))
  ct$pearson_p_BH  <- p.adjust(ct$pearson_p, method = "BH")
  ct$spearman_p_BH <- p.adjust(ct$spearman_p, method = "BH")
  out <- ct
  out[-(1:3)] <- lapply(out[-(1:3)], signif, 3)
  write_tsv(out, paste0("correlations_all_pairs_", sn, ".tsv"))

  # Body size against each lice measure. Lice density divides the count by a
  # function of weight, so a negative correlation with weight can arise from
  # the scaling itself.
  key <- out[(out$var_1 %in% SIZE & out$var_2 %in% LICE) | (out$var_1 %in% LICE & out$var_2 %in% SIZE), ]
  write_tsv(key, paste0("size_vs_lice_", sn, ".tsv"))
  message("Body size against lice measures:")
  print(key[, c("var_1", "var_2", "n", "pearson_r", "pearson_ci_low", "pearson_ci_high", "pearson_p_BH", "spearman_rho", "spearman_p_BH")],
        row.names = FALSE)

  # Lice measures by sex
  if ("Sex" %in% colnames(x) && length(unique(x$Sex)) == 2) {
    sx <- do.call(rbind, lapply(c(SIZE, LICE), function(v) {
      g <- split(x[[v]], x$Sex)
      data.frame(variable = v, group_1 = names(g)[1], median_1 = median(g[[1]], na.rm = TRUE),
                 group_2 = names(g)[2], median_2 = median(g[[2]], na.rm = TRUE),
                 wilcoxon_p = suppressWarnings(wilcox.test(g[[1]], g[[2]])$p.value))
    }))
    sx$wilcoxon_p_BH <- p.adjust(sx$wilcoxon_p, method = "BH")
    sx[c("median_1", "median_2", "wilcoxon_p", "wilcoxon_p_BH")] <- lapply(sx[c("median_1", "median_2", "wilcoxon_p", "wilcoxon_p_BH")], signif, 3)
    write_tsv(sx, paste0("size_and_lice_by_sex_", sn, ".tsv"))
  }

  # Correlation heat maps (lower triangle; * = adjusted p < 0.05)
  heat <- function(r_col, p_col, title) {
    h <- rbind(data.frame(a = ct$var_1, b = ct$var_2, r = ct[[r_col]], p = ct[[p_col]]))
    h$a <- factor(h$a, levels = VARS); h$b <- factor(h$b, levels = rev(VARS))
    h$lab <- paste0(sprintf("%.2f", h$r), ifelse(h$p < 0.05, "*", ""))
    ggplot(h, aes(x = a, y = b, fill = r)) +
      geom_tile(color = "white") +
      geom_text(aes(label = lab), size = 3) +
      scale_fill_gradient2(low = "#2166ac", mid = "#ffffff", high = "#b2182b", midpoint = 0, limits = c(-1, 1), name = "r") +
      labs(x = NULL, y = NULL, title = title, subtitle = paste0(nrow(x), " fish; * adjusted p < 0.05 (Benjamini-Hochberg)")) +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), panel.grid = element_blank())
  }
  p_pe <- heat("pearson_r", "pearson_p_BH", "Pearson correlation")
  p_sp <- heat("spearman_rho", "spearman_p_BH", "Spearman correlation")
  ggsave(file.path(OUT_DIR, paste0("correlation_pearson_", sn, ".png")), p_pe, width = 7.5, height = 6.5, dpi = 300)
  ggsave(file.path(OUT_DIR, paste0("correlation_spearman_", sn, ".png")), p_sp, width = 7.5, height = 6.5, dpi = 300)

  # Each lice measure against length and weight
  long <- do.call(rbind, lapply(LICE, function(lv) do.call(rbind, lapply(intersect(c("Final_Length", "Final_Weight"), SIZE), function(sv) {
    data.frame(lice_measure = lv, size_measure = sv, size = x[[sv]], lice = x[[lv]])
  }))))
  long$lice_measure <- factor(long$lice_measure, levels = LICE)
  p_sc <- ggplot(long, aes(x = size, y = lice)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 0.6, color = "black") +
    facet_grid(lice_measure ~ size_measure, scales = "free") +
    theme_bw(base_size = 12) +
    labs(title = "Lice measures against body size", subtitle = paste0(nrow(x), " fish; line: linear fit with 95% band"),
         x = "Body size (length in cm, weight in g)", y = NULL)
  pdf(file.path(OUT_DIR, paste0("phenotype_figures_", sn, ".pdf")), width = 8, height = 8)
  print(p_pe); print(p_sp); print(p_sc)
  dev.off()
}

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
message("\nPhenotype description: tables and figures in ", OUT_DIR)
