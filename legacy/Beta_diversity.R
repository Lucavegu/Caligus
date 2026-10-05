# ================= NMDS: Bray–Curtis & Jaccard Dissimilarities =================
library(phyloseq)
library(vegan)
library(ggplot2)
library(patchwork)
library(dplyr)

# Use your phyloseq object
ps <- High_Low_microbiome
sample_data(ps)$Burden <- factor(sample_data(ps)$Burden)

# Colorblind-safe palette (Okabe-Ito)
cb2 <- c("Low" = "#0072B2", "High" = "#D55E00")

# Compute distances
bray_dist    <- phyloseq::distance(ps, method = "bray")
jaccard_dist <- phyloseq::distance(ps, method = "jaccard", binary = TRUE)

# NMDS ordination
set.seed(123)
nmds_bray    <- ordinate(ps, method = "NMDS", distance = bray_dist,    trymax = 100)
set.seed(123)
nmds_jaccard <- ordinate(ps, method = "NMDS", distance = jaccard_dist, trymax = 100)

# Extract stress values
get_stress <- function(ord) {
  val <- tryCatch({
    if (!is.null(ord$stress)) ord$stress else NA_real_
  }, error = function(e) NA_real_)
  val
}
stress_bray    <- get_stress(nmds_bray)
stress_jaccard <- get_stress(nmds_jaccard)

# -----------------------
# Helper: NMDS plot
# -----------------------
make_nmds_plot <- function(ps_obj, ord_obj, palette, stress_val) {
  plot_ordination(ps_obj, ord_obj, color = "Burden") +
    geom_point(size = 3, alpha = 0.85) +
    stat_ellipse(
      aes(group = Burden, color = Burden, fill = Burden),
      type = "t", level = 0.95, geom = "polygon",
      alpha = 0.2, linewidth = 0.5
    ) +
    scale_color_manual(values = palette, name = "Burden") +
    scale_fill_manual(values = palette, name = "Burden") +
    labs(x = "NMDS1", y = "NMDS2") +
    theme_classic(base_size = 13) +
    theme(
      legend.position = "right",
      plot.title = element_blank(),
      panel.grid = element_blank()
    ) +
    annotate(
      "text", x = -Inf, y = Inf,
      label = sprintf("Stress = %.3f", stress_val),
      hjust = -0.05, vjust = 1.3, size = 3.5
    )
}

# Create plots
p_bray    <- make_nmds_plot(ps, nmds_bray,    cb2, stress_bray)
p_jaccard <- make_nmds_plot(ps, nmds_jaccard, cb2, stress_jaccard)

# Combine with panel letters
combined_nmds <- p_bray + p_jaccard +
  plot_layout(ncol = 2) &
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 16))

# Show combined plot
combined_nmds

# Save figure
ggsave("nmds_bray_jaccard_ellipses.png", combined_nmds, width = 8, height = 5, dpi = 600)
ggsave("nmds_bray_jaccard_ellipses.pdf", combined_nmds, width = 12, height = 5.5)

# Load required libraries
library(phyloseq)
library(vegan)
library(ggplot2)
library(dplyr)

# Use your phyloseq object
ps <- High_Low_microbiome
meta <- data.frame(sample_data(ps))
meta$Burden <- factor(meta$Burden)

# Colorblind palette
cb2 <- c("Low" = "#0072B2", "High" = "#D55E00")


##Statistics
# -------------------------------
# 1. Calculate distance matrices
# -------------------------------
bray_dist    <- phyloseq::distance(ps, method = "bray")
jaccard_dist <- phyloseq::distance(ps, method = "jaccard", binary = TRUE)

# -------------------------------
# 2. PERMANOVA tests (adonis2)
# -------------------------------
set.seed(123)
permanova_bray <- adonis2(bray_dist ~ Burden, data = meta, permutations = 999)
set.seed(123)
permanova_jaccard <- adonis2(jaccard_dist ~ Burden, data = meta, permutations = 999)

cat("\n### PERMANOVA — Bray-Curtis ###\n")
print(permanova_bray)
cat("\n### PERMANOVA — Jaccard ###\n")
print(permanova_jaccard)

# Save results
write.table(as.data.frame(permanova_bray), "permanova_bray.tsv", sep = "\t", quote = FALSE)
write.table(as.data.frame(permanova_jaccard), "permanova_jaccard.tsv", sep = "\t", quote = FALSE)

# ------------------------------------------------
# 3. Homogeneity of multivariate dispersion check
# ------------------------------------------------
# Bray-Curtis
bray_betadisp <- betadisper(bray_dist, meta$Burden)
anova_bray_disp <- anova(bray_betadisp)
perm_bray_disp  <- permutest(bray_betadisp, permutations = 999)

cat("\n### Homogeneity of dispersion — Bray-Curtis ###\n")
print(anova_bray_disp)
print(perm_bray_disp)

# Jaccard
jaccard_betadisp <- betadisper(jaccard_dist, meta$Burden)
anova_jaccard_disp <- anova(jaccard_betadisp)
perm_jaccard_disp  <- permutest(jaccard_betadisp, permutations = 999)

cat("\n### Homogeneity of dispersion — Jaccard ###\n")
print(anova_jaccard_disp)
print(perm_jaccard_disp)

# ------------------------------------------------
# 4. Visualize dispersion (optional, supplementary figure)
# ------------------------------------------------
p_bray_disp <- ggplot(data.frame(Distance = bray_betadisp$distances,
                                 Burden = meta$Burden),
                      aes(x = Burden, y = Distance, fill = Burden)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.5, colour = "black") +
  geom_jitter(width = 0.15, size = 1.6, alpha = 0.6, aes(color = Burden)) +
  scale_fill_manual(values = cb2) +
  scale_colour_manual(values = cb2) +
  labs(y = "Distance to centroid", x = NULL, title = "Bray-Curtis dispersion") +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "none")

p_jaccard_disp <- ggplot(data.frame(Distance = jaccard_betadisp$distances,
                                    Burden = meta$Burden),
                         aes(x = Burden, y = Distance, fill = Burden)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.5, colour = "black") +
  geom_jitter(width = 0.15, size = 1.6, alpha = 0.6, aes(color = Burden)) +
  scale_fill_manual(values = cb2) +
  scale_colour_manual(values = cb2) +
  labs(y = "Distance to centroid", x = NULL, title = "Jaccard dispersion") +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "none")

# Combine dispersion plots
library(patchwork)
combined_disp <- p_bray_disp + p_jaccard_disp +
  plot_layout(ncol = 2) &
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 16))

ggsave("beta_dispersion_bray_jaccard.png", combined_disp, width = 12, height = 5.5, dpi = 600)
ggsave("beta_dispersion_bray_jaccard.pdf", combined_disp, width = 12, height = 5.5)

# ------------------------------------------------
# 5. Summary table
# ------------------------------------------------

permanova_summary <- data.frame(
  Metric = c("Bray-Curtis", "Jaccard"),
  Df = c(permanova_bray$Df[1], permanova_jaccard$Df[1]),
  F_stat = c(permanova_bray$F[1], permanova_jaccard$F[1]),
  R2 = c(permanova_bray$R2[1], permanova_jaccard$R2[1]),
  PERMANOVA_p = c(permanova_bray$`Pr(>F)`[1], permanova_jaccard$`Pr(>F)`[1]),
  Dispersion_p = c(
    perm_bray_disp$tab[1, "Pr(>F)"],
    perm_jaccard_disp$tab[1, "Pr(>F)"]
  )
)

# Add significance codes for PERMANOVA
permanova_summary$Significance <- cut(
  permanova_summary$PERMANOVA_p,
  breaks = c(-Inf, 0.001, 0.01, 0.05, 1),
  labels = c("***", "**", "*", "ns")
)

# Round values for readability
permanova_summary <- permanova_summary %>%
  mutate(across(c(F_stat, R2, PERMANOVA_p, Dispersion_p), ~ round(., 4)))

# Show final summary table
print(permanova_summary)

# Save as TSV
write.table(permanova_summary,
            "permanova_detailed_summary.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

