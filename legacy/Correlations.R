
############################################
# STATISTICAL SUMMARY REPORT FOR METADATA
############################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(ggthemes)
  library(RColorBrewer)
})

#Convert sample_data -> data.frame 
metadata_df <- as(phyloseq::sample_data(Skin_ps), "data.frame")

#2) Coerce selected columns to numeric
metadata_num <- metadata_df %>%
  dplyr::select(all_of(vars_to_summarize)) %>%
  mutate(across(everything(), ~ suppressWarnings(as.numeric(.))))

# Summarise
desc_stats <- metadata_num %>%
  summarise(across(
    everything(),
    list(
      mean   = ~mean(., na.rm = TRUE),
      sd     = ~sd(., na.rm = TRUE),
      median = ~median(., na.rm = TRUE),
      min    = ~min(., na.rm = TRUE),
      max    = ~max(., na.rm = TRUE)
    ),
    .names = "{.col}_{.fn}"
  ))

# FIX: Use regex to split at the LAST underscore only
desc_stats_tidy <- desc_stats %>%
  tidyr::pivot_longer(
    cols = everything(),
    names_to = c("Variable", "Statistic"),
    names_pattern = "^(.*)_([^_]*)$",  # <-- key fix
    values_to = "Value"
  )

# Check result
head(desc_stats_tidy)

# Save summary statistics
write.table(desc_stats_tidy,
            file = file.path(DIRS$tables, "metadata_summary_stats.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)


###############################
# Prepare numeric metadata again to be safe
metadata_df <- as(phyloseq::sample_data(Skin_ps), "data.frame")
vars <- c("Initial_Weight","Final_Weight","Final_Length",
          "Fin_caligus","Body_caligus","Total_caligus",
          "LogLC","liceD","LogLiceD")
vars <- intersect(vars, colnames(metadata_df))


num_df <- metadata_df %>%
  dplyr::select(all_of(vars)) %>%
  mutate(across(everything(), ~ suppressWarnings(as.numeric(as.character(.)))))


# Pearson correlation matrix
R <- cor(num_df, use = "pairwise.complete.obs", method = "pearson")


# p-values for reference, but we won't blank cells anymore
p.mat <- matrix(NA_real_, ncol = ncol(num_df), nrow = ncol(num_df),
                dimnames = list(colnames(num_df), colnames(num_df)))
for (i in seq_len(ncol(num_df))) {
  for (j in seq_len(ncol(num_df))) {
    if (i == j) {
      p.mat[i, j] <- 0
    } else {
      p.mat[i, j] <- suppressWarnings(cor.test(num_df[[i]], num_df[[j]], method = "pearson")$p.value)
    }
  }
}


# Pearson correlation plot — always show all r values
pearson_plot <- ggcorrplot(R,
                           type = "lower",
                           hc.order = TRUE,
                           lab = TRUE,
                           lab_size = 3,
                           outline.color = "white",
                           colors = c("#2166ac","#ffffff","#b2182b")) +
  ggtitle("") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    panel.grid = element_blank()
  )
pearson_plot

ggsave(file.path(DIRS$tables, "skin_phenotype_correlation_pearson_all.png"),
       pearson_plot, width = 7.5, height = 6.5, dpi = 300)

