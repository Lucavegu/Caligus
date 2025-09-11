# ---- Requirements -------------------------------------------------------------
library(dplyr)
library(car)        # leveneTest, Anova
library(lme4)       # lmer
library(emmeans)    # EMMs
library(DHARMa)     # simulateResiduals

# ---- Input datasets (already created in your script) --------------------------
# Shannon: data.frame with columns Burden, Fish_family, Shannon
# Pielou : data.frame with columns Burden, Fish_family, pielou
# Faith  : data.frame with columns Burden, Fish_family, PD

# ---- Helper: safe Shapiro-by-group -------------------------------------------
shapiro_by_group_fn <- function(df, group = "Burden", y) {
  df %>%
    group_by(.data[[group]]) %>%
    summarise(
      n = dplyr::n(),
      shapiro_p = if (n >= 3 && n <= 5000) {
        # run test on current group's y values
        shapiro.test(.data[[y]])$p.value
      } else NA_real_,
      note = dplyr::case_when(
        n < 3 ~ "n<3: Shapiro not applicable",
        n > 5000 ~ "n>5000: Shapiro not applicable",
        TRUE ~ "OK"
      ),
      .groups = "drop"
    ) %>%
    rename(!!group := .data[[group]])
}

# ---- Helper: run full workflow for one measure --------------------------------
run_workflow <- function(df, y, burden_col = "Burden", rand_col = "Fish_family") {
  # Ensure factors
  df[[burden_col]] <- as.factor(df[[burden_col]])
  df[[rand_col]]   <- as.factor(df[[rand_col]])
  
  # 1) Shapiro by group
  shapiro_tbl <- shapiro_by_group_fn(df, group = burden_col, y = y)
  
  # 2) Levene’s test
  lev_formula <- reformulate(termlabels = burden_col, response = y)   # y ~ Burden
  lev_res <- car::leveneTest(lev_formula, data = df)
  
  # 3) Mixed model
  mm_formula <- as.formula(paste0(y, " ~ ", burden_col, " + (1 | ", rand_col, ")"))
  mm <- lme4::lmer(mm_formula, data = df)
  mm_sum <- summary(mm)
  
  # 4) Type II ANOVA (Wald χ²)
  anova_res <- car::Anova(mm, type = "II")
  
  # 5) EMMs by Burden (optional for plotting)
  emms <- emmeans::emmeans(mm, specs = as.formula(paste0("~ ", burden_col)))
  emms_tbl <- summary(emms) %>%
    as.data.frame() %>%
    rename(
      Burden = !!burden_col,
      emmean = emmean,
      SE = SE,
      df = df,
      Lower_CI = lower.CL,
      Upper_CI = upper.CL
    )
  
  # 6) DHARMa residual diagnostics
  sim_res <- DHARMa::simulateResiduals(mm, n = 1000)
  uni_test <- DHARMa::testUniformity(sim_res)
  disp_test <- DHARMa::testDispersion(sim_res)
  
  # Pack results
  list(
    shapiro = shapiro_tbl,
    levene  = lev_res,
    model   = mm,
    model_summary = mm_sum,
    anova   = anova_res,
    emmeans = emms_tbl,
    dharma  = list(sim = sim_res, uniformity = uni_test, dispersion = disp_test)
  )
}

# ---- Define measures and run loop --------------------------------------------
measures <- list(
  Shannon = list(data = Shannon, y = "Shannon"),
  Pielou  = list(data = Pielou,  y = "pielou"),
  Faith   = list(data = Faith,   y = "PD")
)

results_by_measure <- lapply(names(measures), function(nm) {
  cat("\n==========", nm, "==========\n")
  res <- run_workflow(measures[[nm]]$data, y = measures[[nm]]$y)
  
  # Compact console readouts
  cat("\nShapiro by Burden:\n"); print(res$shapiro)
  cat("\nLevene’s test:\n"); print(res$levene)
  cat("\nType II ANOVA (Wald χ²):\n"); print(res$anova)
  cat("\nEMMs (by Burden):\n"); print(res$emmeans)
  
  # Quick DHARMa plots (comment out if running non-interactive)
  # plot(res$dharma$sim)
  
  invisible(res)
})
names(results_by_measure) <- names(measures)

# Access examples:
# results_by_measure$Shannon$shapiro
# results_by_measure$Pielou$anova
# results_by_measure$Faith$emmeans

### PLOTTING
# Load libraries
library(ggplot2)
library(dplyr)
library(patchwork)

# Colorblind-safe palette (Okabe–Ito)
cb2 <- c("Low" = "#0072B2", "High" = "#D55E00")

# Define datasets and y-variable names
plots_info <- list(
  Shannon = list(data = Shannon, y = "Shannon", ylab = "Shannon Diversity Index"),
  Pielou  = list(data = Pielou,  y = "pielou",  ylab = "Pielou's Evenness"),
  Faith   = list(data = Faith,   y = "PD",      ylab = "Faith's PD")
)

# Function to create individual boxplots (no title)
make_boxplot <- function(df, y_var, y_lab) {
  
  # Ensure Burden is a factor and ordered
  df$Burden <- factor(df$Burden, levels = names(cb2))
  
  # Calculate group means for dotted lines
  group_means <- df %>%
    group_by(Burden) %>%
    summarise(mean_val = mean(.data[[y_var]], na.rm = TRUE), .groups = "drop")
  
  # Create plot
  p <- ggplot(df, aes(x = Burden, y = .data[[y_var]], fill = Burden)) +
    geom_boxplot(
      width = 0.5,
      alpha = 0.7,
      colour = "black",
      linewidth = 0.6,
      outlier.shape = NA
    ) +
    geom_segment(
      data = group_means,
      aes(
        x = as.numeric(Burden) - 0.25,
        xend = as.numeric(Burden) + 0.25,
        y = mean_val,
        yend = mean_val
      ),
      inherit.aes = FALSE,
      colour = "black",
      linetype = "dotted",
      linewidth = 0.7
    ) +
    geom_jitter(
      aes(color = Burden),
      width = 0.15,
      size = 1.5,
      alpha = 0.6,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = cb2) +
    scale_colour_manual(values = cb2) +
    labs(x = NULL, y = y_lab) +
    theme_classic(base_size = 13) +
    theme(
      axis.text.x = element_text(face = "bold"),
      axis.title.x = element_blank(),
      plot.title = element_blank(),
      legend.position = "none"
    )
  
  return(p)
}

# Generate individual plots
p1 <- make_boxplot(Shannon, "Shannon", "Shannon Diversity Index")
p2 <- make_boxplot(Pielou,  "pielou",  "Pielou's Evenness")
p3 <- make_boxplot(Faith,   "PD",      "Faith's PD")

# Combine into a multi-panel figure with letters
combined_plot <- p1 + p2 + p3 +
  plot_layout(ncol = 3) &
  plot_annotation(tag_levels = "A") &           # Labels as A, B, C
  theme(plot.tag = element_text(face = "bold", size = 16))

# Show combined figure
combined_plot

# Save high-resolution figure
ggsave("combined_boxplots.png", combined_plot, width = 7, height = 5, dpi = 600)
ggsave("combined_boxplots.pdf", combined_plot, width = 14, height = 5)
