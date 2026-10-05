#!/usr/bin/env Rscript

###############################################################################
# Skin microbiome & Caligus: alpha and beta diversity against sea lice burden
# Replaces Alpha_diversity.R and Beta_diversity.R. Needs the outputs of
# 02_qc_decontam.R.
#
# What it adds to the original scripts (reviewer comment in brackets):
#   - runs on its own: burden classes and the diversity indices are computed
#     here (the original needed objects created elsewhere)
#   - alpha diversity is computed on all ASVs (before the prevalence filter)
#     as well as on the filtered table, where richness is capped at ~124 taxa
#   - family is a random effect and sequencing run a covariate in the alpha
#     models; run and sequencing depth are covariates in PERMANOVA [R1-6, R2-10]
#   - PERMANOVA table with Df, pseudo-F, R2, exact permutation p-values and
#     the test of homogeneity of dispersion                        [R2-9]
#   - Aitchison distance (CLR) next to Bray-Curtis and Jaccard     [R1-8]
#   - the same tests for the company categories and for lice load as a
#     continuous trait                                             [R1-4, R2-4]
#   - p-values adjusted across the alpha indices (Benjamini-Hochberg) [R1-15]
###############################################################################

## --------------------------- User Parameters ---------------------------------
PROJECT_DIR <- getwd()
RDS_DIR     <- file.path(PROJECT_DIR)
DIV_DIR     <- file.path(PROJECT_DIR, "diversity")
if (!dir.exists(DIV_DIR)) dir.create(DIV_DIR, recursive = TRUE)

FAMILY_COL   <- "Fish_family"
RAREFY_DEPTH <- 10000     # for the all-ASV table (same depth as the filtered table)
RAREFY_SEED  <- 123
N_PERM       <- 9999      # permutations (smallest possible p-value: 0.0001)
PERM_SEED    <- 123

# Company-defined burden categories on the raw lice count
COMPANY_LOW  <- 15        # Low  : fewer than this
COMPANY_HIGH <- 30        # High : more than this

## ------------------------------ Libraries ------------------------------------
suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(ape)
  library(VGAM)
  library(lme4)
  library(car)
  library(ggplot2)
  library(patchwork)
})
HAS_DHARMA <- requireNamespace("DHARMa", quietly = TRUE)
message("vegan version: ", as.character(packageVersion("vegan")))

## ------------------------------- Helpers -------------------------------------
otu_mat <- function(x) {            # counts with ASVs in rows
  m <- as(otu_table(x), "matrix")
  if (!taxa_are_rows(x)) m <- t(m)
  m
}
write_tsv <- function(x, name) {
  write.table(x, file.path(DIV_DIR, name), sep = "\t", quote = FALSE, row.names = FALSE)
}
to_num <- function(v) suppressWarnings(as.numeric(as.character(v)))
cb2 <- c("Low" = "#0072B2", "High" = "#D55E00")     # colour-blind-safe (Okabe-Ito)

###############################################################################
# 1. DATA AND BURDEN DEFINITIONS (same rules as 06_ml_burden_prediction.R)
###############################################################################
ps_rare <- readRDS(file.path(RDS_DIR, "Skin_rare_SILVA.rds"))             # filtered, rarefied
f_unf <- file.path(RDS_DIR, "bacteria_unfiltered_physeq.rds")             # before the prevalence filter
f_clr <- file.path(RDS_DIR, "Skin_clr_SILVA.rds")
f_raw <- file.path(RDS_DIR, "Skin_ps_SILVA.rds")
fish  <- sample_names(ps_rare)

sd0  <- as(sample_data(ps_rare), "data.frame")
meta <- data.frame(
  sample = fish,
  lice   = to_num(sd0$Total_caligus),
  length = to_num(sd0$Final_Length),
  family = if (FAMILY_COL %in% colnames(sd0)) as.character(sd0[[FAMILY_COL]]) else paste0("fish_", fish),
  run    = if ("seq_batch" %in% colnames(sd0)) as.character(sd0$seq_batch) else "run1",
  row.names = fish, stringsAsFactors = FALSE
)
meta$family[is.na(meta$family)] <- paste0("fish_", meta$sample[is.na(meta$family)])
meta$lice_load <- meta$lice / meta$length
meta$log_depth <- if (file.exists(f_raw)) log10(sample_sums(readRDS(f_raw))[fish]) else NA_real_

# Gaussian-mixture classes
fit  <- vglm(meta$lice_load ~ 1, mix2normal(eq.sd = FALSE), trace = FALSE)
pars <- as.vector(coef(fit))
mix  <- list(m1 = pars[2], sd1 = exp(pars[3]), m2 = pars[4], sd2 = exp(pars[5]))
if (mix$m1 > mix$m2) mix <- list(m1 = mix$m2, sd1 = mix$sd2, m2 = mix$m1, sd2 = mix$sd1)
LOW_CUT  <- mix$m1 + mix$sd1
HIGH_CUT <- mix$m2 - mix$sd2
meta$Burden <- ifelse(meta$lice_load <= LOW_CUT, "Low", ifelse(meta$lice_load >= HIGH_CUT, "High", "Moderate"))
meta$Burden_company <- ifelse(meta$lice < COMPANY_LOW, "Low", ifelse(meta$lice > COMPANY_HIGH, "High", "Intermediate"))
message("Gaussian cut-offs: Low <= ", round(LOW_CUT, 4), ", High >= ", round(HIGH_CUT, 4), " lice per cm; ",
        paste(names(table(meta$Burden)), table(meta$Burden), collapse = ", "))
message("Company categories: ", paste(names(table(meta$Burden_company)), table(meta$Burden_company), collapse = ", "))

# Outcome definitions: which fish, and the grouping or continuous variable
OUTCOMES <- list(
  gaussian   = list(fish = fish[meta$Burden != "Moderate"], var = "Burden", type = "group",
                    label = "Gaussian-derived Low vs High (lice per cm)"),
  company    = list(fish = fish[meta$Burden_company != "Intermediate"], var = "Burden_company", type = "group",
                    label = paste0("Company categories (Low < ", COMPANY_LOW, ", High > ", COMPANY_HIGH, " lice)")),
  continuous = list(fish = fish, var = "lice_load", type = "continuous",
                    label = "Lice per cm, all fish")
)

###############################################################################
# 2. COUNT TABLES
###############################################################################
# "filtered": the table used everywhere else (prevalence-filtered, rarefied)
# "all_ASVs": the same fish before the prevalence filter, rarefied to the same depth
tables <- list(filtered = list(counts = t(otu_mat(ps_rare)), tree = phy_tree(ps_rare, errorIfNULL = FALSE)))
if (file.exists(f_unf)) {
  ps_unf <- prune_samples(fish, readRDS(f_unf))
  short  <- sample_names(ps_unf)[sample_sums(ps_unf) < RAREFY_DEPTH]
  if (length(short)) stop("Fish below ", RAREFY_DEPTH, " reads in the unfiltered table: ", paste(short, collapse = ", "))
  ps_unf <- rarefy_even_depth(ps_unf, sample.size = RAREFY_DEPTH, rngseed = RAREFY_SEED,
                              replace = FALSE, trimOTUs = TRUE, verbose = FALSE)
  tables$all_ASVs <- list(counts = t(otu_mat(ps_unf))[fish, , drop = FALSE], tree = phy_tree(ps_unf, errorIfNULL = FALSE))
  message("All-ASV table: ", ntaxa(ps_unf), " ASVs after rarefying to ", RAREFY_DEPTH, " reads")
} else {
  message("bacteria_unfiltered_physeq.rds not found: only the filtered table is analysed.")
}
tables$filtered$counts <- tables$filtered$counts[fish, , drop = FALSE]
clr <- if (file.exists(f_clr)) t(otu_mat(readRDS(f_clr)))[fish, , drop = FALSE] else NULL

###############################################################################
# 3. ALPHA DIVERSITY
###############################################################################
# Faith's phylogenetic diversity: total branch length linking the ASVs present
# in a sample, including the path to the root (as picante::pd with
# include.root = TRUE)
faith_pd <- function(counts, tree) {
  tree <- ape::reorder.phylo(ape::keep.tip(tree, colnames(counts)), "postorder")
  ntip <- length(tree$tip.label)
  pres <- matrix(FALSE, ntip + tree$Nnode, nrow(counts))
  pres[seq_len(ntip), ] <- t(counts[, tree$tip.label, drop = FALSE] > 0)
  for (i in seq_len(nrow(tree$edge))) {             # children are visited before their parents
    pres[tree$edge[i, 1], ] <- pres[tree$edge[i, 1], ] | pres[tree$edge[i, 2], ]
  }
  colSums(pres[tree$edge[, 2], , drop = FALSE] * tree$edge.length)
}

alpha_all <- do.call(rbind, lapply(names(tables), function(tn) {
  cnt <- tables[[tn]]$counts
  obs <- rowSums(cnt > 0)
  sh  <- vegan::diversity(cnt, index = "shannon")
  tr  <- tables[[tn]]$tree
  pd  <- if (!is.null(tr)) {
    if (!ape::is.rooted(tr)) warning("Tree of table '", tn, "' is not rooted; Faith's PD depends on the root.")
    faith_pd(cnt, tr)
  } else NA_real_
  data.frame(table = tn, sample = rownames(cnt), Observed = obs, Shannon = sh,
             Pielou = sh / log(obs), Faith_PD = pd, row.names = NULL)
}))
alpha_all <- cbind(alpha_all, meta[alpha_all$sample, c("Burden", "Burden_company", "lice", "lice_load", "family", "run")])
write_tsv(alpha_all, "alpha_diversity_per_fish.tsv")
INDICES <- c("Observed", "Shannon", "Pielou", "Faith_PD")
INDICES <- INDICES[vapply(INDICES, function(i) !all(is.na(alpha_all[[i]])), logical(1))]

# One index, one table, one outcome definition
alpha_test <- function(tn, on, idx) {
  o <- OUTCOMES[[on]]
  d <- alpha_all[alpha_all$table == tn & alpha_all$sample %in% o$fish, ]
  d$y <- d[[idx]]
  d$family <- factor(d$family)
  use_run <- length(unique(d$run)) > 1
  rhs <- c("x", if (use_run) "run")
  base <- data.frame(table = tn, outcome = on, index = idx, n_fish = nrow(d), n_families = nlevels(d$family))

  if (o$type == "group") {
    d$x <- factor(d[[o$var]], levels = c("Low", "High"))
    desc <- data.frame(mean_low = mean(d$y[d$x == "Low"]), sd_low = sd(d$y[d$x == "Low"]),
                       mean_high = mean(d$y[d$x == "High"]), sd_high = sd(d$y[d$x == "High"]),
                       shapiro_p_low = shapiro.test(d$y[d$x == "Low"])$p.value,
                       shapiro_p_high = shapiro.test(d$y[d$x == "High"])$p.value,
                       levene_p = car::leveneTest(y ~ x, data = d)[1, "Pr(>F)"],
                       wilcoxon_p = suppressWarnings(wilcox.test(y ~ x, data = d)$p.value))
    coef_name <- "xHigh"
  } else {
    d$x <- d[[o$var]]
    ct  <- suppressWarnings(cor.test(d$y, d$x, method = "spearman"))
    desc <- data.frame(spearman_rho = unname(ct$estimate), spearman_p = ct$p.value)
    coef_name <- "x"
  }

  # Mixed model with family as a random effect; a fixed-effects model if every
  # family has a single fish
  mixed <- nlevels(d$family) < nrow(d)
  if (mixed) {
    mm <- suppressMessages(suppressWarnings(
      lme4::lmer(as.formula(paste("y ~", paste(rhs, collapse = " + "), "+ (1 | family)")), data = d)))
    an <- car::Anova(mm, type = "II")
    est <- lme4::fixef(mm)[coef_name]; se <- sqrt(diag(as.matrix(vcov(mm))))[coef_name]
    stat <- an["x", "Chisq"]; p <- an["x", "Pr(>Chisq)"]
    singular <- lme4::isSingular(mm)
    fam_var  <- as.data.frame(lme4::VarCorr(mm))$vcov[1]
  } else {
    mm <- lm(as.formula(paste("y ~", paste(rhs, collapse = " + "))), data = d)
    an <- car::Anova(mm, type = "II")
    est <- coef(mm)[coef_name]; se <- sqrt(diag(vcov(mm)))[coef_name]
    stat <- an["x", "F value"]; p <- an["x", "Pr(>F)"]
    singular <- NA; fam_var <- NA
  }
  dh <- c(NA_real_, NA_real_)
  if (HAS_DHARMA && mixed) {
    dh <- tryCatch({
      sim <- DHARMa::simulateResiduals(mm, n = 1000, plot = FALSE)
      c(DHARMa::testUniformity(sim, plot = FALSE)$p.value, DHARMa::testDispersion(sim, plot = FALSE)$p.value)
    }, error = function(e) c(NA_real_, NA_real_))
  }
  cbind(base, desc,
        data.frame(model = if (mixed) "LMM, family random effect" else "linear model",
                   adjusted_for_run = use_run,
                   estimate = unname(est), se = unname(se),
                   statistic = stat, statistic_type = if (mixed) "Wald chi-square (type II)" else "F (type II)",
                   p_model = p, family_variance = fam_var, singular_fit = singular,
                   dharma_uniformity_p = dh[1], dharma_dispersion_p = dh[2]))
}

alpha_res <- list()
for (on in names(OUTCOMES)) {
  rows <- do.call(rbind, lapply(names(tables), function(tn) {
    r <- do.call(rbind, lapply(INDICES, function(idx) alpha_test(tn, on, idx)))
    r$p_model_BH <- p.adjust(r$p_model, method = "BH")        # across the indices of one table
    r
  }))
  num <- vapply(rows, is.numeric, logical(1))
  rows[num] <- lapply(rows[num], signif, 4)
  alpha_res[[on]] <- rows
  write_tsv(rows, paste0("alpha_tests_", on, ".tsv"))
  message("\nAlpha diversity, ", OUTCOMES[[on]]$label, ":")
  show <- intersect(c("table", "index", "n_fish", "mean_low", "mean_high", "spearman_rho", "estimate", "statistic",
                      "p_model", "p_model_BH", "wilcoxon_p", "singular_fit"), colnames(rows))
  print(rows[, show], row.names = FALSE)
}

## ---- Alpha figures -------------------------------------------------------------
index_lab <- c(Observed = "Observed ASVs", Shannon = "Shannon diversity index",
               Pielou = "Pielou's evenness", Faith_PD = "Faith's PD")
make_boxplot <- function(df, y_var, y_lab, group_var) {
  df$Burden <- factor(df[[group_var]], levels = names(cb2))
  group_means <- aggregate(df[[y_var]], list(Burden = df$Burden), mean)
  ggplot(df, aes(x = Burden, y = .data[[y_var]], fill = Burden)) +
    geom_boxplot(width = 0.5, alpha = 0.7, colour = "black", linewidth = 0.6, outlier.shape = NA) +
    geom_segment(data = group_means,
                 aes(x = as.numeric(Burden) - 0.25, xend = as.numeric(Burden) + 0.25, y = x, yend = x),
                 inherit.aes = FALSE, colour = "black", linetype = "dotted", linewidth = 0.7) +
    geom_jitter(aes(color = Burden), width = 0.15, size = 1.5, alpha = 0.6, show.legend = FALSE) +
    scale_fill_manual(values = cb2) + scale_colour_manual(values = cb2) +
    labs(x = NULL, y = y_lab) +
    theme_classic(base_size = 13) +
    theme(axis.text.x = element_text(face = "bold"), legend.position = "none")
}
for (tn in names(tables)) {
  for (on in c("gaussian", "company")) {
    o  <- OUTCOMES[[on]]
    df <- alpha_all[alpha_all$table == tn & alpha_all$sample %in% o$fish, ]
    pl <- lapply(INDICES, function(idx) make_boxplot(df, idx, index_lab[[idx]], o$var))
    fig <- wrap_plots(pl, nrow = 1) + plot_annotation(tag_levels = "A") &
      theme(plot.tag = element_text(face = "bold", size = 16))
    base <- file.path(DIV_DIR, paste0("alpha_boxplots_", tn, "_", on))
    ggsave(paste0(base, ".pdf"), fig, width = 3.5 * length(INDICES), height = 5)
    ggsave(paste0(base, ".png"), fig, width = 3.5 * length(INDICES), height = 5, dpi = 600)
  }
}

###############################################################################
# 4. BETA DIVERSITY
###############################################################################
# Bray-Curtis and Jaccard (presence/absence) on both rarefied tables; Aitchison
# distance (Euclidean distance between CLR-transformed profiles) on the
# filtered taxa
count_dist <- function(tn, method, binary = FALSE) {
  force(tn); force(method); force(binary)
  function(s) vegdist(tables[[tn]]$counts[s, , drop = FALSE], method = method, binary = binary)
}
beta_sets <- list()
for (tn in names(tables)) {
  beta_sets[[paste0("Bray-Curtis|", tn)]] <- list(metric = "Bray-Curtis", table = tn, dist = count_dist(tn, "bray"))
  beta_sets[[paste0("Jaccard|", tn)]]     <- list(metric = "Jaccard", table = tn, dist = count_dist(tn, "jaccard", binary = TRUE))
}
if (!is.null(clr)) {
  beta_sets[["Aitchison|filtered"]] <- list(metric = "Aitchison", table = "filtered (CLR, unrarefied)",
    dist = function(s) dist(clr[s, , drop = FALSE]))
}

# PERMANOVA for one distance and one outcome: the outcome alone, and the
# outcome after accounting for sequencing run and depth (marginal test)
permanova_rows <- function(bs, on) {
  o <- OUTCOMES[[on]]
  s <- o$fish
  d <- bs$dist(s)
  m <- meta[s, ]
  m$x <- if (o$type == "group") factor(m[[o$var]], levels = c("Low", "High")) else m[[o$var]]
  covs <- c(if (length(unique(m$run)) > 1) "run", if (!anyNA(m$log_depth)) "log_depth")

  set.seed(PERM_SEED)
  a1 <- adonis2(d ~ x, data = m, permutations = N_PERM, by = "margin")
  out <- data.frame(metric = bs$metric, table = bs$table, outcome = on, n_fish = length(s),
                    model = "burden only", Df = a1["x", "Df"], pseudo_F = a1["x", "F"],
                    R2 = a1["x", "R2"], p = a1["x", "Pr(>F)"])
  if (length(covs)) {
    set.seed(PERM_SEED)
    a2 <- adonis2(as.formula(paste("d ~", paste(c(covs, "x"), collapse = " + "))), data = m,
                  permutations = N_PERM, by = "margin")
    out <- rbind(out, data.frame(metric = bs$metric, table = bs$table, outcome = on, n_fish = length(s),
                                 model = paste("burden adjusted for", paste(covs, collapse = " + ")),
                                 Df = a2["x", "Df"], pseudo_F = a2["x", "F"], R2 = a2["x", "R2"], p = a2["x", "Pr(>F)"]))
  }
  # Homogeneity of multivariate dispersion (groups only)
  out$dispersion_F <- NA_real_; out$dispersion_p <- NA_real_
  if (o$type == "group") {
    set.seed(PERM_SEED)
    pt <- permutest(betadisper(d, m$x), permutations = N_PERM)
    out$dispersion_F <- pt$tab[1, "F"]; out$dispersion_p <- pt$tab[1, "Pr(>F)"]
  }
  out$permutations <- N_PERM
  out
}

perm_tab <- do.call(rbind, lapply(names(OUTCOMES), function(on) {
  do.call(rbind, lapply(beta_sets, permanova_rows, on = on))
}))
rownames(perm_tab) <- NULL
perm_tab[c("pseudo_F", "R2", "dispersion_F")] <- lapply(perm_tab[c("pseudo_F", "R2", "dispersion_F")], round, 4)
write_tsv(perm_tab, "beta_permanova.tsv")
message("\nPERMANOVA (", N_PERM, " permutations) and homogeneity of dispersion:")
print(perm_tab[, c("metric", "table", "outcome", "n_fish", "model", "Df", "pseudo_F", "R2", "p", "dispersion_p")], row.names = FALSE)

## ---- Ordination and dispersion figures (Gaussian and company classes) ----------
ord_plot <- function(pts, group, xlab, ylab, note) {
  df <- data.frame(A1 = pts[, 1], A2 = pts[, 2], Burden = factor(group, levels = names(cb2)))
  ggplot(df, aes(x = A1, y = A2, color = Burden)) +
    geom_point(size = 3, alpha = 0.85) +
    stat_ellipse(aes(group = Burden, fill = Burden), type = "t", level = 0.95, geom = "polygon", alpha = 0.2, linewidth = 0.5) +
    scale_color_manual(values = cb2, name = "Burden") + scale_fill_manual(values = cb2, name = "Burden") +
    labs(x = xlab, y = ylab) +
    theme_classic(base_size = 13) +
    annotate("text", x = -Inf, y = Inf, label = note, hjust = -0.05, vjust = 1.3, size = 3.5)
}
for (on in c("gaussian", "company")) {
  o <- OUTCOMES[[on]]
  grp <- meta[o$fish, o$var]
  ords <- list(); disps <- list(); stress <- NULL
  for (bn in names(beta_sets)) {
    bs <- beta_sets[[bn]]
    if (bs$table == "all_ASVs") next                     # figures: the filtered table
    d <- bs$dist(o$fish)
    if (bs$metric == "Aitchison") {
      pc <- cmdscale(d, k = 2, eig = TRUE)
      ve <- round(100 * pc$eig[1:2] / sum(pc$eig[pc$eig > 0]), 1)
      ords[[bs$metric]] <- ord_plot(pc$points, grp, paste0("PC1 (", ve[1], "%)"), paste0("PC2 (", ve[2], "%)"), "Aitchison")
    } else {
      set.seed(PERM_SEED)
      nm <- metaMDS(d, k = 2, trymax = 100, trace = FALSE)
      stress <- rbind(stress, data.frame(outcome = on, metric = bs$metric, stress = round(nm$stress, 4), converged = nm$converged))
      ords[[bs$metric]] <- ord_plot(nm$points, grp, "NMDS1", "NMDS2", sprintf("%s, stress = %.3f", bs$metric, nm$stress))
    }
    bd <- betadisper(d, factor(grp, levels = names(cb2)))
    disps[[bs$metric]] <- ggplot(data.frame(Distance = bd$distances, Burden = factor(grp, levels = names(cb2))),
                                 aes(x = Burden, y = Distance, fill = Burden)) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.5, colour = "black") +
      geom_jitter(width = 0.15, size = 1.6, alpha = 0.6, aes(color = Burden)) +
      scale_fill_manual(values = cb2) + scale_colour_manual(values = cb2) +
      labs(y = "Distance to centroid", x = NULL, title = bs$metric) +
      theme_classic(base_size = 13) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "none")
  }
  if (!is.null(stress)) write_tsv(stress, paste0("nmds_stress_", on, ".tsv"))
  fig_o <- wrap_plots(ords, nrow = 1, guides = "collect") + plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 16))
  fig_d <- wrap_plots(disps, nrow = 1) + plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 16))
  w <- 5 * length(ords)
  ggsave(file.path(DIV_DIR, paste0("beta_ordination_", on, ".pdf")), fig_o, width = w, height = 5)
  ggsave(file.path(DIV_DIR, paste0("beta_ordination_", on, ".png")), fig_o, width = w, height = 5, dpi = 600)
  ggsave(file.path(DIV_DIR, paste0("beta_dispersion_", on, ".pdf")), fig_d, width = 4 * length(disps), height = 5)
  ggsave(file.path(DIV_DIR, paste0("beta_dispersion_", on, ".png")), fig_d, width = 4 * length(disps), height = 5, dpi = 600)
}

###############################################################################
# 5. RECORD
###############################################################################
writeLines(c(
  paste("date:", format(Sys.time())),
  paste("fish:", length(fish), "; Gaussian cut-offs:", round(LOW_CUT, 4), "and", round(HIGH_CUT, 4), "lice per cm"),
  paste("company categories: Low <", COMPANY_LOW, ", High >", COMPANY_HIGH, "lice"),
  paste("tables:", paste(names(tables), collapse = ", "), "; all-ASV table rarefied to", RAREFY_DEPTH, "reads, seed", RAREFY_SEED),
  "alpha: Observed ASVs, Shannon, Pielou (Shannon / ln Observed), Faith's PD (root included)",
  "alpha tests: linear mixed model, index ~ burden + sequencing run + (1 | family); type II Wald chi-square;",
  "  Benjamini-Hochberg adjustment across the indices of one table; Wilcoxon test as a non-parametric check",
  "beta: Bray-Curtis, Jaccard (presence/absence), Aitchison (Euclidean on CLR, pseudocount as in the QC script)",
  paste("PERMANOVA: adonis2, marginal tests,", N_PERM, "permutations, seed", PERM_SEED,
        "; dispersion: betadisper + permutest"),
  capture.output(sessionInfo())
), file.path(DIV_DIR, "diversity_run_record.txt"))
message("\nDiversity: tables and figures in ", DIV_DIR)
