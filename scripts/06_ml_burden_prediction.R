#!/usr/bin/env Rscript

###############################################################################
# Skin microbiome & Caligus: prediction of sea lice burden from the skin
# microbiome (machine learning with and without TaxaHFE feature selection)
# Needs the outputs of 02_qc_decontam.R.
#
# How this differs from the original analysis (reviewer comment in brackets):
#   - TaxaHFE is run on the training part of every fold only; the test fish
#     never influence which features are selected              [R1-1]
#   - centring/scaling use training-fold values only            [R1-1, R1-8]
#   - the only tuned settings (glmnet penalties) are tuned inside the training
#     fold; all other settings are fixed in advance and listed  [R1-16, R2 minor]
#   - folds are stratified by burden class and keep each family together
#                                                               [R1-6, R1-16, R2-10]
#   - AUC is computed within each test fold with a fixed direction and
#     averaged; its 95% CI and the comparisons between feature sets use the
#     corrected resampled t-test for repeated cross-validation (Nadeau &
#     Bengio 2003); sensitivity, specificity, Brier score and calibration are
#     reported                                                  [R1-2, R1-16]
#   - baseline models from host and technical variables, and combined models
#                                                               [R1-7, R2-11]
#   - how often each HFE feature is selected across folds       [R1-13]
#   - sensitivity analyses: without the second sequencing run, unrarefied
#     proportions, CLR, and folds that ignore family            [R1-8, R1-9]
#   - alternative definitions of the outcome: the Gaussian-derived extremes
#     (primary), the company categories on raw counts, all fish split at the
#     median, and lice load as a continuous trait (per cm and raw count),
#     each of which can also be run with TaxaHFE                [R1-4, R2-4, R2-5]
#   - the old procedure (HFE on all fish, then cross-validation) is kept as
#     "hfe_leaky", for comparison only, to show how optimistic it was
#
# Two-step use, as for SourceTracker: run this script; it writes the TaxaHFE
# inputs and ml/run_hfe.sh. Run that shell script (Docker), then run this
# script again. Models that do not need TaxaHFE are fitted on the first run.
###############################################################################

## --------------------------- User Parameters ---------------------------------
PROJECT_DIR <- getwd()
RDS_DIR     <- file.path(PROJECT_DIR)
ML_DIR      <- file.path(PROJECT_DIR, "ml")

# Cross-validation
N_FOLDS   <- 5
N_REPEATS <- 20        # each repeat needs N_FOLDS TaxaHFE runs per HFE variant
SEED      <- 1000

# Metadata columns
FAMILY_COL <- "Fish_family"
HOST_VARS  <- c("Final_Length", "Final_Weight", "Sex")

# Models (the eight of the original analysis) and the one used for headline
# comparisons and plots
MODELS        <- c("rf", "xgboost", "svm", "lasso", "ridge", "enet", "knn", "lda")
PRIMARY_MODEL <- "rf"

# TaxaHFE (Docker), command line of version 2.4: the image starts TaxaHFE
# itself, reads the inputs from /data and writes <output dir>/taxahfe_sf.csv.
# Pin a specific image version instead of "latest" for the paper, and report
# it. On Apple-silicon Macs you may need
# HFE_DOCKER_ARGS <- c("--platform", "linux/amd64")
HFE_DOCKER_IMAGE <- "aoliver44/taxa_hfe:latest"
HFE_DOCKER_ARGS  <- character(0)
HFE_LOWEST_LEVEL <- 3
HFE_NCORES       <- 4
HFE_SEED         <- 123
RUN_HFE_FROM_R   <- FALSE   # TRUE runs ml/run_hfe.sh from this script

# Company-defined burden categories on the raw lice count:
# Low < COMPANY_LOW, High > COMPANY_HIGH, Intermediate in between
COMPANY_LOW  <- 15
COMPANY_HIGH <- 30

# Analyses to run. Only those listed in HFE_VARIANTS use TaxaHFE (each costs
# N_FOLDS x N_REPEATS + 1 Docker runs). "clr" cannot use it. To answer the
# outcome-definition comment in full, add "company", "all_fish_median" and
# "regression" to HFE_VARIANTS.
RUN_VARIANTS <- c("primary", "company", "all_fish_median", "regression", "regression_count",
                  "no_run2", "proportions", "clr", "ungrouped")
HFE_VARIANTS <- c("primary", "company", "all_fish_median", "regression")

## ------------------------------ Libraries ------------------------------------
suppressPackageStartupMessages({
  library(phyloseq)
  library(VGAM)
  library(randomForest)
  library(e1071)
  library(glmnet)
  library(MASS)
  library(class)
  library(pROC)
  library(ggplot2)
})
HAS_XGB <- requireNamespace("xgboost", quietly = TRUE)
if ("xgboost" %in% MODELS && !HAS_XGB) {
  message("Package xgboost not installed: that model is skipped.")
  MODELS <- setdiff(MODELS, "xgboost")
}
for (d in c(ML_DIR, file.path(ML_DIR, "hfe"))) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

## ------------------------------- Helpers -------------------------------------
otu_mat <- function(x) {            # counts with ASVs in rows
  m <- as(otu_table(x), "matrix")
  if (!taxa_are_rows(x)) m <- t(m)
  m
}
write_tsv <- function(x, name) {
  write.table(x, file.path(ML_DIR, name), sep = "\t", quote = FALSE, row.names = FALSE)
}
to_num <- function(v) suppressWarnings(as.numeric(as.character(v)))
load_ps <- function(file, required = FALSE) {
  f <- file.path(RDS_DIR, file)
  if (!file.exists(f)) {
    if (required) stop("File not found: ", f)
    return(NULL)
  }
  readRDS(f)
}
clean_name <- function(x) gsub("^_+|_+$", "", gsub("[^a-z0-9]+", "_", tolower(x)))
alnum_only <- function(x) gsub("[^a-z0-9]", "", tolower(x))
HFE_OUT    <- file.path("out", "taxahfe_sf.csv")       # TaxaHFE result inside each input folder

###############################################################################
# 1. DATA AND PHENOTYPE
###############################################################################
ps_rare <- load_ps("Skin_rare_SILVA.rds", required = TRUE)   # primary dataset
ps_skin <- load_ps("Skin_ps_SILVA.rds")                      # unrarefied counts (depth)
ps_prop <- load_ps("Skin_proportions_SILVA.rds")
ps_clr  <- load_ps("Skin_clr_SILVA.rds")

# One metadata table for every fish available in any dataset
ps_meta <- if (!is.null(ps_skin)) ps_skin else ps_rare
sd_all  <- as(sample_data(ps_meta), "data.frame")
meta <- data.frame(
  sample    = sample_names(ps_meta),
  lice      = to_num(sd_all$Total_caligus),
  length    = to_num(sd_all$Final_Length),
  weight    = to_num(sd_all$Final_Weight),
  sex       = as.character(sd_all$Sex),
  family    = if (FAMILY_COL %in% colnames(sd_all)) as.character(sd_all[[FAMILY_COL]]) else NA_character_,
  run       = if ("seq_batch" %in% colnames(sd_all)) as.character(sd_all$seq_batch) else "run1",
  log_depth = if (!is.null(ps_skin)) log10(sample_sums(ps_skin)[sample_names(ps_meta)]) else NA_real_,
  row.names = sample_names(ps_meta), stringsAsFactors = FALSE
)
meta$lice_load <- meta$lice / meta$length       # lice per cm, as in the original analysis
if (all(is.na(meta$family))) {
  warning("Column '", FAMILY_COL, "' not found: folds cannot keep families together.", immediate. = TRUE)
}

## ---- Burden classes from a two-Gaussian mixture (as in the original) ---------
# Fitted on the fish of the primary (rarefied) dataset; the same two cut-offs
# are then applied to every dataset.
fit_fish <- sample_names(ps_rare)
tc <- meta[fit_fish, "lice_load"]
fit  <- vglm(tc ~ 1, mix2normal(eq.sd = FALSE), trace = FALSE)
pars <- as.vector(coef(fit))
mix  <- list(w = plogis(pars[1]), m1 = pars[2], sd1 = exp(pars[3]), m2 = pars[4], sd2 = exp(pars[5]))
if (mix$m1 > mix$m2) mix <- list(w = 1 - mix$w, m1 = mix$m2, sd1 = mix$sd2, m2 = mix$m1, sd2 = mix$sd1)
LOW_CUT  <- mix$m1 + mix$sd1     # Low  : lice load <= mean 1 + 1 SD
HIGH_CUT <- mix$m2 - mix$sd2     # High : lice load >= mean 2 - 1 SD
if (LOW_CUT >= HIGH_CUT) stop("The mixture fit gives overlapping cut-offs (", round(LOW_CUT, 3), " >= ", round(HIGH_CUT, 3), ").")
meta$burden <- ifelse(meta$lice_load <= LOW_CUT, "Low", ifelse(meta$lice_load >= HIGH_CUT, "High", "Moderate"))
MEDIAN_CUT <- median(tc)
meta$burden_median <- ifelse(meta$lice_load > MEDIAN_CUT, "High", "Low")

thr <- data.frame(
  item  = c("fish used for the mixture fit", "mixture mean 1", "mixture SD 1", "mixture mean 2", "mixture SD 2",
            "mixing weight of component 1", "Low cut-off (lice per cm, <=)", "High cut-off (lice per cm, >=)",
            "median lice per cm (median-split analysis)",
            "Low fish", "Moderate fish (excluded from the main analysis)", "High fish"),
  value = c(length(tc), round(c(mix$m1, mix$sd1, mix$m2, mix$sd2, mix$w, LOW_CUT, HIGH_CUT, MEDIAN_CUT), 4),
            sum(meta[fit_fish, "burden"] == "Low"), sum(meta[fit_fish, "burden"] == "Moderate"),
            sum(meta[fit_fish, "burden"] == "High")))
write_tsv(thr, "burden_thresholds.tsv")
print(thr, row.names = FALSE)

# Company categories (raw counts), and how the two definitions agree
meta$burden_company <- ifelse(meta$lice < COMPANY_LOW, "Low", ifelse(meta$lice > COMPANY_HIGH, "High", "Intermediate"))
xt <- table(factor(meta[fit_fish, "burden"], levels = c("Low", "Moderate", "High")),
            factor(meta[fit_fish, "burden_company"], levels = c("Low", "Intermediate", "High")))
xt_df <- data.frame(gaussian_class = rownames(xt), as.data.frame.matrix(xt), row.names = NULL)
colnames(xt_df)[-1] <- paste0("company_", colnames(xt))
write_tsv(xt_df, "burden_definitions_crosstab.tsv")
message("Gaussian-derived classes (rows) against company categories (columns: Low < ", COMPANY_LOW,
        ", High > ", COMPANY_HIGH, " lice):")
print(xt_df, row.names = FALSE)

## ---- Families, and host/technical variables against burden class -------------
prim <- meta[fit_fish, ]
prim_hl <- prim[prim$burden != "Moderate", ]
if (!all(is.na(prim_hl$family))) {
  fam <- as.data.frame.matrix(table(family = prim_hl$family, burden = prim_hl$burden))
  fam <- data.frame(family = rownames(fam), fam, n_fish = rowSums(fam))
  write_tsv(fam[order(-fam$n_fish), ], "family_representation.tsv")
  message("Families among the High/Low fish: ", nrow(fam), "; families with more than one fish: ",
          sum(fam$n_fish > 1), "; largest family: ", max(fam$n_fish), " fish")
}
assoc <- do.call(rbind, lapply(c("length", "weight", "log_depth", "lice"), function(v) {
  x <- prim_hl[[v]]
  if (all(is.na(x))) return(NULL)
  data.frame(variable = v, test = "Wilcoxon",
             low = round(median(x[prim_hl$burden == "Low"], na.rm = TRUE), 3),
             high = round(median(x[prim_hl$burden == "High"], na.rm = TRUE), 3),
             p = signif(suppressWarnings(wilcox.test(x ~ prim_hl$burden))$p.value, 3))
}))
for (v in c("sex", "run")) {
  tb <- table(prim_hl[[v]], prim_hl$burden)
  if (nrow(tb) > 1) {
    assoc <- rbind(assoc, data.frame(variable = v, test = "Fisher",
                                     low = paste(rownames(tb), tb[, "Low"], collapse = "; "),
                                     high = paste(rownames(tb), tb[, "High"], collapse = "; "),
                                     p = signif(fisher.test(tb)$p.value, 3)))
  }
}
write_tsv(assoc, "host_and_technical_variables_vs_burden.tsv")
message("Host and technical variables against burden class (medians or counts in Low and High):")
print(assoc, row.names = FALSE)

###############################################################################
# 2. FEATURES
###############################################################################
# Taxonomy paths with a placeholder for missing ranks, so that every ASV has
# all seven levels (the original pasted prefixes onto the non-missing ranks
# only, which shifted the ranks of ASVs with a missing genus).
tax <- as(tax_table(ps_rare), "matrix")
RANKS    <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
PREFIXES <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__")
tax_filled <- t(apply(tax[, RANKS, drop = FALSE], 1, function(r) {
  last <- "root"
  for (i in seq_along(r)) {
    if (is.na(r[i]) || r[i] == "") r[i] <- paste0("unclassified_", last) else last <- r[i]
  }
  r
}))
tax_filled <- gsub("[|\t ]+", "_", tax_filled)
tax_full   <- cbind(tax_filled, Species = rownames(tax))
paths <- sapply(seq_len(7), function(k) {
  apply(tax_full[, seq_len(k), drop = FALSE], 1, function(r) paste0(PREFIXES[seq_len(k)], r, collapse = "|"))
})
rownames(paths) <- rownames(tax)

# Relative abundance (samples x ASVs) from a count-like phyloseq object
rel_abund <- function(ps) {
  m <- t(otu_mat(ps))[, rownames(tax), drop = FALSE]
  m / rowSums(m)
}
# Summed relative abundance of every clade at every rank (samples x clades)
clade_matrix <- function(rel) {
  out <- lapply(seq_len(7), function(k) {
    cm <- t(rowsum(t(rel), group = paths[colnames(rel), k]))
    cm
  })
  mat  <- do.call(cbind, out)
  info <- data.frame(clade = colnames(mat),
                     rank  = rep(c(RANKS, "ASV"), times = vapply(out, ncol, integer(1))),
                     name  = sub("^.*\\|", "", colnames(mat)), stringsAsFactors = FALSE)
  info$name <- sub("^[a-z]__", "", info$name)
  list(mat = mat, info = info)
}

###############################################################################
# 3. FOLDS: stratified by class, families kept together
###############################################################################
make_folds <- function(strat, groups, k, repeats, seed) {
  set.seed(seed)
  strat <- factor(strat)
  n <- length(strat)
  out <- matrix(NA_integer_, n, repeats)
  for (r in seq_len(repeats)) {
    g_ids <- sample(unique(groups))
    g_ids <- g_ids[order(-as.vector(table(groups)[g_ids]))]     # large families first
    counts <- matrix(0, k, nlevels(strat))
    for (g in g_ids) {
      idx <- which(groups == g)
      cnt <- as.vector(table(strat[idx]))
      # fold with the fewest fish of the classes this family carries; ties: smallest fold, then random
      score <- as.vector(counts %*% (cnt > 0)) + rowSums(counts) / (10 * n) + runif(k) / (1000 * n)
      f <- which.min(score)
      out[idx, r] <- f
      counts[f, ] <- counts[f, ] + cnt
    }
  }
  out
}

###############################################################################
# 4. MODELS
###############################################################################
# Settings fixed in advance (not tuned): random forest 500 trees, default mtry;
# xgboost 10 rounds, default eta and depth; SVM radial kernel, cost 1, default
# gamma; k-NN k = 5; LDA with tol = 0. Tuned inside the training fold by an
# inner 5-fold cross-validation: the glmnet penalty (lasso, ridge, elastic
# net) and the elastic-net mixing parameter (0.25, 0.5, 0.75).
prep_xy <- function(xtr, xte) {
  for (j in seq_len(ncol(xtr))) {                       # missing values: training median
    med <- median(xtr[, j], na.rm = TRUE)
    xtr[is.na(xtr[, j]), j] <- med
    xte[is.na(xte[, j]), j] <- med
  }
  mu <- colMeans(xtr)
  s  <- apply(xtr, 2, sd)
  ok <- is.finite(s) & s > 0
  if (!any(ok)) return(NULL)
  list(xtr = scale(xtr[, ok, drop = FALSE], center = mu[ok], scale = s[ok]),
       xte = scale(xte[, ok, drop = FALSE], center = mu[ok], scale = s[ok]))
}

glmnet_prob <- function(xtr, ytr, xte, alphas, family) {
  if (ncol(xtr) < 2) {                                  # glmnet needs two columns
    d  <- data.frame(y = ytr, x = xtr[, 1])
    fm <- glm(y ~ x, data = d, family = if (family == "binomial") binomial() else gaussian())
    return(as.vector(predict(fm, data.frame(x = xte[, 1]), type = "response")))
  }
  nf <- min(5, if (family == "binomial") min(table(ytr)) else length(ytr))
  foldid <- sample(rep(seq_len(nf), length.out = length(ytr)))
  fits <- lapply(alphas, function(a) cv.glmnet(xtr, ytr, family = family, alpha = a, foldid = foldid))
  best <- fits[[which.min(vapply(fits, function(f) min(f$cvm), numeric(1)))]]
  as.vector(predict(best, xte, s = "lambda.min", type = "response"))
}

fit_predict <- function(model, xtr, ytr, xte) {
  yf <- factor(ytr, levels = c(0, 1))
  switch(model,
    rf = predict(randomForest(xtr, y = yf), xte, type = "prob")[, "1"],
    xgboost = {
      b <- xgboost::xgb.train(params = list(objective = "binary:logistic", nthread = 1),
                              data = xgboost::xgb.DMatrix(xtr, label = ytr), nrounds = 10, verbose = 0)
      predict(b, xgboost::xgb.DMatrix(xte))
    },
    svm = {
      m <- svm(xtr, yf, probability = TRUE, scale = FALSE)
      attr(predict(m, xte, probability = TRUE), "probabilities")[, "1"]     # by name: P(High)
    },
    lasso = glmnet_prob(xtr, ytr, xte, 1, "binomial"),
    ridge = glmnet_prob(xtr, ytr, xte, 0, "binomial"),
    enet  = glmnet_prob(xtr, ytr, xte, c(0.25, 0.5, 0.75), "binomial"),
    knn = {
      kk  <- min(5, nrow(xtr) - 1)
      res <- class::knn(xtr, xte, yf, k = kk, prob = TRUE)
      p   <- attr(res, "prob")
      ifelse(res == "1", p, 1 - p)
    },
    lda = predict(MASS::lda(xtr, grouping = yf, tol = 0), xte)$posterior[, "1"],
    stop("Unknown model: ", model))
}
fit_predict_reg <- function(model, xtr, ytr, xte) {
  switch(model,
    rf   = predict(randomForest(xtr, y = ytr), xte),
    enet = glmnet_prob(xtr, ytr, xte, c(0.25, 0.5, 0.75), "gaussian"),
    stop("Unknown model: ", model))
}

###############################################################################
# 5. TAXAHFE INSIDE THE FOLDS
###############################################################################
hfe_cmds <- character(0)

# Writes the TaxaHFE input for one set of training fish and queues the command
hfe_prepare <- function(dir, rel_train, labels, numeric_label = FALSE) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  otu <- data.frame(clade_name = paths[colnames(rel_train), 7], t(rel_train), check.names = FALSE)
  write.table(otu, file.path(dir, "otu.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(data.frame(Sample_name = rownames(rel_train), Burden = labels),
              file.path(dir, "meta.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  if (!file.exists(file.path(dir, HFE_OUT))) {
    hfe_cmds <<- c(hfe_cmds, paste(
      "docker run --rm", paste(HFE_DOCKER_ARGS, collapse = " "), paste0("--cpus=", HFE_NCORES),
      "-v", paste0(shQuote(normalizePath(dir)), ":/data"), HFE_DOCKER_IMAGE,
      "meta.txt otu.txt -o out",
      "--subject_identifier Sample_name --label Burden",
      if (numeric_label) "--feature_type numeric",
      "--lowest_level", HFE_LOWEST_LEVEL, "--ncores", HFE_NCORES, "--seed", HFE_SEED))
  }
}

# Reads one TaxaHFE result and returns the clades it selected. Each output
# column is matched to a clade by its values in the training fish, so the
# result does not depend on how TaxaHFE names its columns. The selected clades
# can then be computed for the test fish from their own abundances.
hfe_read <- function(dir, clades_train, info) {
  res <- read.csv(file.path(dir, HFE_OUT), check.names = FALSE, stringsAsFactors = FALSE)
  ids <- as.character(res[[1]])
  tr  <- rownames(clades_train)
  pos <- match(tr, ids)
  if (anyNA(pos)) pos <- match(clean_name(tr), clean_name(ids))
  if (anyNA(pos)) {
    if (nrow(res) != length(tr)) stop("Cannot align TaxaHFE output with the training fish in ", dir)
    warning("TaxaHFE sample IDs not recognised in ", dir, "; rows are assumed to be in input order.")
    pos <- seq_along(tr)
  }
  feat <- res[pos, -c(1, 2), drop = FALSE]
  sel <- character(0); log <- NULL
  for (j in seq_len(ncol(feat))) {
    v <- to_num(feat[[j]])
    if (anyNA(v) || sd(v) == 0) next
    r <- suppressWarnings(as.vector(cor(v, clades_train)))
    r[is.na(r)] <- -1
    cand <- which(r >= max(r) - 1e-9)
    hint <- which(alnum_only(info$clade[cand]) == alnum_only(colnames(feat)[j]))        # TaxaHFE names columns by the full path
    best <- if (length(hint)) cand[hint[1]] else cand[1]          # ties: name match, else highest rank
    ok   <- r[best] >= 0.999
    log  <- rbind(log, data.frame(hfe_column = colnames(feat)[j], clade = info$clade[best], rank = info$rank[best],
                                  name = info$name[best], correlation = round(r[best], 5), matched = ok))
    if (ok) sel <- c(sel, info$clade[best])
  }
  list(selected = unique(sel), log = log)
}

###############################################################################
# 6. ONE ANALYSIS (variant)
###############################################################################
VARIANTS <- list(
  primary         = list(data = "rarefied",    label = "extremes", drop_run2 = FALSE, grouped = TRUE,
                         note = "High vs Low, rarefied, families kept together"),
  no_run2         = list(data = "rarefied",    label = "extremes", drop_run2 = TRUE,  grouped = TRUE,
                         note = "as primary, without the fish of the second sequencing run"),
  all_fish_median = list(data = "rarefied",    label = "median",   drop_run2 = FALSE, grouped = TRUE,
                         note = "all fish, split at the median lice load (includes intermediate fish)"),
  proportions     = list(data = "proportions", label = "extremes", drop_run2 = FALSE, grouped = TRUE,
                         note = "as primary, unrarefied relative abundance (keeps low-depth fish)"),
  clr             = list(data = "clr",         label = "extremes", drop_run2 = FALSE, grouped = TRUE,
                         note = "as primary, CLR-transformed abundance"),
  ungrouped       = list(data = "rarefied",    label = "extremes", drop_run2 = FALSE, grouped = FALSE,
                         note = "as primary, folds ignore family (shows the effect of grouping)"),
  regression      = list(data = "rarefied",    label = "continuous", drop_run2 = FALSE, grouped = TRUE,
                         note = "all fish, lice per cm as a continuous trait"),
  company         = list(data = "rarefied",    label = "company",  drop_run2 = FALSE, grouped = TRUE,
                         note = "company categories on raw counts, High vs Low, intermediate fish excluded"),
  regression_count = list(data = "rarefied",   label = "continuous_count", drop_run2 = FALSE, grouped = TRUE,
                         note = "all fish, raw lice count as a continuous trait")
)

auc_fixed <- function(y, p) {
  if (anyNA(p) || length(unique(y)) < 2) return(NA_real_)
  as.numeric(pROC::auc(pROC::roc(y, p, levels = c(0, 1), direction = "<", quiet = TRUE)))
}
# Mean of fold-wise values with a 95% CI from the corrected resampled t-test
# (Nadeau & Bengio 2003). Test folds of repeated cross-validation overlap, so
# the ordinary standard error would be too small; the correction adds
# n_test / n_train = 1 / (k - 1) to the variance factor.
nb_summary <- function(x, k) {
  x <- x[!is.na(x)]
  J <- length(x)
  if (J < 2) return(c(mean = mean(x), low = NA, high = NA, se = NA, df = NA))
  se <- sqrt((1 / J + 1 / (k - 1)) * var(x))
  tq <- qt(0.975, J - 1)
  c(mean = mean(x), low = mean(x) - tq * se, high = mean(x) + tq * se, se = se, df = J - 1)
}
# Performance is computed inside each test fold and then averaged. Pooling the
# predictions of all folds before computing AUC is biased downwards: a fold
# with more Low fish in the test part has fewer in the training part, which
# shifts all its predictions towards High.

run_variant <- function(vname) {
  v  <- VARIANTS[[vname]]
  ps <- switch(v$data, rarefied = ps_rare, proportions = ps_prop, clr = ps_clr)
  if (is.null(ps)) { message("Variant ", vname, ": input object not found, skipped."); return(NULL) }
  is_reg  <- v$label %in% c("continuous", "continuous_count")
  use_hfe <- vname %in% HFE_VARIANTS && v$data != "clr"

  # Fish and outcome
  m <- meta[sample_names(ps), ]
  if (v$drop_run2) m <- m[m$run != "run2", ]
  if (v$label == "extremes") m <- m[m$burden != "Moderate", ]
  if (v$label == "company")  m <- m[m$burden_company != "Intermediate", ]
  y <- switch(v$label, extremes = as.numeric(m$burden == "High"),
              company = as.numeric(m$burden_company == "High"),
              median = as.numeric(m$burden_median == "High"),
              continuous = m$lice_load, continuous_count = m$lice)
  fish <- m$sample
  n <- length(fish)
  if (!is_reg && min(table(factor(y, levels = c(0, 1)))) < 2 * N_FOLDS) {
    message("Variant ", vname, ": too few fish in one class (", paste(table(factor(y, levels = c(0, 1))), collapse = " vs "),
            ") for ", N_FOLDS, "-fold cross-validation, skipped.")
    return(NULL)
  }

  # Features
  if (v$data == "clr") {
    x_asv <- t(otu_mat(ps))[fish, , drop = FALSE]
    clades <- NULL
  } else {
    rel    <- rel_abund(ps)[fish, , drop = FALSE]
    x_asv  <- rel
    clades <- clade_matrix(rel)
  }
  x_host <- cbind(length = m$length, weight = m$weight)
  if (length(unique(m$sex)) > 1) {
    sx <- model.matrix(~ sex, data = data.frame(sex = factor(m$sex)))[, -1, drop = FALSE]
    x_host <- cbind(x_host, sx)
  }
  x_tech <- cbind(log_depth = m$log_depth, run2 = as.numeric(m$run == "run2"))
  x_tech <- x_tech[, colSums(!is.na(x_tech)) > 0, drop = FALSE]
  rownames(x_host) <- rownames(x_tech) <- fish

  # Folds: created once and stored, so that TaxaHFE results stay valid between runs
  vdir <- file.path(ML_DIR, "hfe", vname)
  if (!dir.exists(vdir)) dir.create(vdir, recursive = TRUE)
  groups <- if (v$grouped && !all(is.na(m$family))) ifelse(is.na(m$family), paste0("single_", fish), m$family) else fish
  strat  <- if (is_reg) cut(y, quantile(y, 0:4 / 4), include.lowest = TRUE) else y
  ffile  <- file.path(vdir, "folds.rds")
  key    <- list(fish = fish, y = y, k = N_FOLDS, repeats = N_REPEATS, grouped = v$grouped)
  folds  <- NULL
  if (file.exists(ffile)) {
    old <- readRDS(ffile)
    if (identical(old$key, key)) folds <- old$folds else if (use_hfe) {
      stop("Variant ", vname, ": the fish, classes or cross-validation settings changed since the TaxaHFE inputs were written.\n",
           "Delete the folder ", vdir, " and run TaxaHFE again.")
    }
  }
  if (is.null(folds)) {
    folds <- make_folds(strat, groups, N_FOLDS, N_REPEATS, SEED)
    saveRDS(list(key = key, folds = folds), ffile)
  }
  if (!is_reg) {
    fold_tab <- do.call(rbind, lapply(seq_len(N_REPEATS), function(r) {
      tb <- table(factor(folds[, r], levels = seq_len(N_FOLDS)), factor(y, levels = c(0, 1)))
      data.frame(variant = vname, rep = r, fold = seq_len(N_FOLDS), n_low = tb[, 1], n_high = tb[, 2])
    }))
  } else fold_tab <- NULL

  # TaxaHFE: inputs for every training set and for all fish (the leaky version)
  hfe_sel <- NULL; hfe_full <- NULL; hfe_log <- NULL
  if (use_hfe) {
    lab <- if (is_reg) y else ifelse(y == 1, "High", "Low")
    jobs <- c(list(list(dir = file.path(vdir, "all_fish"), train = seq_len(n))),
              unlist(lapply(seq_len(N_REPEATS), function(r) lapply(seq_len(N_FOLDS), function(f) {
                list(dir = file.path(vdir, sprintf("rep%02d_fold%d", r, f)), train = which(folds[, r] != f), r = r, f = f)
              })), recursive = FALSE))
    for (j in jobs) hfe_prepare(j$dir, rel[j$train, , drop = FALSE], lab[j$train], numeric_label = is_reg)
    have <- vapply(jobs, function(j) file.exists(file.path(j$dir, HFE_OUT)), logical(1))
    if (all(have)) {
      hfe_sel <- matrix(list(), N_REPEATS, N_FOLDS)
      for (j in jobs) {
        res <- hfe_read(j$dir, clades$mat[j$train, , drop = FALSE], clades$info)
        if (!is.null(res$log)) hfe_log <- rbind(hfe_log, data.frame(variant = vname, run = basename(j$dir), res$log))
        if (is.null(j$r)) hfe_full <- res$selected else hfe_sel[[j$r, j$f]] <- res$selected
      }
    } else {
      message("Variant ", vname, ": TaxaHFE results missing for ", sum(!have), " of ", length(have),
              " training sets; HFE feature sets are left out of this run.")
      use_hfe <- FALSE
    }
  }

  # Feature sets
  fsets <- c("host", if (ncol(x_tech) > 0) "technical", "asv", "host_asv")
  if (use_hfe) fsets <- c(fsets, "hfe", "host_hfe", "hfe_leaky")
  models <- if (is_reg) c("rf", "enet") else MODELS
  get_x <- function(fs, r, f) {
    switch(fs,
      host      = x_host,
      technical = x_tech,
      asv       = x_asv,
      host_asv  = cbind(x_host, x_asv),
      hfe       = clades$mat[, hfe_sel[[r, f]], drop = FALSE],
      host_hfe  = cbind(x_host, clades$mat[, hfe_sel[[r, f]], drop = FALSE]),
      hfe_leaky = clades$mat[, hfe_full, drop = FALSE])
  }

  # Cross-validation
  pred <- array(NA_real_, c(n, N_REPEATS, length(fsets), length(models)),
                dimnames = list(fish, NULL, fsets, models))
  fails <- list()
  set.seed(SEED + 1)
  message("Variant ", vname, ": ", n, " fish; feature sets: ", paste(fsets, collapse = ", "))
  for (r in seq_len(N_REPEATS)) {
    for (f in seq_len(N_FOLDS)) {
      te <- which(folds[, r] == f); tr <- which(folds[, r] != f)
      for (fs in fsets) {
        x <- get_x(fs, r, f)
        if (ncol(x) == 0) next
        xy <- prep_xy(x[tr, , drop = FALSE], x[te, , drop = FALSE])
        if (is.null(xy)) next
        for (mod in models) {
          p <- tryCatch(
            suppressWarnings(if (is_reg) fit_predict_reg(mod, xy$xtr, y[tr], xy$xte) else fit_predict(mod, xy$xtr, y[tr], xy$xte)),
            error = function(e) { fails[[length(fails) + 1]] <<- paste(fs, mod, conditionMessage(e)); rep(NA_real_, length(te)) })
          pred[te, r, fs, mod] <- as.vector(p)
        }
      }
    }
    message("  repeat ", r, "/", N_REPEATS)
  }
  if (length(fails)) {
    ft <- as.data.frame(table(failure = unlist(fails)))
    write_tsv(ft, paste0("model_failures_", vname, ".tsv"))
    message("  ", length(fails), " model fits failed (see model_failures_", vname, ".tsv)")
  }

  # Performance: one value per test fold (rows: repeats, columns: folds)
  fold_stat <- function(fs, mod, fun) {
    t(sapply(seq_len(N_REPEATS), function(r) sapply(seq_len(N_FOLDS), function(f) {
      te <- which(folds[, r] == f)
      fun(y[te], pred[te, r, fs, mod])
    })))
  }
  fold_auc <- list()
  perf <- do.call(rbind, lapply(fsets, function(fs) do.call(rbind, lapply(models, function(mod) {
    P <- pred[, , fs, mod]
    if (all(is.na(P))) return(NULL)
    if (is_reg) {
      # A constant prediction (a model that found nothing to use) carries no ranking: correlation 0
      rho  <- fold_stat(fs, mod, function(yy, pp) {
        if (anyNA(pp)) NA else if (sd(pp) == 0) 0 else suppressWarnings(cor(pp, yy, method = "spearman"))
      })
      rmse <- fold_stat(fs, mod, function(yy, pp) if (anyNA(pp)) NA else sqrt(mean((yy - pp)^2)))
      fold_auc[[paste(fs, mod)]] <<- rho
      nb <- nb_summary(as.vector(rho), N_FOLDS)
      return(data.frame(variant = vname, feature_set = fs, model = mod, n_fish = n,
                        spearman_mean = round(nb["mean"], 3), spearman_ci95_low = round(nb["low"], 3),
                        spearman_ci95_high = round(nb["high"], 3),
                        rmse_mean = round(mean(rmse, na.rm = TRUE), 3), sd_of_trait = round(sd(y), 3),
                        folds_ok = sum(!is.na(rho)), row.names = NULL))
    }
    A <- fold_stat(fs, mod, auc_fixed)
    fold_auc[[paste(fs, mod)]] <<- A
    nb  <- nb_summary(as.vector(A), N_FOLDS)
    cls <- P >= 0.5
    sens <- colMeans(cls[y == 1, , drop = FALSE]); spec <- colMeans(!cls[y == 0, , drop = FALSE])
    pm  <- rowMeans(P, na.rm = TRUE)                                  # each fish: mean of its out-of-fold predictions
    pc  <- pmin(pmax(pm, 1e-4), 1 - 1e-4)
    cal <- tryCatch(suppressWarnings(coef(glm(y ~ qlogis(pc), family = binomial()))), error = function(e) c(NA, NA))
    data.frame(variant = vname, feature_set = fs, model = mod, n_fish = n, n_low = sum(y == 0), n_high = sum(y == 1),
               auc_mean = round(nb["mean"], 3), auc_ci95_low = round(max(0, nb["low"]), 3), auc_ci95_high = round(min(1, nb["high"]), 3),
               auc_sd_between_repeats = round(sd(rowMeans(A, na.rm = TRUE)), 3),
               auc_lowest_repeat = round(min(rowMeans(A, na.rm = TRUE)), 3),
               auc_highest_repeat = round(max(rowMeans(A, na.rm = TRUE)), 3),
               sensitivity = round(mean(sens, na.rm = TRUE), 3), specificity = round(mean(spec, na.rm = TRUE), 3),
               balanced_accuracy = round(mean((sens + spec) / 2, na.rm = TRUE), 3),
               brier = round(mean(colMeans((P - y)^2), na.rm = TRUE), 3),
               calibration_intercept = round(cal[1], 2), calibration_slope = round(cal[2], 2),
               folds_ok = sum(!is.na(A)), row.names = NULL)
  }))))

  list(name = vname, note = v$note, y = y, fish = fish, pred = pred, perf = perf, fold_tab = fold_tab,
       is_reg = is_reg, fsets = fsets, models = models, hfe_sel = hfe_sel, hfe_full = hfe_full,
       hfe_log = hfe_log, clades = clades, used_hfe = use_hfe, folds = folds, fold_auc = fold_auc)
}

###############################################################################
# 7. RUN
###############################################################################
unknown <- setdiff(RUN_VARIANTS, names(VARIANTS))
if (length(unknown)) stop("Unknown variant(s): ", paste(unknown, collapse = ", "))
results <- list()
for (vn in RUN_VARIANTS) results[[vn]] <- run_variant(vn)
results <- results[!vapply(results, is.null, logical(1))]

# TaxaHFE commands still to run
sh_file <- file.path(ML_DIR, "run_hfe.sh")
writeLines(c("#!/bin/bash", "set -e", hfe_cmds), sh_file)
if (length(hfe_cmds) > 0 && RUN_HFE_FROM_R) {
  message("Running ", length(hfe_cmds), " TaxaHFE jobs (log: ml/run_hfe_log.txt)...")
  st <- system2("bash", shQuote(sh_file), stdout = file.path(ML_DIR, "run_hfe_log.txt"), stderr = file.path(ML_DIR, "run_hfe_log.txt"))
  if (st != 0) stop("TaxaHFE stopped with an error; see ", file.path(ML_DIR, "run_hfe_log.txt"))
  hfe_cmds <- character(0)
  for (vn in intersect(RUN_VARIANTS, HFE_VARIANTS)) results[[vn]] <- run_variant(vn)
}

###############################################################################
# 8. TABLES
###############################################################################
cls_res <- results[!vapply(results, function(r) r$is_reg, logical(1))]
reg_res <- results[vapply(results, function(r) r$is_reg, logical(1))]

perf_all <- do.call(rbind, lapply(cls_res, `[[`, "perf"))
if (!is.null(perf_all)) {
  write_tsv(perf_all, "performance_all_models.tsv")
  message("\nAUC of the ", PRIMARY_MODEL, " model (mean of fold-wise AUCs; corrected 95% CI):")
  print(perf_all[perf_all$model == PRIMARY_MODEL,
                 c("variant", "feature_set", "n_fish", "auc_mean", "auc_ci95_low", "auc_ci95_high", "balanced_accuracy")],
        row.names = FALSE)
}
if (length(reg_res)) {
  perf_reg <- do.call(rbind, lapply(reg_res, `[[`, "perf"))
  write_tsv(perf_reg, "performance_regression.tsv")
  message("\nContinuous lice load (fold-wise Spearman correlation between predicted and observed):")
  print(perf_reg, row.names = FALSE)
}
write_tsv(data.frame(variant = names(VARIANTS), description = vapply(VARIANTS, `[[`, character(1), "note"),
                     run = names(VARIANTS) %in% names(results),
                     with_hfe = names(VARIANTS) %in% names(results)[vapply(results, `[[`, logical(1), "used_hfe")]),
          "variants.tsv")
fold_all <- do.call(rbind, lapply(cls_res, `[[`, "fold_tab"))
if (!is.null(fold_all)) write_tsv(fold_all, "class_counts_per_fold.tsv")

# Out-of-fold predictions of the primary model (mean over repeats), per fish
pred_tab <- do.call(rbind, lapply(cls_res, function(r) {
  if (!PRIMARY_MODEL %in% r$models) return(NULL)
  do.call(rbind, lapply(r$fsets, function(fs) {
    data.frame(variant = r$name, feature_set = fs, sample = r$fish, observed_high = r$y,
               predicted_prob_high = round(rowMeans(r$pred[, , fs, PRIMARY_MODEL], na.rm = TRUE), 4))
  }))
}))
if (!is.null(pred_tab)) write_tsv(pred_tab, "predictions_primary_model.tsv")

# Paired comparisons of feature sets: difference of fold-wise AUCs on the same
# folds, tested with the corrected resampled t-test
pairs <- list(c("host", "host_hfe"), c("host", "host_asv"), c("hfe", "hfe_leaky"), c("asv", "hfe"),
              c("technical", "hfe"), c("host", "hfe"), c("host", "asv"))
cmp <- do.call(rbind, lapply(results, function(r) {
  mod <- if (PRIMARY_MODEL %in% r$models) PRIMARY_MODEL else r$models[1]
  do.call(rbind, lapply(pairs, function(pr) {
    k1 <- paste(pr[1], mod); k2 <- paste(pr[2], mod)
    if (!all(c(k1, k2) %in% names(r$fold_auc))) return(NULL)
    d  <- as.vector(r$fold_auc[[k2]] - r$fold_auc[[k1]])
    nb <- nb_summary(d, N_FOLDS)
    data.frame(variant = r$name, model = mod, metric = if (r$is_reg) "Spearman" else "AUC", set_1 = pr[1], set_2 = pr[2],
               value_1 = round(mean(r$fold_auc[[k1]], na.rm = TRUE), 3), value_2 = round(mean(r$fold_auc[[k2]], na.rm = TRUE), 3),
               difference = round(nb["mean"], 3), diff_ci95_low = round(nb["low"], 3), diff_ci95_high = round(nb["high"], 3),
               p_corrected_t = signif(2 * pt(-abs(nb["mean"] / nb["se"]), nb["df"]), 3), row.names = NULL)
  }))
}))
if (!is.null(cmp)) {
  write_tsv(cmp, "feature_set_comparisons.tsv")
  message("\nFeature-set comparisons (set_2 minus set_1):")
  print(cmp, row.names = FALSE)
}

# HFE: how often each feature is selected, and the features chosen on all fish
stab_all <- NULL
for (r in results) {
  if (!r$used_hfe) next
  sel   <- unlist(r$hfe_sel)
  if (length(sel) == 0) next
  n_run <- N_REPEATS * N_FOLDS
  freq  <- sort(table(sel), decreasing = TRUE)
  info  <- r$clades$info[match(names(freq), r$clades$info$clade), ]
  cm    <- r$clades$mat[, names(freq), drop = FALSE]
  # Association with the outcome on all fish. Descriptive only: the features
  # were chosen using the outcome, so these p-values are not independent tests.
  if (r$is_reg) {
    eff <- apply(cm, 2, function(z) suppressWarnings(cor(z, r$y, method = "spearman")))
    pv  <- apply(cm, 2, function(z) suppressWarnings(cor.test(z, r$y, method = "spearman")$p.value))
  } else {
    eff <- 100 * (colMeans(cm[r$y == 1, , drop = FALSE]) - colMeans(cm[r$y == 0, , drop = FALSE]))
    pv  <- apply(cm, 2, function(z) suppressWarnings(wilcox.test(z ~ r$y)$p.value))
  }
  stab  <- data.frame(variant = r$name, rank = info$rank, taxon = info$name,
                      selected_in_pct_of_training_sets = round(100 * as.vector(freq) / n_run, 1),
                      selected_on_all_fish = names(freq) %in% r$hfe_full,
                      mean_abund_pct = round(100 * colMeans(cm), 4),
                      prevalence_pct = round(100 * colMeans(cm > 0), 1),
                      association = if (r$is_reg) "Spearman rho with the trait" else "High minus Low, mean abundance (% points)",
                      effect = signif(eff, 3), p = signif(pv, 3),
                      clade = names(freq), stringsAsFactors = FALSE, row.names = NULL)
  # Benjamini-Hochberg adjustment over the features chosen on all fish
  stab$p_BH_all_fish_features <- NA_real_
  af <- stab$selected_on_all_fish
  stab$p_BH_all_fish_features[af] <- signif(p.adjust(stab$p[af], method = "BH"), 3)
  stab_all <- rbind(stab_all, stab)
  nsel <- vapply(r$hfe_sel, length, integer(1))
  message("\nVariant ", r$name, ": TaxaHFE selected ", round(mean(nsel), 1), " features per training set (range ",
          min(nsel), "-", max(nsel), "); ", length(r$hfe_full), " on all fish.")
}
if (!is.null(stab_all)) {
  write_tsv(stab_all, "hfe_feature_stability.tsv")
  for (vn in unique(stab_all$variant)) {
    print(head(stab_all[stab_all$variant == vn, c("variant", "rank", "taxon", "selected_in_pct_of_training_sets", "selected_on_all_fish")], 15),
          row.names = FALSE)
  }
}
hfe_log_all <- do.call(rbind, lapply(results, `[[`, "hfe_log"))
if (!is.null(hfe_log_all)) {
  write_tsv(hfe_log_all, "hfe_column_matching.tsv")
  if (any(!hfe_log_all$matched)) {
    warning(sum(!hfe_log_all$matched), " TaxaHFE output columns could not be matched to a clade and were left out; ",
            "see ml/hfe_column_matching.tsv", immediate. = TRUE)
  }
}

###############################################################################
# 9. FIGURES
###############################################################################
pdf(file.path(ML_DIR, "ml_results.pdf"), width = 10, height = 7)

# (1) Lice load and burden classes
xs <- seq(min(tc), max(tc), length.out = 400)
dens <- data.frame(x = xs, y = mix$w * dnorm(xs, mix$m1, mix$sd1) + (1 - mix$w) * dnorm(xs, mix$m2, mix$sd2))
print(ggplot(data.frame(lice_load = tc), aes(x = lice_load)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "grey80", color = "grey50") +
  geom_line(data = dens, aes(x = x, y = y), color = "red") +
  geom_vline(xintercept = c(LOW_CUT, HIGH_CUT), linetype = "dashed") +
  geom_vline(xintercept = MEDIAN_CUT, linetype = "dotted") +
  theme_minimal() +
  labs(title = "Lice load and burden classes",
       subtitle = paste0("Dashed: Low <= ", round(LOW_CUT, 3), " and High >= ", round(HIGH_CUT, 3),
                         " lice per cm; dotted: median (", round(MEDIAN_CUT, 3), ")"),
       x = "Lice per cm", y = "Density"))

if (!is.null(perf_all)) {
  # (2) Primary analysis: every model and feature set
  pp <- perf_all[perf_all$variant == RUN_VARIANTS[1], ]
  pp$feature_set <- factor(pp$feature_set, levels = unique(pp$feature_set))
  print(ggplot(pp, aes(x = model, y = auc_mean, color = feature_set)) +
    geom_hline(yintercept = 0.5, linetype = "dashed") +
    geom_pointrange(aes(ymin = auc_ci95_low, ymax = auc_ci95_high), position = position_dodge(width = 0.7)) +
    theme_minimal() +
    labs(title = paste0("Cross-validated AUC, analysis '", RUN_VARIANTS[1], "'"),
         subtitle = paste0("Mean fold-wise AUC over ", N_REPEATS, " repeats of ", N_FOLDS, "-fold cross-validation; bar: corrected 95% CI"),
         x = "Model", y = "AUC", color = "Feature set"))

  # (3) All analyses, primary model
  pv <- perf_all[perf_all$model == PRIMARY_MODEL, ]
  pv$variant <- factor(pv$variant, levels = names(cls_res))
  print(ggplot(pv, aes(x = variant, y = auc_mean, color = feature_set)) +
    geom_hline(yintercept = 0.5, linetype = "dashed") +
    geom_pointrange(aes(ymin = auc_ci95_low, ymax = auc_ci95_high), position = position_dodge(width = 0.7)) +
    theme_minimal() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
    labs(title = paste0("Sensitivity analyses, model: ", PRIMARY_MODEL),
         subtitle = "Mean fold-wise AUC with corrected 95% CI",
         x = "Analysis", y = "AUC", color = "Feature set"))

  # (4) ROC curves and (5) calibration, primary analysis and model
  r <- cls_res[[1]]
  if (PRIMARY_MODEL %in% r$models) {
    # ROC curve of every test fold, averaged at fixed false-positive rates
    grid <- seq(0, 1, by = 0.02)
    pa <- perf_all[perf_all$variant == r$name & perf_all$model == PRIMARY_MODEL, ]
    roc_df <- do.call(rbind, lapply(r$fsets, function(fs) {
      tprs <- NULL
      for (rr in seq_len(N_REPEATS)) for (f in seq_len(N_FOLDS)) {
        te <- which(r$folds[, rr] == f)
        pp <- r$pred[te, rr, fs, PRIMARY_MODEL]
        if (anyNA(pp) || length(unique(r$y[te])) < 2) next
        ro <- pROC::roc(r$y[te], pp, levels = c(0, 1), direction = "<", quiet = TRUE)
        fpr <- rev(1 - ro$specificities); tpr <- rev(ro$sensitivities)
        tprs <- rbind(tprs, vapply(grid, function(g) max(tpr[fpr <= g]), numeric(1)))
      }
      if (is.null(tprs)) return(NULL)
      i <- which(pa$feature_set == fs)
      data.frame(fs = sprintf("%s: %.2f (%.2f-%.2f)", fs, pa$auc_mean[i], pa$auc_ci95_low[i], pa$auc_ci95_high[i]),
                 fpr = grid, tpr = colMeans(tprs))
    }))
    print(ggplot(roc_df, aes(x = fpr, y = tpr, color = fs)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
      geom_line() + coord_equal() + theme_minimal() +
      labs(title = paste0("ROC curves, analysis '", r$name, "', model: ", PRIMARY_MODEL),
           subtitle = "Average of the test-fold ROC curves; legend: mean AUC (corrected 95% CI)",
           x = "1 - specificity", y = "Sensitivity", color = "Feature set"))

    cal_df <- do.call(rbind, lapply(r$fsets, function(fs) {
      pm <- rowMeans(r$pred[, , fs, PRIMARY_MODEL], na.rm = TRUE)
      if (anyNA(pm)) return(NULL)
      b <- cut(pm, breaks = unique(quantile(pm, 0:5 / 5)), include.lowest = TRUE)
      data.frame(fs = fs, predicted = as.vector(tapply(pm, b, mean)), observed = as.vector(tapply(r$y, b, mean)))
    }))
    print(ggplot(cal_df, aes(x = predicted, y = observed)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
      geom_line() + geom_point() + facet_wrap(~fs) + coord_equal(xlim = c(0, 1), ylim = c(0, 1)) + theme_minimal() +
      labs(title = paste0("Calibration, analysis '", r$name, "', model: ", PRIMARY_MODEL),
           subtitle = "Fish grouped into fifths of predicted probability",
           x = "Mean predicted probability of High burden", y = "Observed proportion of High burden"))
  }
}

# (6) HFE feature stability
if (!is.null(stab_all)) {
  st <- head(stab_all[stab_all$variant == stab_all$variant[1], ], 30)
  st$label <- factor(paste0(st$taxon, " (", st$rank, ")"), levels = rev(unique(paste0(st$taxon, " (", st$rank, ")"))))
  print(ggplot(st, aes(x = label, y = selected_in_pct_of_training_sets, fill = selected_on_all_fish)) +
    geom_col() + coord_flip() + theme_minimal() +
    labs(title = "How often TaxaHFE selects each feature",
         subtitle = paste0("Across ", N_REPEATS * N_FOLDS, " training sets (", N_REPEATS, " repeats x ", N_FOLDS, " folds)"),
         x = NULL, y = "Selected in % of training sets", fill = "Also selected\non all fish"))
}
invisible(dev.off())

###############################################################################
# 10. RECORD
###############################################################################
saveRDS(lapply(results, function(r) r[c("name", "note", "y", "fish", "folds", "pred", "perf", "fold_auc", "hfe_sel", "hfe_full")]),
        file.path(ML_DIR, "ml_results.rds"))
writeLines(c(
  paste("date:", format(Sys.time())),
  paste("cross-validation:", N_FOLDS, "folds x", N_REPEATS, "repeats; seed", SEED),
  "folds: stratified by class; fish of the same family kept in the same fold (except variant 'ungrouped')",
  paste("models:", paste(MODELS, collapse = ", "), "; primary:", PRIMARY_MODEL),
  "fixed settings: randomForest 500 trees, default mtry; xgboost 10 rounds, default eta and depth;",
  "  SVM radial kernel, cost 1, default gamma; k-NN k = 5; LDA tol = 0",
  "tuned inside each training fold (inner 5-fold CV): glmnet lambda (lasso, ridge, elastic net); elastic-net alpha in 0.25, 0.5, 0.75",
  "features centred and scaled with training-fold means and SDs",
  "AUC: computed within each test fold (levels Low = 0, High = 1; higher prediction = High) and averaged;",
  "  95% CI and feature-set comparisons: corrected resampled t-test (Nadeau & Bengio 2003)",
  paste("TaxaHFE image:", HFE_DOCKER_IMAGE, "; lowest_level", HFE_LOWEST_LEVEL, "; seed", HFE_SEED),
  paste("variants run:", paste(names(results), collapse = ", ")),
  paste("variants with TaxaHFE:", paste(names(results)[vapply(results, `[[`, logical(1), "used_hfe")], collapse = ", ")),
  capture.output(sessionInfo())
), file.path(ML_DIR, "ml_run_record.txt"))

if (length(hfe_cmds) > 0) {
  message("\nTaxaHFE still has to run for ", length(hfe_cmds), " training sets. In a terminal with Docker running:\n  bash ",
          shQuote(sh_file), "\nthen run this script again to add the HFE feature sets.")
}
message("ML: tables and figures in ", ML_DIR)
