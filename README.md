# Skin mucus microbiota and sea lice burden in Atlantic salmon

Analysis code for the study of the skin mucus microbiota of Atlantic salmon
(*Salmo salar*), sampled before challenge, as a predictor of subsequent sea
lice (*Caligus rogercresseyi*) burden.

Collaborative PhD project: University of Chile, Laval University and Benchmark
Genetics. Contact: Lucas Venegas (lucas.venegas@ug.uchile.cl).

## Repository layout

```
scripts/   the analysis pipeline, numbered in the order it is run
legacy/    the scripts deposited with the first submission (superseded)
```

## Pipeline

Scripts 01 and 02 must be run first, in that order. Scripts 03 to 07 each read
the files written by 02 and can be run in any order.

| Script | Purpose | Main inputs | Writes |
|---|---|---|---|
| `01_dada2_asv_pipeline.R` | Primer removal (cutadapt), DADA2 denoising per sequencing run, SILVA taxonomy, tree, phyloseq object | Paired FASTQ files, metadata table | `phyloseq_caligus_microbiome_SILVA.rds`, ASV tables, `logs/` |
| `02_qc_decontam.R` | Contaminant removal (decontam), taxonomy and prevalence filters (skin and water separately), rarefaction, normalization alternatives, control and batch checks | Output of 01 | `Skin_rare_SILVA.rds`, `Skin_ps_SILVA.rds`, `Water_ps_SILVA.rds`, `bacteria_unfiltered_physeq.rds`, other `.rds`, `logs/qc/` |
| `03_alpha_beta_diversity.R` | Alpha diversity (mixed models) and beta diversity (PERMANOVA, dispersion, ordination) against burden | Output of 02 | `diversity/` |
| `04_relative_abundance.R` | Relative abundance figures and tables, water and skin | Output of 02 | `abundance/` |
| `05_phenotype_correlations.R` | Fish description, data checks, correlations between body size and lice measures | Output of 02 | `phenotype/` |
| `06_ml_burden_prediction.R` | Prediction of burden with eight algorithms; TaxaHFE inside cross-validation folds; baselines and sensitivity analyses | Output of 02; Docker | `ml/` |
| `07_source_tracking.R` | Overlap between skin, water and kit controls; SourceTracker2 inputs and summaries | Outputs of 01 and 02; SourceTracker2 | `sourcetracker/` |

## How to run

Each script uses the current working directory as the project folder: it reads
its inputs from there and writes its outputs there. Start R in the folder that
holds the FASTQ files and the metadata table, then run for example:

```r
source("path/to/Caligus/scripts/01_dada2_asv_pipeline.R")
```

or, from a terminal in that folder:

```bash
Rscript path/to/Caligus/scripts/01_dada2_asv_pipeline.R
```

All settings are in the "User Parameters" block at the top of each script.
Before the first run, edit these machine-specific paths:

- `01_dada2_asv_pipeline.R`: `CUTADAPT_BIN`, `SILVA_TRAIN`, `SILVA_SPEC`,
  `METADATA_TSV`
- `07_source_tracking.R`: `SOURCETRACKER_BIN`

Both executables must be the program file inside the conda environment (for
example `.../envs/cutadapt/bin/cutadapt`), not the environment folder.

### Two-step scripts

`06_ml_burden_prediction.R` and `07_source_tracking.R` call an external tool.
Each is run in three steps:

1. Run the R script. It writes the tool's input files and a shell script
   (`ml/run_hfe.sh` or `sourcetracker/run_sourcetracker.sh`).
2. Run that shell script in a terminal (`bash ml/run_hfe.sh`). It can be
   interrupted and restarted; finished jobs are kept.
3. Run the R script again to read the results.

For script 06, the TaxaHFE results are tied to the cross-validation folds. If
the fish, the burden definition, `N_FOLDS` or `N_REPEATS` change, delete the
`ml/hfe` folder and repeat the three steps.

## Input data

- **FASTQ files:** paired-end 16S rRNA V4 reads (515F/806R), named
  `<sample>_V..._R1_001.fastq.gz` and `..._R2_001.fastq.gz`.
- **Metadata table:** tab-delimited, one row per sample, sample IDs in the
  first column. Required columns: `Pittag`, `Sex`, `Total_caligus`,
  `Final_Weight`, `Final_Length`, `Fish_family`. Water samples have
  `Pittag = CTRL`; extraction-kit negative controls have `Sex = Kit`.
- **Sequencing runs:** the samples of the second run are listed in
  `BATCH2_SAMPLES` at the top of scripts 01 and 02.

Sequence data: [add accession number].

## Software

R (version 4.2 or later) with these packages:

| Script | Packages |
|---|---|
| 01 | dada2, ShortRead, Biostrings, DECIPHER, phangorn, phyloseq, dplyr, ggplot2 |
| 02 | phyloseq, decontam, vegan, ape, phangorn, ggplot2; optional: metagenomeSeq |
| 03 | phyloseq, vegan, ape, VGAM, lme4, car, ggplot2, patchwork; optional: DHARMa |
| 04 | phyloseq, ggplot2, patchwork |
| 05 | phyloseq, ggplot2 |
| 06 | phyloseq, VGAM, randomForest, e1071, glmnet, MASS, class, pROC, ggplot2; optional: xgboost |
| 07 | phyloseq, biomformat, ggplot2 |

External tools:

- cutadapt (script 01)
- SILVA v138.2 reference files formatted for DADA2 (script 01)
- TaxaHFE, run through Docker (script 06). The command line is that of
  TaxaHFE 2.4. Replace `latest` in `HFE_DOCKER_IMAGE` by the version used.
- SourceTracker2 (script 07)

Each script records the R session, package versions and its settings in its
output folder. To record the full environment, run from the project folder:

```bash
conda env export -n cutadapt > cutadapt_env.yml
conda env export -n st2 > st2_env.yml
```

```r
renv::init(); renv::snapshot()
```

and add the resulting files to the repository.

## Notes on the analysis

- **Burden classes** are defined from lice per cm of body length with a
  two-component Gaussian mixture; alternative definitions (company categories
  on raw counts, median split, continuous lice load) are analysed in scripts
  03 and 06.
- **Cross-validation** in script 06 is five-fold, repeated, stratified by
  burden class, with fish of the same family kept in the same fold. Feature
  selection (TaxaHFE), scaling and hyperparameter tuning are done within each
  training fold. The feature set `hfe_leaky` reproduces the procedure of the
  first submission (selection on all fish) and is reported for comparison
  only.
- **Random seeds** are fixed in every script and listed in its parameter block.

## Citation

[Add the reference once published.]
