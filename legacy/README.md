# Original scripts

These are the six scripts as deposited with the first submission of the
manuscript. They are kept unchanged for reference and are **superseded by the
scripts in `../scripts/`**. Do not use them for new analyses.

| Original script | Replaced by |
|---|---|
| `Preprocessing_16S.R` | `01_dada2_asv_pipeline.R` |
| `QC&DECONTAM.R` | `02_qc_decontam.R` |
| `Alpha_diversity.R`, `Beta_diversity.R` | `03_alpha_beta_diversity.R` |
| `Script_CMS.R`, abundance plots | `04_relative_abundance.R` |
| `Correlations.R` | `05_phenotype_correlations.R` |
| `Script_CMS.R`, machine learning and TaxaHFE | `06_ml_burden_prediction.R` |
| `Script_CMS.R`, SourceTracker | `07_source_tracking.R` |

Main differences in the revised scripts:

- TaxaHFE feature selection, scaling and tuning are done inside each
  cross-validation training fold. In `Script_CMS.R` features were selected on
  all fish before cross-validation.
- AUC is computed with a fixed direction, and SVM probabilities are read by
  class name.
- Chloroplast and mitochondrial ASVs are removed at the Order and Family ranks
  where SILVA places them.
- Controls are identified by metadata labels, not by row position.
- The two sequencing runs are denoised with separate DADA2 error models.
- Skin and water samples are prevalence-filtered separately.
- Every script runs on its own from the files written by the previous ones.
