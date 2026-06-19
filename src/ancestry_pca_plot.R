## =============================================================================
## ancestry_pca_plot.R
## =============================================================================
##
## PURPOSE
##   Assigns superpopulation ancestry labels to all study samples using a
##   random forest classifier trained on 1000G samples with known labels.
##   PCA is computed on 1000G only; all study samples (related and unrelated)
##   are projected onto the reference PC space via PLINK2 variant weights,
##   giving stable ancestry axes independent of study composition.
##
## USAGE
##   Rscript ancestry_pca_plot.R <eigenvec> <sscore> <psam> <n_pcs> <output_prefix>
##
## ARGUMENTS
##   eigenvec        PLINK2 .eigenvec file (1000G reference PCA)
##   sscore          PLINK2 .sscore file   (ALL study samples projected onto reference PCs)
##   psam            1000G .psam file with SuperPop column (tab-separated)
##   n_pcs           Number of PCs to use for classification (e.g. 10)
##   output_prefix   Prefix for output files
##   eigenval        PLINK2 .eigenval file (PC eigenvalues, one per line)
##
## OUTPUTS
##   <output_prefix>_ancestry_pca.png         PC1 vs PC2 scatter plot
##   <output_prefix>_ancestry_assignments.tsv Per-sample ancestry assignments
## =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(randomForest)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop("Usage: Rscript ancestry_pca_plot.R <eigenvec> <sscore> <psam> <n_pcs> <output_prefix> <eigenval>")
}

eigenvec_file  <- args[1]
sscore_file    <- args[2]
psam_file      <- args[3]
n_pcs          <- as.integer(args[4])
output_prefix  <- args[5]
eigenval_file  <- args[6]

pc_cols <- paste0("PC", seq_len(n_pcs))

# -----------------------------------------------------------------------------
# Load 1000G sample metadata
# .psam format: #IID PAT MAT SEX SuperPop Population (tab-separated)
# -----------------------------------------------------------------------------
psam <- read_tsv(psam_file, comment = "") |>
  rename(IID = `#IID`) |>
  select(IID, SuperPop)

# -----------------------------------------------------------------------------
# Extract 1000G reference PC scores from the eigenvec.
# The eigenvec contains 1000G samples only (PCA was reference-only).
# -----------------------------------------------------------------------------
pcs_joint <- read_table(eigenvec_file, comment = "") |>
  rename_with(~ sub("^#", "", .)) |>
  mutate(across(c(FID, IID), as.character))

# The 1000G sample name may be in IID or FID depending on how the reference FAM
# was built (some FAMs use FID=sample_name, IID=0 or vice versa). Check which
# column has more overlap with psam IIDs and normalise to IID.
n_iid <- sum(pcs_joint$IID %in% psam$IID)
n_fid <- sum(pcs_joint$FID %in% psam$IID)
if (n_fid > n_iid) {
  message(sprintf("Using FID column as sample ID (%d matches vs %d for IID).", n_fid, n_iid))
  pcs_joint <- mutate(pcs_joint, IID = FID)
}

ref <- pcs_joint |>
  inner_join(psam, join_by(IID)) |>
  filter(!is.na(SuperPop))

if (nrow(ref) == 0) stop("No 1000G reference samples found in eigenvec — check that eigenvec IIDs match psam IIDs.")

# -----------------------------------------------------------------------------
# Load projected PC scores for ALL study samples (including related).
# PLINK2 .sscore format (when --score-col-nums is used):
#   FID  IID  ALLELE_CT  NAMED_ALLELE_DOSAGE_SUM  SCORE1_AVG  SCORE2_AVG ...
#
# SCORE_AVG and eigenvec values are on different scales:
#   eigenvec_k = SCORE_AVG_k * 2 / sqrt(eigenvalue_k)
# This arises because SCORE_SUM = sqrt(eigenvalue_k) * M_pca * eigenvec_k
# and SCORE_AVG = SCORE_SUM / ALLELE_CT = SCORE_SUM / (2 * M_scored),
# which reduces to SCORE_AVG = sqrt(eigenvalue_k) * eigenvec_k / 2
# (the M terms cancel for any missingness level).
# -----------------------------------------------------------------------------
eigenvalues   <- scan(eigenval_file, quiet = TRUE)[seq_len(n_pcs)]
scale_factors <- setNames(2 / sqrt(eigenvalues), paste0("SCORE", seq_len(n_pcs), "_AVG"))

sscore_raw <- read_table(sscore_file, comment = "") |>
  rename_with(~ sub("^#", "", .)) |>
  mutate(FID = as.character(FID), IID = as.character(IID))

score_avg_cols <- paste0("SCORE", seq_len(n_pcs), "_AVG")
missing_cols   <- setdiff(score_avg_cols, names(sscore_raw))
if (length(missing_cols) > 0) {
  stop(paste("Expected score columns not found in sscore file:",
             paste(missing_cols, collapse = ", "),
             "\nFound:", paste(names(sscore_raw), collapse = ", ")))
}

study <- sscore_raw |>
  mutate(across(all_of(score_avg_cols), ~ . * scale_factors[cur_column()])) |>
  rename_with(~ pc_cols, all_of(score_avg_cols)) |>
  select(FID, IID, all_of(pc_cols)) |>
  anti_join(psam, join_by(IID))

if (nrow(study) == 0) stop("No study samples found in sscore file — check IID matching.")

message(sprintf("Reference samples (1000G): %d", nrow(ref)))
message(sprintf("Study samples (total):     %d", nrow(study)))

# -----------------------------------------------------------------------------
# Train random forest on 1000G PCs
# -----------------------------------------------------------------------------
set.seed(42)
rf_model <- randomForest(
  x = ref |> select(all_of(pc_cols)),
  y = as.factor(ref$SuperPop),
  ntree = 500,
  importance = FALSE
)

# -----------------------------------------------------------------------------
# Predict ancestry for study samples
# -----------------------------------------------------------------------------
pred_probs <- predict(rf_model, newdata = study |> select(all_of(pc_cols)), type = "prob")

study <- study |>
  mutate(
    predicted_pop  = predict(rf_model, newdata = study |> select(all_of(pc_cols))),
    predicted_prob = apply(pred_probs, 1, max),
    # Assign ancestry only if probability > 0.8; otherwise mark as unassigned
    superpop_final = ifelse(predicted_prob > 0.8, as.character(predicted_pop), "unassigned")
  )

# -----------------------------------------------------------------------------
# Write ancestry assignments
# Ancestry is assigned only if Random Forest prediction probability > 0.8
# (high confidence). Samples with probability ≤ 0.8 are marked "unassigned"
# for manual review or alternative assignment methods.
# -----------------------------------------------------------------------------
assignments <- study |>
  select(FID, IID, superpop_final, predicted_prob) |>
  mutate(FID = ifelse(is.na(FID), IID, FID)) |>
  rename(superpop = superpop_final, probability = predicted_prob)

write_tsv(assignments, paste0(output_prefix, "_ancestry_assignments.tsv"))

message(sprintf("Ancestry assignments written to %s_ancestry_assignments.tsv", output_prefix))
message("Superpopulation assignments (probability > 0.8):")
assignments |> count(superpop) |> print()
n_unassigned <- sum(assignments$superpop == "unassigned")
if (n_unassigned > 0) {
  message(sprintf("Note: %d samples marked 'unassigned' (probability ≤ 0.8)", n_unassigned))
}

# -----------------------------------------------------------------------------
# Plot: study samples coloured by predicted ancestry
# -----------------------------------------------------------------------------
superpop_colours <- c(
  AFR = "#E41A1C",      # Red
  AMR = "#FF7F00",      # Orange
  EAS = "#4DAF4A",      # Green
  EUR = "#377EB8",      # Blue
  SAS = "#984EA3",      # Purple
  unassigned = "#CCCCCC" # Light grey
)

plot_data <- study |> mutate(pop = as.character(superpop_final))

p <- ggplot(plot_data, aes(x = PC1, y = PC2, colour = pop)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_colour_manual(values = superpop_colours, name = "Superpopulation") +
  labs(
    title = "Ancestry PCA: study samples projected onto 1000G reference",
    x = "PC1",
    y = "PC2"
  ) +
  theme_bw() +
  theme(legend.position = "right")

ggsave(paste0(output_prefix, "_ancestry_pca.png"), plot = p, width = 8, height = 6, dpi = 150)
message(sprintf("PCA plot saved to %s_ancestry_pca.png", output_prefix))
