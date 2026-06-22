## =============================================================================
## ancestry_pca_plot.R
## =============================================================================
##
## PURPOSE
##   Assigns superpopulation ancestry labels to all study samples using a
##   random forest classifier trained on 1000G reference samples with known
##   labels. PCA is computed on the 1000G reference panel only; all study
##   samples (including related pairs) are projected into the reference PC
##   space using PLINK2 variant weights, giving stable axes independent of
##   study composition.
##
## USAGE
##   Rscript ancestry_pca_plot.R \
##       <eigenvec> <sscore> <psam> <n_pcs> <output_prefix> \
##       <eigenval> <ancestry_prob_threshold>
##
## ARGUMENTS
##   eigenvec                 PLINK2 .eigenvec  — reference-panel PC scores
##   sscore                   PLINK2 .sscore    — study samples projected via --score
##   psam                     1000G  .psam      — sample metadata with SuperPop column
##   n_pcs                    Number of PCs to use for classification (e.g. 10)
##   output_prefix            Prefix for all output files
##   eigenval                 PLINK2 .eigenval  — one eigenvalue per line
##   ancestry_prob_threshold  RF probability threshold for assignment (0–1, e.g. 0.5)
##                            Samples below this are labelled "unassigned".
##
## SCALE CONVERSION
##   PLINK2 --score outputs SCORE_AVG, which differs from eigenvec values by a
##   per-PC factor: eigenvec_k = SCORE_AVG_k × 2 / sqrt(eigenvalue_k).
##   This script applies that conversion so projected study scores are on the
##   same scale as the reference eigenvec before classification.
##
## OUTPUTS
##   <output_prefix>_ancestry_pca.png         Grid of PC pair scatter plots
##                                            (PC1/2, PC3/4, PC5/6, PC7/8, PC9/10)
##   <output_prefix>_ancestry_assignments.tsv Per-sample assignments:
##                                            FID, IID, superpop, probability
## =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(randomForest)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
  stop("Usage: Rscript ancestry_pca_plot.R <eigenvec> <sscore> <psam> <n_pcs> <output_prefix> <eigenval> <ancestry_prob_threshold>")
}

eigenvec_file  <- args[1]
sscore_file    <- args[2]
psam_file      <- args[3]
n_pcs          <- as.integer(args[4])
output_prefix  <- args[5]
eigenval_file  <- args[6]
ancestry_prob_threshold <- as.numeric(args[7])

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
    # Assign ancestry only if probability > predicted_prob; otherwise mark as unassigned
    superpop_final = ifelse(predicted_prob > ancestry_prob_threshold, as.character(predicted_pop), "unassigned")
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
message(sprintf("Superpopulation assignments (probability > %.2f):", ancestry_prob_threshold))
assignments |> count(superpop) |> print()
n_unassigned <- sum(assignments$superpop == "unassigned")
if (n_unassigned > 0) {
  message(sprintf("Note: %d samples marked 'unassigned' (probability ≤ %.2f)", n_unassigned, ancestry_prob_threshold))
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

# Create multiple PC plots (PC1 vs PC2, PC3 vs PC4, etc.)
p1 <- ggplot(plot_data, aes(x = PC1, y = PC2, colour = pop)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_colour_manual(values = superpop_colours, name = "Superpopulation") +
  labs(x = "PC1", y = "PC2") +
  theme_bw() +
  theme(legend.position = "right")

p2 <- ggplot(plot_data, aes(x = PC3, y = PC4, colour = pop)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_colour_manual(values = superpop_colours, name = "Superpopulation") +
  labs(x = "PC3", y = "PC4") +
  theme_bw() +
  theme(legend.position = "none")

p3 <- ggplot(plot_data, aes(x = PC5, y = PC6, colour = pop)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_colour_manual(values = superpop_colours, name = "Superpopulation") +
  labs(x = "PC5", y = "PC6") +
  theme_bw() +
  theme(legend.position = "none")

p4 <- ggplot(plot_data, aes(x = PC7, y = PC8, colour = pop)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_colour_manual(values = superpop_colours, name = "Superpopulation") +
  labs(x = "PC7", y = "PC8") +
  theme_bw() +
  theme(legend.position = "none")

p5 <- ggplot(plot_data, aes(x = PC9, y = PC10, colour = pop)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_colour_manual(values = superpop_colours, name = "Superpopulation") +
  labs(x = "PC9", y = "PC10") +
  theme_bw() +
  theme(legend.position = "none")

# Combine all plots
p <- (p1 + p2 + p3) / (p4 + p5) +
  plot_annotation(
    title = "Ancestry PCA: study samples projected onto 1000G reference",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 14))
  )

ggsave(paste0(output_prefix, "_ancestry_pca.png"), plot = p, width = 8, height = 6, dpi = 150)
message(sprintf("PCA plot saved to %s_ancestry_pca.png", output_prefix))
