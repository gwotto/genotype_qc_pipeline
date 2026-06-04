## =============================================================================
## ancestry_pca_plot.R
## =============================================================================
##
## PURPOSE
##   Projects study samples onto a 1000 Genomes Phase 3 PCA, assigns
##   superpopulation ancestry labels using a random forest classifier trained
##   on 1000G samples with known labels, and produces a scatter plot.
##
## USAGE
##   Rscript ancestry_pca_plot.R <eigenvec> <psam> <n_pcs> <output_prefix>
##
## ARGUMENTS
##   eigenvec        PLINK2 .eigenvec file (study + 1000G samples, joint PCA)
##   psam            1000G .psam file with SuperPop column (tab-separated)
##   n_pcs           Number of PCs to use for classification (e.g. 10)
##   output_prefix   Prefix for output files
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
if (length(args) != 4) {
  stop("Usage: Rscript ancestry_pca_plot.R <eigenvec> <psam> <n_pcs> <output_prefix>")
}

eigenvec_file  <- args[1]
psam_file      <- args[2]
n_pcs          <- as.integer(args[3])
output_prefix  <- args[4]

# -----------------------------------------------------------------------------
# Load PCA results
# PLINK2 .eigenvec format: #FID IID PC1 PC2 ... PCn
# -----------------------------------------------------------------------------
pcs <- read_table(eigenvec_file, comment = "") |>
  rename(FID = `#FID`)

pc_cols <- paste0("PC", seq_len(n_pcs))

# -----------------------------------------------------------------------------
# Load 1000G sample metadata
# .psam format: #IID PAT MAT SEX SuperPop Population (tab-separated)
# -----------------------------------------------------------------------------
psam <- read_tsv(psam_file, comment = "") |>
  rename(IID = `#IID`) |>
  select(IID, SuperPop)

# -----------------------------------------------------------------------------
# Split into 1000G reference (labelled) and study samples (unlabelled)
# -----------------------------------------------------------------------------
ref <- pcs |>
  inner_join(psam, join_by(IID)) |>
  filter(!is.na(SuperPop))

study <- pcs |>
  anti_join(psam, join_by(IID))

if (nrow(ref) == 0) stop("No 1000G reference samples found in eigenvec file — check IID matching.")
if (nrow(study) == 0) stop("No study samples found in eigenvec file — check IID matching.")

message(sprintf("Reference samples (1000G): %d", nrow(ref)))
message(sprintf("Study samples:             %d", nrow(study)))

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
    predicted_prob = apply(pred_probs, 1, max)
  )

# -----------------------------------------------------------------------------
# Write ancestry assignments
# -----------------------------------------------------------------------------
assignments <- study |>
  select(FID, IID, predicted_pop, predicted_prob) |>
  rename(superpop = predicted_pop, probability = predicted_prob)

write_tsv(assignments, paste0(output_prefix, "_ancestry_assignments.tsv"))

message(sprintf("Ancestry assignments written to %s_ancestry_assignments.tsv", output_prefix))
message("Predicted superpopulation counts:")
assignments |> count(superpop) |> print()

# -----------------------------------------------------------------------------
# Plot: study samples overlaid on 1000G reference clusters
# -----------------------------------------------------------------------------
superpop_colours <- c(
  AFR = "#E41A1C",
  AMR = "#FF7F00",
  EAS = "#4DAF4A",
  EUR = "#377EB8",
  SAS = "#984EA3"
)

plot_data <- bind_rows(
  ref   |> mutate(dataset = "1000G",  pop = SuperPop),
  study |> mutate(dataset = "Study",  pop = as.character(predicted_pop))
)

p <- ggplot(plot_data, aes(x = PC1, y = PC2, colour = pop, shape = dataset, alpha = dataset, size = dataset)) +
  geom_point() +
  scale_colour_manual(values = superpop_colours, name = "Superpopulation") +
  scale_shape_manual(values = c("1000G" = 16, "Study" = 17), name = "Dataset") +
  scale_alpha_manual(values = c("1000G" = 0.3, "Study" = 0.9), name = "Dataset") +
  scale_size_manual(values  = c("1000G" = 1.0, "Study" = 2.0), name = "Dataset") +
  labs(
    title = "Ancestry PCA: study samples vs 1000 Genomes Phase 3",
    x = "PC1",
    y = "PC2"
  ) +
  theme_bw() +
  theme(legend.position = "right")

ggsave(paste0(output_prefix, "_ancestry_pca.png"), plot = p, width = 8, height = 6, dpi = 150)
message(sprintf("PCA plot saved to %s_ancestry_pca.png", output_prefix))
