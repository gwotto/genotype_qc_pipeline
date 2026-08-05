#!/usr/bin/env Rscript
# Create comprehensive sample QC status table
# Single table with ancestry and relatedness information
# Users can derive their own filtering logic as needed

library(tidyverse)

# Arguments
args <- commandArgs(trailingOnly = TRUE)
ancestry_file <- args[1]      # ancestry_assignments.tsv
relatedness_file <- args[2]   # relatedness_flagged_samples.txt
output_prefix <- args[3]      # prefix for output
pca_file <- if (length(args) >= 4) args[4] else "NA"  # covariate PCs (FID IID PC1..PCn); "NA" = none

# Load ancestry assignments
ancestry <- read_tsv(ancestry_file, col_types = cols(FID = col_character(), IID = col_character())) |>
  select(FID, IID, superpop, probability) |>
  rename(ancestry = superpop, ancestry_prob = probability)

# Load relatedness flagged samples
relatedness_flagged <- tryCatch(
  {
    read_tsv(relatedness_file, col_names = c("FID", "IID")) |>
      mutate(related = TRUE) |>
      distinct()
  },
  error = function(e) {
    # If file is empty or doesn't exist, return empty tibble
    tibble(FID = character(), IID = character(), related = logical())
  }
)

# Create QC status table
qc_status <- ancestry |>
  left_join(relatedness_flagged, by = c("FID", "IID")) |>
  mutate(
    related = ifelse(is.na(related), FALSE, related)
  ) |>
  select(FID, IID, ancestry, ancestry_prob, related)

# Optionally append within-cohort covariate PCs (FID IID PC1..PCn).
# Coerce join keys to character on both sides to avoid type-mismatch errors.
if (!is.na(pca_file) && pca_file != "NA" && file.exists(pca_file)) {
  pcs <- read_tsv(pca_file, col_types = cols(FID = col_character(), IID = col_character()))
  qc_status <- qc_status |>
    left_join(pcs, by = c("FID", "IID"))
  cat(sprintf("Appended %d covariate PCs from: %s\n", ncol(pcs) - 2, pca_file))
}

# Print summary statistics
cat("\n=== QC Summary ===\n")
cat(sprintf("Total samples: %d\n", nrow(qc_status)))
cat(sprintf("Ancestry assigned: %d (%.1f%%)\n", 
  sum(qc_status$ancestry != "unassigned"), 
  100 * mean(qc_status$ancestry != "unassigned")))
cat(sprintf("Unrelated samples: %d (%.1f%%)\n", 
  sum(!qc_status$related), 
  100 * mean(!qc_status$related)))
cat(sprintf("\nQC status table written to: %s_sample_qc_status.tsv\n\n", output_prefix))

# Write the QC table
output_file <- paste0(output_prefix, "_sample_qc_status.tsv")
write_tsv(qc_status, output_file)
