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

# Load ancestry assignments
ancestry <- read_tsv(ancestry_file) |>
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
