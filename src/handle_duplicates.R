library(dplyr)

# File path passed from WDL as command-line argument
args     <- commandArgs(trailingOnly = TRUE)
fam_file <- args[1]

fam <- read.table(
    fam_file,
    header    = FALSE,
    col.names = c("FID", "IID", "PAT", "MAT", "SEX", "PHENO")
)

# Identify duplicate IIDs and append a numeric suffix to make them unique
duplicates <- fam |>
    group_by(IID) |>
    tally() |>
    filter(n > 1) |>
    pull(IID)

if (length(duplicates) > 0) {
    for (dup in duplicates) {
        indices           <- which(fam$IID == dup)
        fam[indices, "IID"] <- paste0(dup, "_", seq_along(indices))
    }
}

write.table(fam, fam_file, quote = FALSE, row.names = FALSE, col.names = FALSE)
