library(ggplot2)
library(dplyr)
library(patchwork)

# Thresholds passed from WDL as command-line arguments
args           <- commandArgs(trailingOnly = TRUE)
geno_threshold <- as.numeric(args[1])
mind_threshold <- as.numeric(args[2])

geno <- read.table("geno_sweep.txt", header = TRUE) |>
    mutate(threshold = as.numeric(threshold))

mind <- read.table("mind_sweep.txt", header = TRUE) |>
    mutate(threshold = as.numeric(threshold))

p_geno <- ggplot(geno, aes(x = threshold, y = n_snps)) +
    geom_line(colour = "#2166ac", linewidth = 0.8) +
    geom_point(colour = "#2166ac", size = 3) +
    geom_vline(
        xintercept = geno_threshold,
        linetype   = "dashed",
        colour     = "grey40"
    ) +
    annotate(
        "text",
        x      = geno_threshold,
        y      = max(geno$n_snps),
        label  = paste0("chosen\n(", geno_threshold * 100, "%)"),
        hjust  = -0.1,
        size   = 3,
        colour = "grey40"
    ) +
    scale_x_continuous(labels = scales::percent_format()) +
    labs(
        title    = "SNP retention across --geno thresholds",
        subtitle = "SNPs retained as missingness threshold is relaxed",
        x        = "Maximum missingness per SNP (--geno)",
        y        = "Number of SNPs retained"
    ) +
    theme_bw()

p_mind <- ggplot(mind, aes(x = threshold, y = n_samples)) +
    geom_line(colour = "#d6604d", linewidth = 0.8) +
    geom_point(colour = "#d6604d", size = 3) +
    geom_vline(
        xintercept = mind_threshold,
        linetype   = "dashed",
        colour     = "grey40"
    ) +
    annotate(
        "text",
        x      = mind_threshold,
        y      = max(mind$n_samples),
        label  = paste0("chosen\n(", mind_threshold * 100, "%)"),
        hjust  = -0.1,
        size   = 3,
        colour = "grey40"
    ) +
    scale_x_continuous(labels = scales::percent_format()) +
    labs(
        title    = "Sample retention across --mind thresholds",
        subtitle = "Samples retained as missingness threshold is relaxed",
        x        = "Maximum missingness per sample (--mind)",
        y        = "Number of samples retained"
    ) +
    theme_bw()

p_geno / p_mind +
    plot_annotation(
        title   = "Call rate threshold sweep",
        caption = "Dashed line indicates the threshold applied in the QC pipeline"
    )

ggsave("threshold_sweep.png", width = 8, height = 8, dpi = 150)
