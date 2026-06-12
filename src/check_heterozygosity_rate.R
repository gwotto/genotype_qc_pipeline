library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)

# ── Load and compute heterozygosity rates ──────────────────────────────────────
het <- read.table("R_check.het", header = TRUE) |>
    mutate(HET_RATE = (N.NM. - O.HOM.) / N.NM.)

het_mean <- mean(het$HET_RATE)
het_sd   <- sd(het$HET_RATE)
lo       <- het_mean - 3 * het_sd
hi       <- het_mean + 3 * het_sd

het <- het |>
    mutate(outlier = HET_RATE < lo | HET_RATE > hi)

n_outlier <- sum(het$outlier)
n_total   <- nrow(het)

# ── Panel 1: Histogram ────────────────────────────────────────────────────────
p_hist <- ggplot(het, aes(x = HET_RATE)) +
    geom_histogram(
        aes(fill = outlier),
        bins      = 60,
        colour    = "white",
        linewidth = 0.2
    ) +
    geom_vline(xintercept = het_mean, colour = "grey30",  linewidth = 0.8) +
    geom_vline(xintercept = lo,       colour = "#d6604d", linewidth = 0.8, linetype = "dashed") +
    geom_vline(xintercept = hi,       colour = "#d6604d", linewidth = 0.8, linetype = "dashed") +
    annotate("text", x = het_mean, y = Inf,
             label = "mean", hjust = -0.2, vjust = 1.5, size = 3, colour = "grey30") +
    annotate("text", x = lo, y = Inf,
             label = "mean − 3 SD", hjust = 1.05, vjust = 1.5, size = 3, colour = "#d6604d") +
    annotate("text", x = hi, y = Inf,
             label = "mean + 3 SD", hjust = -0.05, vjust = 1.5, size = 3, colour = "#d6604d") +
    scale_fill_manual(
        values = c("FALSE" = "#2166ac", "TRUE" = "#d6604d"),
        labels = c("FALSE" = "Retained", "TRUE" = "Outlier (±3 SD)"),
        name   = NULL
    ) +
    scale_x_continuous(labels = percent_format(accuracy = 0.1)) +
    scale_y_continuous(labels = comma) +
    labs(
        title    = "Heterozygosity rate distribution",
        subtitle = paste0(
            "Samples: ", comma(n_total), "  |  ",
            "Outliers (> ±3 SD): ", n_outlier
        ),
        x = "Observed heterozygosity rate",
        y = "Number of samples"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")

# ── Panel 2: Per-sample sorted dot plot ───────────────────────────────────────
het_sorted <- het |>
    arrange(HET_RATE) |>
    mutate(rank = row_number())

p_dot <- ggplot(het_sorted, aes(x = rank, y = HET_RATE, colour = outlier)) +
    geom_point(size = 0.6, alpha = 0.6) +
    geom_hline(yintercept = het_mean, colour = "grey30",  linewidth = 0.6) +
    geom_hline(yintercept = lo,       colour = "#d6604d", linewidth = 0.8, linetype = "dashed") +
    geom_hline(yintercept = hi,       colour = "#d6604d", linewidth = 0.8, linetype = "dashed") +
    scale_colour_manual(
        values = c("FALSE" = "#2166ac", "TRUE" = "#d6604d"),
        guide  = "none"
    ) +
    scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
    scale_x_continuous(labels = comma) +
    labs(
        title = "Per-sample heterozygosity (sorted)",
        x     = "Sample rank",
        y     = "Observed heterozygosity rate"
    ) +
    theme_bw()

# ── Combine and save ───────────────────────────────────────────────────────────
p_hist / p_dot +
    plot_annotation(
        title   = "Heterozygosity outlier detection",
        caption = paste0(
            "Outlier threshold: mean ± 3 SD  |  ",
            "Mean: ",       round(het_mean * 100, 2), "%  |  ",
            "SD: ",         round(het_sd   * 100, 3), "%  |  ",
            "Thresholds: [", round(lo * 100, 2), "%, ", round(hi * 100, 2), "%]"
        )
    )

ggsave("heterozygosity.png", width = 8, height = 10, dpi = 150)
