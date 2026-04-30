library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)

args          <- commandArgs(trailingOnly = TRUE)
maf_threshold <- as.numeric(args[1])

# ── Read frequency files ───────────────────────────────────────────────────────
before <- read.table("freq_before.frq", header = TRUE) |>
    mutate(MAF = as.numeric(MAF))

after <- read.table("freq_after.frq", header = TRUE) |>
    mutate(MAF = as.numeric(MAF))

# ── Barplot: SNPs removed vs retained ─────────────────────────────────────────
n_total       <- nrow(before)
n_monomorphic <- sum(before$MAF == 0, na.rm = TRUE)
n_low_maf     <- sum(before$MAF > 0 & before$MAF < maf_threshold, na.rm = TRUE)
n_retained    <- nrow(after)

counts <- tibble(
    category = factor(
        c("Monomorphic\n(MAF = 0)",
          paste0("Low MAF\n(0 < MAF < ", maf_threshold, ")"),
          "Retained"),
        levels = c(
            "Monomorphic\n(MAF = 0)",
            paste0("Low MAF\n(0 < MAF < ", maf_threshold, ")"),
            "Retained"
        )
    ),
    n    = c(n_monomorphic, n_low_maf, n_retained),
    fill = c("removed", "removed", "retained")
)

p_bar <- ggplot(counts, aes(x = category, y = n, fill = fill)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = comma(n)), vjust = -0.4, size = 3.5) +
    scale_fill_manual(
        values = c("removed" = "#d6604d", "retained" = "#2166ac"),
        guide  = "none"
    ) +
    scale_y_continuous(
        labels = comma,
        expand = expansion(mult = c(0, 0.12))
    ) +
    labs(
        title    = "SNPs removed by MAF filter",
        subtitle = paste0(
            "Total input: ", comma(n_total), " SNPs  |  ",
            "Retained: ",    comma(n_retained), " SNPs  |  ",
            "Removed: ",     comma(n_total - n_retained)
        ),
        x = NULL,
        y = "Number of SNPs"
    ) +
    theme_bw()

# ── MAF distribution before and after ─────────────────────────────────────────
maf_combined <- bind_rows(
    before |> filter(MAF > 0) |> mutate(stage = "Before filtering"),
    after  |> filter(MAF > 0) |> mutate(stage = "After filtering")
) |>
    mutate(stage = factor(stage, levels = c("Before filtering", "After filtering")))

p_dist <- ggplot(maf_combined, aes(x = MAF, fill = stage)) +
    geom_histogram(binwidth = 0.01, boundary = 0, colour = "white", linewidth = 0.2) +
    geom_vline(
        xintercept = maf_threshold,
        linetype   = "dashed",
        colour     = "grey30"
    ) +
    annotate(
        "text",
        x      = maf_threshold,
        y      = Inf,
        label  = paste0("threshold\n(", maf_threshold * 100, "%)"),
        hjust  = -0.1,
        vjust  = 1.3,
        size   = 3,
        colour = "grey30"
    ) +
    scale_fill_manual(values = c(
        "Before filtering" = "#92c5de",
        "After filtering"  = "#2166ac"
    )) +
    scale_x_continuous(labels = percent_format()) +
    scale_y_continuous(labels = comma) +
    facet_wrap(~stage, ncol = 1) +
    labs(
        title = "MAF distribution before and after filtering",
        x     = "Minor allele frequency",
        y     = "Number of SNPs",
        fill  = NULL
    ) +
    theme_bw() +
    theme(legend.position = "none")

# ── Combine and save ───────────────────────────────────────────────────────────
p_bar / p_dist +
    plot_annotation(
        title   = "MAF quality control",
        caption = paste0(
            "Monomorphic SNPs removed first (MAF = 0), ",
            "then SNPs with MAF < ", maf_threshold * 100, "% removed"
        )
    )

ggsave("maf_plot.png", width = 8, height = 10, dpi = 150)