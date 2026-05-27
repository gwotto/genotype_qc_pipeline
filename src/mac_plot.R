library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)

args          <- commandArgs(trailingOnly = TRUE)
mac_threshold <- as.integer(args[1])

# ── Read frequency files ───────────────────────────────────────────────────────
before <- read.table("freq_before.frq", header = TRUE) |>
    mutate(
        MAF     = as.numeric(MAF),
        NCHROBS = as.integer(NCHROBS),
        MAC     = as.integer(round(MAF * NCHROBS))
    )

after <- read.table("freq_after.frq", header = TRUE) |>
    mutate(
        MAF     = as.numeric(MAF),
        NCHROBS = as.integer(NCHROBS),
        MAC     = as.integer(round(MAF * NCHROBS))
    )

# ── Barplot: SNPs removed vs retained ─────────────────────────────────────────
n_total       <- nrow(before)
n_monomorphic <- sum(before$MAC == 0, na.rm = TRUE)
n_low_mac     <- sum(before$MAC > 0 & before$MAC < mac_threshold, na.rm = TRUE)
n_retained    <- nrow(after)

counts <- tibble(
    category = factor(
        c("Monomorphic\n(MAC = 0)",
          paste0("Low MAC\n(0 < MAC < ", mac_threshold, ")"),
          "Retained"),
        levels = c(
            "Monomorphic\n(MAC = 0)",
            paste0("Low MAC\n(0 < MAC < ", mac_threshold, ")"),
            "Retained"
        )
    ),
    n    = c(n_monomorphic, n_low_mac, n_retained),
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
        title    = "SNPs removed by MAC filter",
        subtitle = paste0(
            "Total input: ", comma(n_total), " SNPs  |  ",
            "Retained: ",    comma(n_retained), " SNPs  |  ",
            "Removed: ",     comma(n_total - n_retained)
        ),
        x = NULL,
        y = "Number of SNPs"
    ) +
    theme_bw()

# ── MAC distribution before and after ─────────────────────────────────────────
mac_combined <- bind_rows(
    before |> filter(MAC > 0) |> mutate(stage = "Before filtering"),
    after  |> filter(MAC > 0) |> mutate(stage = "After filtering")
) |>
    mutate(stage = factor(stage, levels = c("Before filtering", "After filtering")))

# Sensible binwidth: 1 for small thresholds, wider for larger datasets
bw <- max(1L, as.integer(mac_threshold / 5L))

p_dist <- ggplot(mac_combined, aes(x = MAC, fill = stage)) +
    geom_histogram(binwidth = bw, boundary = 0, colour = "white", linewidth = 0.2) +
    geom_vline(
        xintercept = mac_threshold,
        linetype   = "dashed",
        colour     = "grey30"
    ) +
    annotate(
        "text",
        x      = mac_threshold,
        y      = Inf,
        label  = paste0("threshold\n(MAC = ", mac_threshold, ")"),
        hjust  = -0.1,
        vjust  = 1.3,
        size   = 3,
        colour = "grey30"
    ) +
    scale_fill_manual(values = c(
        "Before filtering" = "#92c5de",
        "After filtering"  = "#2166ac"
    )) +
    scale_x_continuous(labels = comma) +
    scale_y_continuous(labels = comma) +
    facet_wrap(~stage, ncol = 1) +
    labs(
        title = "MAC distribution before and after filtering",
        x     = "Minor allele count",
        y     = "Number of SNPs",
        fill  = NULL
    ) +
    theme_bw() +
    theme(legend.position = "none")

# ── Combine and save ───────────────────────────────────────────────────────────
p_bar / p_dist +
    plot_annotation(
        title   = "MAC quality control",
        caption = paste0(
            "Monomorphic SNPs removed first (MAC = 0), ",
            "then SNPs with MAC < ", mac_threshold, " removed"
        )
    )

ggsave("mac_plot.png", width = 8, height = 10, dpi = 150)
