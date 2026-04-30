library(ggplot2)
library(dplyr)
library(patchwork)

args       <- commandArgs(trailingOnly = TRUE)
eigenvec_f <- args[1]
eigenval_f <- args[2]
n_pcs      <- as.integer(args[3])

# ── Read eigenvectors ──────────────────────────────────────────────────────────
# PLINK2 eigenvec has a header: #FID IID PC1 PC2 ...
eigenvec <- read.table(eigenvec_f, header = TRUE, comment.char = "") |>
    rename(FID = 1, IID = 2)

pc_cols <- paste0("PC", seq_len(n_pcs))
colnames(eigenvec)[3:(n_pcs + 2)] <- pc_cols

# ── Read eigenvalues and compute variance explained ───────────────────────────
eigenval <- read.table(eigenval_f, header = FALSE) |>
    mutate(
        PC             = paste0("PC", row_number()),
        variance       = V1 / sum(V1) * 100
    )

# ── Scree plot ─────────────────────────────────────────────────────────────────
p_scree <- ggplot(eigenval, aes(x = reorder(PC, -variance), y = variance)) +
    geom_col(fill = "#2166ac", width = 0.7) +
    geom_text(
        aes(label = sprintf("%.1f%%", variance)),
        vjust = -0.4, size = 3
    ) +
    scale_y_continuous(
        expand = expansion(mult = c(0, 0.12)),
        labels = function(x) paste0(x, "%")
    ) +
    labs(
        title    = "Scree plot",
        subtitle = "Variance explained by each principal component",
        x        = "Principal component",
        y        = "Variance explained (%)"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ── PC scatter plots ───────────────────────────────────────────────────────────
pct <- function(pc) {
    sprintf("%.1f%%", eigenval$variance[eigenval$PC == pc])
}

make_pc_plot <- function(x_pc, y_pc) {
    ggplot(eigenvec, aes(x = .data[[x_pc]], y = .data[[y_pc]])) +
        geom_point(alpha = 0.4, size = 1.2, colour = "#2166ac") +
        labs(
            x = paste0(x_pc, " (", pct(x_pc), ")"),
            y = paste0(y_pc, " (", pct(y_pc), ")")
        ) +
        theme_bw()
}

p_pc12 <- make_pc_plot("PC1", "PC2") + labs(title = "PC1 vs PC2")
p_pc13 <- make_pc_plot("PC1", "PC3") + labs(title = "PC1 vs PC3")
p_pc23 <- make_pc_plot("PC2", "PC3") + labs(title = "PC2 vs PC3")
p_pc14 <- make_pc_plot("PC1", "PC4") + labs(title = "PC1 vs PC4")

# ── Combine ────────────────────────────────────────────────────────────────────
p_scree / (p_pc12 + p_pc13) / (p_pc23 + p_pc14) +
    plot_annotation(
        title    = "PCA — pre-imputation diagnostic",
        subtitle = "Inspect for population outliers. No samples are removed at this step.",
        caption  = paste0(
            "LD-pruned SNPs used. ",
            "Points far from the main cluster may indicate ancestry outliers."
        )
    )

ggsave("pca_plot.png", width = 10, height = 12, dpi = 150)