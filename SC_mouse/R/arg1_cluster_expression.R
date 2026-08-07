library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)

options(stringsAsFactors = FALSE, future.globals.maxSize = 4 * 1024^3)

# ============================================================
# Arg1 expression across clusters: WT vs KO, by EAE condition
#
# For each cluster (stable_20_clusters) x EAE level x genotype:
#   - Mean SCT expression (log1p-normalised)
#   - Percent cells expressing (> 0)
#   - Cell count
#   - Wilcoxon rank-sum test (WT vs KO), BH-adjusted
#
# Output: CSV table + dot plot + violin plot + bar plot + README
# ============================================================

local({
  ca <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", grep("^--file=", ca, value = TRUE))
  paths_file <- if (length(f) == 1) {
    file.path(dirname(normalizePath(f)), "00_paths.R")
  } else {
    file.path("00_paths.R")
  }
  if (!file.exists(paths_file)) {
    for (cand in c("R/00_paths.R", "00_paths.R")) {
      if (file.exists(cand)) {
        paths_file <- cand
        break
      }
    }
  }
  source(paths_file)
})

OUTPUT_DIR <- ensure_out_dir("arg1_expression")
PROCESSED_DIR <- OUTPUT_DIR

GENE <- "Arg1"
CLUSTER_COL <- "stable_20_clusters"
EAE_LEVELS <- c("NO EAE", "EAE")

source_methods("load_annotated_object.R")

# ============================================================
# Load data
# ============================================================

cat("Loading annotated Seurat object...\n")
so <- load_annotated_object()

available_genes <- rownames(LayerData(so, assay = "SCT", layer = "data"))
if (!(GENE %in% available_genes)) {
    stop(sprintf("Gene '%s' not found in SCT assay. Available genes: %d", GENE, length(available_genes)))
}
cat(sprintf("Verified: %s is present in SCT assay (%d total features)\n", GENE, length(available_genes)))

# ============================================================
# Extract Arg1 expression + metadata
# ============================================================

arg1_expr <- as.numeric(LayerData(so, assay = "SCT", layer = "data")[GENE, ])
meta <- so[[]]
meta$arg1_expr <- arg1_expr
meta$cluster <- as.character(meta[[CLUSTER_COL]])
meta$genotype <- ifelse(meta$condition == "ctrl", "WT", "KO")

rm(so)
gc()

# ============================================================
# Compute per-group summary statistics
# ============================================================

cat("Computing per-group statistics...\n")

summary_df <- meta %>%
    filter(eae %in% EAE_LEVELS, !is.na(cluster)) %>%
    group_by(eae, cluster, genotype) %>%
    summarise(
        mean_expr = mean(arg1_expr),
        pct_expressing = mean(arg1_expr > 0) * 100,
        n_cells = n(),
        .groups = "drop"
    )

# All clusters present in the data (ensures complete grid for BH correction)
all_clusters <- sort(unique(as.character(meta$cluster[!is.na(meta$cluster)])))

summary_wide <- summary_df %>%
    pivot_wider(
        id_cols = c(eae, cluster),
        names_from = genotype,
        values_from = c(mean_expr, pct_expressing, n_cells),
        names_glue = "{.value}_{genotype}"
    ) %>%
    rename(
        mean_WT = mean_expr_WT,
        mean_KO = mean_expr_KO,
        pct_WT = pct_expressing_WT,
        pct_KO = pct_expressing_KO,
        n_WT = n_cells_WT,
        n_KO = n_cells_KO
    ) %>%
    # Ensure all cluster x EAE combos exist (even if 0 cells in both genotypes)
    complete(eae = EAE_LEVELS, cluster = all_clusters)

# Replace NA counts with 0 (clusters absent in one genotype)
summary_wide <- summary_wide %>%
    mutate(
        across(starts_with("n_"), ~ replace_na(.x, 0L)),
        across(starts_with("mean_"), ~ replace_na(.x, 0)),
        across(starts_with("pct_"), ~ replace_na(.x, 0))
    )

# ============================================================
# Wilcoxon tests: WT vs KO per cluster x EAE
# ============================================================

cat("Running Wilcoxon rank-sum tests...\n")

cell_data <- meta %>%
    filter(eae %in% EAE_LEVELS, !is.na(cluster))

wilcox_results <- summary_wide %>%
    rowwise() %>%
    mutate(pvalue = {
        cells <- cell_data %>%
            filter(eae == .env$eae, cluster == .env$cluster)
        wt_vals <- cells$arg1_expr[cells$genotype == "WT"]
        ko_vals <- cells$arg1_expr[cells$genotype == "KO"]

        if (length(wt_vals) == 0 || length(ko_vals) == 0) {
            NA_real_
        } else if (all(wt_vals == 0) && all(ko_vals == 0)) {
            NA_real_
        } else {
            tryCatch(
                wilcox.test(wt_vals, ko_vals, exact = FALSE)$p.value,
                error = function(e) NA_real_
            )
        }
    }) %>%
    ungroup()

# BH correction within each EAE level
wilcox_results <- wilcox_results %>%
    group_by(eae) %>%
    mutate(padj = p.adjust(pvalue, method = "BH")) %>%
    ungroup()

# Sort: by EAE, then by max mean expression descending
wilcox_results <- wilcox_results %>%
    mutate(max_mean = pmax(mean_WT, mean_KO, na.rm = TRUE)) %>%
    arrange(eae, desc(max_mean)) %>%
    select(-max_mean)

# Ensure cluster is numeric-sorted within eae for the final table
wilcox_results$cluster <- as.integer(wilcox_results$cluster)
wilcox_results <- wilcox_results %>%
    arrange(eae, desc(pmax(mean_WT, mean_KO)))

# ============================================================
# Save CSV
# ============================================================

csv_path <- file.path(OUTPUT_DIR, "arg1_expression_by_cluster.csv")
write.csv(wilcox_results, csv_path, row.names = FALSE)
file.copy(csv_path, file.path(PROCESSED_DIR, "arg1_expression_by_cluster.csv"), overwrite = TRUE)
cat(sprintf("Saved: %s\n", csv_path))

# Print top clusters
cat("\n=== Top Arg1-expressing clusters ===\n")
for (eae_lev in EAE_LEVELS) {
    cat(sprintf("\n%s:\n", eae_lev))
    sub <- wilcox_results %>% filter(eae == eae_lev) %>% head(5)
    for (i in seq_len(nrow(sub))) {
        row <- sub[i, ]
        cat(sprintf("  Cluster %d: WT=%.3f (%0.1f%%), KO=%.3f (%0.1f%%), n_WT=%d, n_KO=%d, padj=%s\n",
            row$cluster, row$mean_WT, row$pct_WT, row$mean_KO, row$pct_KO,
            row$n_WT, row$n_KO,
            ifelse(is.na(row$padj), "NA", sprintf("%.4g", row$padj))
        ))
    }
}

# ============================================================
# Prepare plotting data
# ============================================================

# Numeric cluster ordering: 1, 2, 3, ... 20
cluster_order <- as.character(sort(unique(as.integer(cell_data$cluster[!is.na(cell_data$cluster)]))))

plot_data <- cell_data %>%
    filter(!is.na(cluster)) %>%
    mutate(
        cluster = factor(as.integer(cluster), levels = as.integer(cluster_order)),
        genotype = factor(genotype, levels = c("WT", "KO")),
        eae = factor(eae, levels = EAE_LEVELS)
    )

# ============================================================
# Visualisation 1: Dot plot
# ============================================================

cat("\nGenerating dot plot...\n")

dot_data <- summary_df %>%
    mutate(
        cluster = factor(as.integer(cluster), levels = as.integer(cluster_order)),
        genotype = factor(genotype, levels = c("WT", "KO")),
        eae = factor(eae, levels = EAE_LEVELS)
    ) %>%
    filter(!is.na(cluster))

p_dot <- ggplot(dot_data, aes(
    x = cluster,
    y = genotype,
    size = pct_expressing,
    colour = mean_expr
)) +
    geom_point() +
    facet_wrap(~ eae, ncol = 1, scales = "free_x") +
    scale_colour_viridis_c(option = "magma", direction = -1, name = "Mean SCT\nexpression") +
    scale_size_continuous(range = c(1, 10), name = "% expressing") +
    labs(
        title = paste0(GENE, " expression across clusters (WT vs KO)"),
        x = "Cluster",
        y = ""
    ) +
    theme_classic(base_size = 14) +
    theme(
        strip.text = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 16, face = "bold"),
        legend.position = "right"
    )

dot_pdf <- file.path(OUTPUT_DIR, "arg1_dotplot.pdf")
dot_png <- file.path(OUTPUT_DIR, "arg1_dotplot.png")
cairo_pdf(dot_pdf, width = 14, height = 7)
print(p_dot)
dev.off()
ggsave(dot_png, p_dot, width = 14, height = 7, dpi = 300)
file.copy(dot_pdf, file.path(PROCESSED_DIR, "arg1_dotplot.pdf"), overwrite = TRUE)
file.copy(dot_png, file.path(PROCESSED_DIR, "arg1_dotplot.png"), overwrite = TRUE)
cat(sprintf("Saved: %s\n", dot_pdf))

# ============================================================
# Visualisation 2: Violin + boxplot per cluster (WT vs KO)
# ============================================================

cat("Generating violin plot...\n")

# Significance stars helper
sig_stars <- function(p) {
    ifelse(is.na(p), "",
    ifelse(p < 0.001, "***",
    ifelse(p < 0.01, "**",
    ifelse(p < 0.05, "*", "ns"))))
}

# Annotation: significance above each cluster
annot <- wilcox_results %>%
    mutate(
        cluster = factor(as.integer(cluster), levels = as.integer(cluster_order)),
        eae = factor(eae, levels = EAE_LEVELS),
        label = sig_stars(padj)
    ) %>%
    filter(!is.na(cluster), label != "")

# y-position for annotations (above the highest violin in that cluster)
annot_y <- plot_data %>%
    group_by(eae, cluster) %>%
    summarise(ymax = max(arg1_expr, na.rm = TRUE), .groups = "drop")

annot <- annot %>%
    left_join(annot_y, by = c("eae", "cluster")) %>%
    mutate(y_label = ymax + 0.05)

p_violin <- ggplot(plot_data, aes(
    x = cluster,
    y = arg1_expr,
    fill = genotype
)) +
    geom_violin(
        position = position_dodge(width = 0.9),
        alpha = 0.7,
        trim = TRUE,
        scale = "width",
        linewidth = 0.3,
        na.rm = TRUE
    ) +
    geom_boxplot(
        position = position_dodge(width = 0.9),
        width = 0.15,
        outlier.size = 0.3,
        outlier.alpha = 0.4,
        alpha = 0.85,
        na.rm = TRUE
    ) +
    geom_text(
        data = annot,
        aes(x = cluster, y = y_label, label = label),
        inherit.aes = FALSE,
        size = 3, vjust = 0
    ) +
    facet_wrap(~ eae, nrow = 2, scales = "free_y") +
    scale_fill_manual(
        values = c("WT" = "#4393C3", "KO" = "#D6604D"),
        name = "Genotype"
    ) +
    labs(
        title = paste0(GENE, " expression: WT vs KO per cluster"),
        x = "Cluster",
        y = paste0(GENE, " SCT expression (log1p)")
    ) +
    theme_classic(base_size = 13) +
    theme(
        strip.text = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 11),
        plot.title = element_text(size = 16, face = "bold"),
        legend.position = "top"
    )

violin_pdf <- file.path(OUTPUT_DIR, "arg1_violins.pdf")
violin_png <- file.path(OUTPUT_DIR, "arg1_violins.png")
cairo_pdf(violin_pdf, width = 18, height = 10)
print(p_violin)
dev.off()
ggsave(violin_png, p_violin, width = 18, height = 10, dpi = 300)
file.copy(violin_pdf, file.path(PROCESSED_DIR, "arg1_violins.pdf"), overwrite = TRUE)
file.copy(violin_png, file.path(PROCESSED_DIR, "arg1_violins.png"), overwrite = TRUE)
cat(sprintf("Saved: %s\n", violin_pdf))

# ============================================================
# Visualisation 3: Bar plot (mean expression + SE by sample)
# ============================================================

cat("Generating bar plot...\n")

bar_summary <- plot_data %>%
    filter(!is.na(arg1_expr)) %>%
    group_by(eae, cluster, genotype) %>%
    summarise(
        mean_expr = mean(arg1_expr, na.rm = TRUE),
        se_expr = sd(arg1_expr, na.rm = TRUE) / sqrt(n()),
        .groups = "drop"
    )

p_bar <- ggplot(bar_summary, aes(
    x = cluster,
    y = mean_expr,
    fill = genotype
)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.8) +
    geom_errorbar(
        aes(ymin = mean_expr - se_expr, ymax = mean_expr + se_expr),
        position = position_dodge(width = 0.8), width = 0.25
    ) +
    facet_wrap(~ eae, nrow = 2) +
    scale_fill_manual(
        values = c("WT" = "#4393C3", "KO" = "#D6604D"),
        name = "Genotype"
    ) +
    labs(
        title = paste0(GENE, " mean expression per cluster (WT vs KO)"),
        x = "Cluster",
        y = paste0("Mean ", GENE, " SCT expression (log1p)")
    ) +
    theme_classic(base_size = 14) +
    theme(
        strip.text = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 11),
        plot.title = element_text(size = 16, face = "bold"),
        legend.position = "top"
    )

bar_pdf <- file.path(OUTPUT_DIR, "arg1_barplot.pdf")
bar_png <- file.path(OUTPUT_DIR, "arg1_barplot.png")
cairo_pdf(bar_pdf, width = 14, height = 10)
print(p_bar)
dev.off()
ggsave(bar_png, p_bar, width = 14, height = 10, dpi = 300)
file.copy(bar_pdf, file.path(PROCESSED_DIR, "arg1_barplot.pdf"), overwrite = TRUE)
file.copy(bar_png, file.path(PROCESSED_DIR, "arg1_barplot.png"), overwrite = TRUE)
cat(sprintf("Saved: %s\n", bar_pdf))

# ============================================================
# README.md
# ============================================================

readme_path <- file.path(OUTPUT_DIR, "README.md")
writeLines(c(
    "# Arg1 Expression Analysis Across Clusters",
    "",
    sprintf("**Date**: %s", Sys.Date()),
    "",
    "## Background",
    "",
    "This analysis quantifies arginase-1 (Arg1) expression across all 20 stable",
    "clusters (`stable_20_clusters`) in the Krzak et al. 2025 scRNA-seq dataset,",
    "stratified by genotype (WT vs Sucnr1 KO) and EAE condition (NO EAE, EAE).",
    "",
    "## Data",
    "",
    "- **Seurat object**: `seurat_object.rds`",
    "  (SCT-normalised, Harmony-corrected, protein-coding genes)",
    "- **Clustering**: `stable_20_clusters` derived from ClustAssess stability",
    "  (`clustassess_object.rds`)",
    "- **Expression layer**: SCT `data` (log1p-normalised corrected counts)",
    "",
    "## Methods",
    "",
    "### Expression quantification",
    "",
    "For each combination of EAE condition (NO EAE, EAE), cluster (1-20), and",
    "genotype (WT, KO), the following metrics are computed:",
    "",
    "- **Mean expression**: mean of SCT log1p-normalised values across all cells",
    "  in the group (same approach as used in the bubble plot and heatmap scripts)",
    "- **Percent expressing**: percentage of cells with Arg1 expression > 0",
    "- **Cell count**: number of cells in the group",
    "",
    "### Statistical test",
    "",
    "- **Wilcoxon rank-sum test** comparing Arg1 SCT expression between WT and KO",
    "  cells within each cluster, separately for each EAE condition",
    "- This is the same test used by Seurat `FindMarkers()` in the genotype DEG",
    "  analysis (`deg_ctrl_vs_ko.R`)",
    "- P-values are BH-adjusted within each EAE condition (20 tests per condition)",
    "- Clusters with 0 cells in either genotype: test not performed (p = NA)",
    "",
    "## Output Files",
    "",
    "| File | Description |",
    "|------|-------------|",
    "| `arg1_expression_by_cluster.csv` | Summary table with expression metrics and test results |",
    "| `arg1_dotplot.pdf/png` | Dot plot: all 20 clusters, dot size = % expressing, colour = mean expression |",
    "| `arg1_violins.pdf/png` | Violin + boxplot: WT (blue) vs KO (red) side by side per cluster, significance stars |",
    "| `arg1_barplot.pdf/png` | Bar plot (mean ± SEM) for all clusters |",
    "",
    "## How to read the plots",
    "",
    "### Dot plot",
    "- Each dot represents one genotype within one cluster",
    "- Dot size encodes the percentage of cells expressing Arg1",
    "- Dot colour encodes the mean SCT expression level",
    "- Clusters are ordered left to right by overall mean Arg1 expression (highest first)",
    "- Two panels: NO EAE (top) and EAE (bottom)",
    "",
    "### Violin plot",
    "- WT (blue) and KO (red) violins shown side by side per cluster",
    "- Box inside each violin shows median and IQR",
    "- Significance stars above each cluster: *** p<0.001, ** p<0.01, * p<0.05, ns = not significant",
    "- Two panels: NO EAE (top) and EAE (bottom)",
    "- Note: NO EAE panel is nearly empty — Arg1 is not expressed outside of EAE",
    "",
    "### Bar plot",
    "- All 20 clusters shown in numerical order",
    "- Bar height = mean SCT expression, error bars = standard error of the mean",
    "- WT (blue) and KO (red) side by side per cluster",
    "",
    "## Script",
    "",
    "Generated by the corresponding script in `R/`",
    "",
    "```bash",
    "Rscript R/<script>.R",
    "```"
), con = readme_path)
file.copy(readme_path, file.path(PROCESSED_DIR, "README.md"), overwrite = TRUE)
cat(sprintf("Saved: %s\n", readme_path))



cat("\nDone.\n")
