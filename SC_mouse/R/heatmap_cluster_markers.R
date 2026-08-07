library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(grid)

options(future.globals.maxSize = 4 * 1024^3)

# ============================================================
# Heatmap of top cluster marker genes (cluster vs complement)
#
# 1. Load top-N marker genes per cluster from deg_cluster_markers output
# 2. Compute mean (or sum) SCT expression per gene per cluster
# 3. Z-score each gene across clusters
# 4. Plot heatmap with diagonal ordering
#
# Output: PDF + PNG heatmaps per EAE filter level.
# ============================================================

# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

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

MARKERS_DIR <- file.path(OUT_DIR, "deg_cluster_markers")
OUTPUT_DIR <- ensure_out_dir(file.path("deg_cluster_markers", "heatmaps"))
PROCESSED_DIR <- OUTPUT_DIR

# Top-N genes CSV from deg_cluster_markers.R
TOP_GENES_PATH <- file.path(MARKERS_DIR, "markers_top_genes.csv")

# Clustering column
CLUSTER_COL <- "stable_20_clusters"

# Aggregation function: "mean" (default, standard) or "sum"
AGGREGATION_FN <- "mean"

# Top N per cluster (must match what was used in deg_cluster_markers.R)
TOP_N <- 5

# EAE filter levels to plot
EAE_FILTERS <- c("ALL", "EAE", "NO EAE")

# Z-score clamp range (values beyond this are clamped for plotting)
Z_CLAMP <- 2

# Plot dimensions
PLOT_WIDTH_PER_GENE <- 0.3    # inches per gene on x-axis
PLOT_HEIGHT <- 10              # inches
PLOT_MIN_WIDTH <- 12

# Sanitise filter name for filenames
sanitise_filter <- function(filter_name) gsub(" ", "", filter_name)

# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------

source_methods("load_annotated_object.R")
cat("Loading annotated Seurat object...\n")
so <- load_annotated_object(add_marker_celltypes = FALSE)
cat(sprintf("  %d cells x %d features\n", ncol(so), nrow(so)))
gc()

# ------------------------------------------------------------
# Helper: compute aggregated expression matrix (clusters x genes)
# ------------------------------------------------------------

compute_cluster_expression <- function(so_sub, genes, cluster_col, agg_fn) {
    # Get SCT data matrix for selected genes
    expr <- GetAssayData(so_sub, assay = "SCT", layer = "data")
    genes_present <- intersect(genes, rownames(expr))

    if (length(genes_present) < length(genes)) {
        missing <- setdiff(genes, genes_present)
        cat(sprintf("    WARNING: %d genes not found in SCT assay: %s\n",
            length(missing), paste(head(missing, 5), collapse = ", ")))
    }

    expr <- expr[genes_present, , drop = FALSE]
    clusters <- so_sub[[cluster_col, drop = TRUE]]

    # Compute mean or sum per cluster
    agg_func <- if (agg_fn == "sum") Matrix::colSums else Matrix::colMeans
    cluster_ids <- sort(unique(clusters))

    result <- sapply(cluster_ids, function(cl) {
        cells <- which(clusters == cl)
        if (length(cells) == 0) return(rep(NA, length(genes_present)))
        sub_expr <- expr[, cells, drop = FALSE]
        if (agg_fn == "sum") {
            Matrix::rowSums(sub_expr)
        } else {
            Matrix::rowMeans(sub_expr)
        }
    })

    # Ensure result is always a matrix (sapply can simplify to vector)
    if (!is.matrix(result)) {
        result <- matrix(result, nrow = length(genes_present),
                         ncol = length(cluster_ids))
    }
    colnames(result) <- cluster_ids
    rownames(result) <- genes_present
    result
}

# ------------------------------------------------------------
# Helper: z-score each gene (row) across clusters (columns)
# ------------------------------------------------------------

zscore_rows <- function(mat) {
    t(apply(mat, 1, function(x) {
        s <- sd(x, na.rm = TRUE)
        if (is.na(s) || s == 0) return(rep(0, length(x)))
        (x - mean(x, na.rm = TRUE)) / s
    }))
}

# ------------------------------------------------------------
# Helper: order genes for diagonal pattern
# ------------------------------------------------------------

order_genes_diagonal <- function(zscore_mat, gene_cluster_map, cluster_order) {
    # gene_cluster_map: named vector, gene -> owning cluster
    # cluster_order: desired cluster order (top to bottom on y-axis)
    # We want genes ordered left-to-right matching clusters top-to-bottom

    ordered_genes <- c()
    for (cl in cluster_order) {
        cl_genes <- names(gene_cluster_map[gene_cluster_map == cl])
        cl_genes <- intersect(cl_genes, rownames(zscore_mat))
        if (length(cl_genes) == 0) next

        if (length(cl_genes) > 1) {
            # Sort by z-score in this cluster (descending)
            cl_col <- as.character(cl)
            if (cl_col %in% colnames(zscore_mat)) {
                cl_zscores <- zscore_mat[cl_genes, cl_col]
                cl_genes <- cl_genes[order(-cl_zscores)]
            }
        }
        ordered_genes <- c(ordered_genes, cl_genes)
    }
    ordered_genes
}

# ------------------------------------------------------------
# Helper: plot heatmap
# ------------------------------------------------------------

plot_heatmap <- function(zscore_mat, gene_order, cluster_order,
                         filter_label, z_clamp, agg_fn) {
    # Clamp z-scores
    zscore_clamped <- pmin(pmax(zscore_mat, -z_clamp), z_clamp)

    # Convert to long format
    df <- as.data.frame(zscore_clamped)
    df$gene <- rownames(df)
    df_long <- pivot_longer(df, cols = -gene, names_to = "cluster", values_to = "value")

    # Set factor levels for ordering
    df_long$gene <- factor(df_long$gene, levels = gene_order)
    df_long$cluster <- factor(df_long$cluster, levels = rev(cluster_order))

    p <- ggplot(df_long, aes(x = gene, y = cluster, fill = value)) +
        geom_tile(color = NA) +
        scale_fill_gradientn(
            colors = c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
                       "#F7F7F7",
                       "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"),
            limits = c(-z_clamp, z_clamp),
            name = "value"
        ) +
        labs(
            x = NULL,
            y = "Clusters",
            title = NULL
        ) +
        theme_minimal(base_size = 18) +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 18,
                                       face = "italic"),
            axis.text.y = element_text(size = 22, face = "bold"),
            axis.title.y = element_text(size = 24, face = "bold"),
            axis.ticks = element_blank(),
            panel.grid = element_blank(),
            panel.background = element_rect(fill = "white", color = NA),
            plot.background = element_rect(fill = "white", color = NA),
            legend.position = "right",
            legend.text = element_text(size = 18),
            legend.title = element_text(size = 20),
            legend.key.height = unit(3, "cm"),
            plot.margin = margin(15, 15, 15, 15)
        )

    p
}

# ============================================================
# Main loop: one heatmap per EAE filter
# ============================================================

for (eae_filter in EAE_FILTERS) {
    filter_tag <- sanitise_filter(eae_filter)

    cat(sprintf("\n##########################################################\n"))
    cat(sprintf("# Heatmap: %s (aggregation: %s)\n", eae_filter, AGGREGATION_FN))
    cat(sprintf("##########################################################\n"))

    # Get top genes for this filter
    filter_genes <- top_genes[top_genes$filter == eae_filter, ]
    if (nrow(filter_genes) == 0) {
        cat("  No top genes found for this filter, skipping.\n")
        next
    }

    # Build gene -> owning cluster map
    # If a gene appears in multiple clusters, assign to highest avg_log2FC
    gene_cluster_map <- filter_genes %>%
        group_by(gene) %>%
        slice_max(avg_log2FC, n = 1, with_ties = FALSE) %>%
        ungroup()
    gene_map <- setNames(as.character(gene_cluster_map$cluster),
                         gene_cluster_map$gene)

    gene_universe <- unique(filter_genes$gene)
    cat(sprintf("  Gene universe: %d unique genes from %d clusters\n",
        length(gene_universe), length(unique(filter_genes$cluster))))

    # Subset Seurat object by EAE filter
    if (eae_filter == "ALL") {
        so_sub <- so
    } else {
        so_sub <- subset(so, eae == eae_filter)
    }
    cat(sprintf("  Cells: %d\n", ncol(so_sub)))

    # Compute expression matrix
    cat(sprintf("  Computing %s expression per cluster...\n", AGGREGATION_FN))
    expr_mat <- compute_cluster_expression(so_sub, gene_universe, CLUSTER_COL,
                                            AGGREGATION_FN)

    # Z-score across clusters
    cat("  Z-scoring across clusters...\n")
    z_mat <- zscore_rows(expr_mat)

    # Determine cluster order (descending numeric)
    cluster_ids <- colnames(z_mat)
    cluster_nums <- suppressWarnings(as.numeric(cluster_ids))
    if (!any(is.na(cluster_nums))) {
        cluster_order <- as.character(sort(cluster_nums, decreasing = TRUE))
    } else {
        cluster_order <- sort(cluster_ids, decreasing = TRUE)
    }

    # Order genes for diagonal
    gene_order <- order_genes_diagonal(z_mat, gene_map, cluster_order)
    cat(sprintf("  Gene order: %d genes arranged for diagonal pattern\n",
        length(gene_order)))

    # Plot
    cat("  Plotting...\n")
    p <- plot_heatmap(z_mat, gene_order, cluster_order,
                      eae_filter, Z_CLAMP, AGGREGATION_FN)

    # Save
    plot_width <- max(PLOT_MIN_WIDTH, length(gene_order) * PLOT_WIDTH_PER_GENE + 3)

    fname_pdf <- sprintf("heatmap_markers_%s.pdf", filter_tag)
    fname_png <- sprintf("heatmap_markers_%s.png", filter_tag)

    pdf_path <- file.path(OUTPUT_DIR, fname_pdf)
    ggsave(pdf_path, p, width = plot_width, height = PLOT_HEIGHT,
           device = cairo_pdf, limitsize = FALSE)
    cat(sprintf("  Saved: %s (%.0f x %d in)\n", basename(pdf_path),
        plot_width, PLOT_HEIGHT))

    png_path <- file.path(OUTPUT_DIR, fname_png)
    ggsave(png_path, p, width = plot_width, height = PLOT_HEIGHT, dpi = 300,
           limitsize = FALSE)
    cat(sprintf("  Saved: %s\n", basename(png_path)))

    if (!identical(normalizePath(PROCESSED_DIR), normalizePath(OUTPUT_DIR))) {
        file.copy(pdf_path, file.path(PROCESSED_DIR, fname_pdf), overwrite = TRUE)
        file.copy(png_path, file.path(PROCESSED_DIR, fname_png), overwrite = TRUE)
    }

    if (eae_filter != "ALL") {
        rm(so_sub)
        gc()
    }
}

rm(so)
gc()



cat("\nDone.\n")
