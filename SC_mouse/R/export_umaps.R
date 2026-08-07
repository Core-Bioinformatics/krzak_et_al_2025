library(Seurat)
library(ClustAssess)
library(dittoSeq)
library(ggplot2)
library(grid)

options(future.globals.maxSize = 4 * 1024^3)

# ============================================================
# Export standardized UMAP PDFs with fixed panel geometry.
#
# Why this script exists:
# - keep all exported UMAPs on the same embedding limits
# - reserve a fixed legend column so the panel size does not shrink
#   when a legend has more entries
# - provide a reproducible export path outside interactive notebooks
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

OUTPUT_DIR <- ensure_out_dir("umaps")
PROCESSED_DIR <- OUTPUT_DIR

CLUSTASSESS_PATH <- file.path(DATA_DIR, CLUSTASSESS_RDS)
SEURAT_PATH <- file.path(DATA_DIR, SEURAT_RDS)

FEATURE_TYPE <- "Most_Abundant"
FEATURE_SIZE <- 1950
CLUSTERING_METHOD <- "SLM"
N_CLUSTERS <- 20
CLUSTER_COL <- paste0("stable_", N_CLUSTERS, "_clusters")

PDF_WIDTH <- 20
PDF_HEIGHT <- 20
LEGEND_WIDTH_CM <- 6.5
LEGEND_GAP_CM <- 0.5
BARPLOT_CLUSTER_WIDTH <- 6
BARPLOT_CLUSTER_HEIGHT <- 10
BARPLOT_CELLTYPE_BY_CLUSTER_WIDTH <- 10
BARPLOT_CELLTYPE_BY_CLUSTER_HEIGHT <- 10
BARPLOT_CELLTYPE_BY_GROUP_WIDTH <- 6
BARPLOT_CELLTYPE_BY_GROUP_HEIGHT <- 10

make_legacy_cluster_palette <- function(cluster_levels, base_palette = umap21_fixed) {
  cluster_levels <- sort_feature_levels(as.character(cluster_levels))
  alpha_levels <- sort(cluster_levels)

  stopifnot(length(alpha_levels) <= length(base_palette))

  alpha_palette <- base_palette[seq_along(alpha_levels)]
  names(alpha_palette) <- alpha_levels

  alpha_palette[cluster_levels]
}

source_methods("umap_plotter.R")
source_methods("marker_genes.R")
if (!file.exists(SEURAT_PATH)) stop("Missing Seurat RDS in DATA_DIR: ", SEURAT_PATH)
if (!file.exists(CLUSTASSESS_PATH)) stop("Missing ClustAssess RDS in DATA_DIR: ", CLUSTASSESS_PATH)


cat("Loading cached objects...\n")
test_automm <- readRDS(CLUSTASSESS_PATH)
so_subset_clean <- readRDS(SEURAT_PATH)

chosen_clusters <- get_clusters_from_clustassess_object(
  test_automm,
  feature_type = FEATURE_TYPE,
  feature_size = FEATURE_SIZE,
  clustering_method = CLUSTERING_METHOD,
  nclusters = N_CLUSTERS
)

so_subset_clean[[CLUSTER_COL]] <- factor(
  chosen_clusters[[as.character(N_CLUSTERS)]]$partitions[[1]]$mb,
  levels = seq_len(N_CLUSTERS)
)

umap_mat <- test_automm[[FEATURE_TYPE]][[as.character(FEATURE_SIZE)]]$umap
stopifnot(nrow(umap_mat) == ncol(so_subset_clean))
colnames(umap_mat) <- c("umap_1", "umap_2")
so_subset_clean@reductions$umap@cell.embeddings <- umap_mat
Idents(so_subset_clean) <- CLUSTER_COL

so_subset_clean$eae_condition <- ifelse(
  is.na(so_subset_clean$eae) | is.na(so_subset_clean$condition),
  NA_character_,
  paste0(so_subset_clean$eae, "_", so_subset_clean$condition)
)
so_subset_clean$eae_condition <- factor(so_subset_clean$eae_condition)

cat("Adding celltype marker calls...\n")
so_subset_clean <- add_celltype_markers(
  so_subset_clean,
  assay = "SCT",
  layer = "data",
  default_min_ratio = 0.3
)

fixed_limits <- extract_umap_limits(so_subset_clean, embedding_name = "umap")
assigned_celltype_palette <- so_subset_clean@misc$celltype_markers_palette
assigned_celltype_palette <- assigned_celltype_palette[
  names(assigned_celltype_palette) != "Unassigned"
]

cl_chr <- as.character(so_subset_clean[[]][[CLUSTER_COL]])
lvls_num <- sort(unique(as.integer(cl_chr)))
lvls <- as.character(lvls_num)

so_subset_clean[[paste0(CLUSTER_COL, "_ordered")]] <- factor(cl_chr, levels = lvls)

stopifnot(length(lvls) <= length(umap21_fixed))
cluster_palette <- make_legacy_cluster_palette(lvls)
stopifnot(identical(names(cluster_palette), lvls))

so_ct_plot <- subset(
  so_subset_clean,
  subset = !is.na(celltype_markers) & celltype_markers != "Unassigned"
)

so_ct_plot$celltype_markers <- factor(
  as.character(so_ct_plot$celltype_markers),
  levels = names(assigned_celltype_palette)
)

x_levels_alpha <- levels(as.factor(as.character(so_subset_clean[[CLUSTER_COL, drop = TRUE]])))
x_reorder <- order(as.integer(x_levels_alpha))

plot_specs <- list(
  list(
    file = "umap_stable_clusters.pdf",
    feature = CLUSTER_COL,
    palette_colors = cluster_palette,
    plot_title = "UMAP: Stable Clusters",
    legend_title = "Cluster"
  ),
  list(
    file = "umap_condition.pdf",
    feature = "condition",
    palette_colors = c("ctrl" = "#0C7FB0", "ko" = "#D81B60"),
    plot_title = "UMAP: Condition",
    legend_title = "Condition"
  ),
  list(
    file = "umap_condition_no_eae.pdf",
    feature = "condition",
    filters = list(eae = c("NO EAE")),
    palette_colors = c("ctrl" = "#0C7FB0", "ko" = "#D81B60"),
    plot_title = "UMAP: Condition - NO EAE",
    legend_title = "Condition"
  ),
  list(
    file = "umap_condition_eae.pdf",
    feature = "condition",
    filters = list(eae = c("EAE")),
    palette_colors = c("ctrl" = "#0C7FB0", "ko" = "#D81B60"),
    plot_title = "UMAP: Condition - EAE",
    legend_title = "Condition"
  ),
  list(
    file = "umap_celltype_markers.pdf",
    feature = "celltype_markers",
    filters = list(celltype_markers = names(assigned_celltype_palette)),
    palette_colors = assigned_celltype_palette,
    plot_title = "UMAP: Cell Type Markers",
    legend_title = "Cell type"
  )
)

for (spec in plot_specs) {
  out_path <- file.path(OUTPUT_DIR, spec$file)
  cat(sprintf("Exporting %s\n", spec$file))

  plot_args <- list(
    seurat_obj = so_subset_clean,
    feature = spec$feature,
    filters = if (!is.null(spec$filters)) spec$filters else list(),
    palette_colors = spec$palette_colors,
    palette_mode = if (!is.null(spec$palette_mode)) spec$palette_mode else "auto",
    add_title = TRUE,
    plot_title = spec$plot_title,
    legend_title = spec$legend_title,
    alpha = 0.7,
    point_size = 0.8,
    legend_text_size = 16,
    legend_key_point_size = 7,
    legend_title_size = 18,
    plot_title_size = 18,
    pdf_file_out = out_path,
    pdf_width = PDF_WIDTH,
    pdf_height = PDF_HEIGHT,
    pdf_autoscale_style_sizes = FALSE,
    fixed_limits = fixed_limits,
    legend_width_cm = LEGEND_WIDTH_CM,
    legend_gap_cm = LEGEND_GAP_CM
  )

  do.call(plot_umap_metadata, plot_args)
  file.copy(out_path, file.path(PROCESSED_DIR, spec$file), overwrite = TRUE)
}

barplot_cluster_path <- file.path(OUTPUT_DIR, "barplot_clusters_by_eae_condition.pdf")
cat(sprintf("Exporting %s\n", basename(barplot_cluster_path)))
grDevices::pdf(
  file = barplot_cluster_path,
  width = BARPLOT_CLUSTER_WIDTH,
  height = BARPLOT_CLUSTER_HEIGHT,
  useDingbats = FALSE
)
print(
  dittoSeq::dittoBarPlot(
    so_subset_clean,
    paste0(CLUSTER_COL, "_ordered"),
    group.by = "eae_condition",
    main = "Distribution of Stable Clusters by EAE/Condition",
    ylab = "Fraction of cells",
    xlab = "EAE & Condition",
    legend.title = "Cluster",
    retain.factor.levels = TRUE,
    var.labels.reorder = seq_along(lvls),
    color.panel = cluster_palette,
    colors = seq_along(cluster_palette)
  )
)
grDevices::dev.off()
file.copy(
  barplot_cluster_path,
  file.path(PROCESSED_DIR, basename(barplot_cluster_path)),
  overwrite = TRUE
)

barplot_celltype_cluster_path <- file.path(OUTPUT_DIR, "barplot_celltype_by_cluster.pdf")
cat(sprintf("Exporting %s\n", basename(barplot_celltype_cluster_path)))
grDevices::pdf(
  file = barplot_celltype_cluster_path,
  width = BARPLOT_CELLTYPE_BY_CLUSTER_WIDTH,
  height = BARPLOT_CELLTYPE_BY_CLUSTER_HEIGHT,
  useDingbats = FALSE
)
print(
  dittoSeq::dittoBarPlot(
    so_ct_plot,
    var = "celltype_markers",
    group.by = CLUSTER_COL,
    main = "Cell-type composition by stable 20 clusters",
    ylab = "Fraction of cells",
    xlab = "Cluster",
    legend.title = "Cell type",
    color.panel = assigned_celltype_palette,
    x.reorder = x_reorder,
    colors = seq_along(assigned_celltype_palette)
  )
)
grDevices::dev.off()
file.copy(
  barplot_celltype_cluster_path,
  file.path(PROCESSED_DIR, basename(barplot_celltype_cluster_path)),
  overwrite = TRUE
)

barplot_celltype_group_path <- file.path(OUTPUT_DIR, "barplot_celltype_by_eae_condition.pdf")
cat(sprintf("Exporting %s\n", basename(barplot_celltype_group_path)))
grDevices::pdf(
  file = barplot_celltype_group_path,
  width = BARPLOT_CELLTYPE_BY_GROUP_WIDTH,
  height = BARPLOT_CELLTYPE_BY_GROUP_HEIGHT,
  useDingbats = FALSE
)
print(
  dittoSeq::dittoBarPlot(
    so_ct_plot,
    var = "celltype_markers",
    group.by = "eae_condition",
    main = "Cell-type composition by EAE & Condition",
    ylab = "Fraction of cells",
    xlab = "EAE & Condition",
    legend.title = "Cell type",
    color.panel = assigned_celltype_palette,
    colors = seq_along(assigned_celltype_palette)
  )
)
grDevices::dev.off()
file.copy(
  barplot_celltype_group_path,
  file.path(PROCESSED_DIR, basename(barplot_celltype_group_path)),
  overwrite = TRUE
)

cat(sprintf("\nCopied PDFs to %s\n", PROCESSED_DIR))


cat("\nDone.\n")
