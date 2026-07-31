library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(Matrix)
library(tidyverse)
library(ggplot2)
library(ggnewscale)
library(viridis)

base_dir <- "/iss-corescratch/lp488/monica/luca/data"

source("/iss-scratch/CoreBioinformatics/rk720/human_sucnr1/helper_functions.R")
output_dir <- "/iss-scratch/CoreBioinformatics/rk720/human_sucnr1/expression_umaps_v2"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

genes_friends <- c(
  "ARG1", "FIZZ1", "MGL1", "MGL2", "HIF1A", "IL6", "IL1B", "IL2", "CCR7", "CCR5",
  "NLRP1", "P2RY12", "FUCA1", "KLF4", "TLR7", "HRH1", "MMP12", "CCL7", "MMP1", "CD1E",
  "IL1R2", "CLEC5A", "ADORA1", "AI182371", "AI464131", "AQP11", "DLL1", "ESRP2",
  "FAM131B", "FUZ", "GDPD4", "GGACT", "H1FX", "IGSF11", "KCNA6", "KCTD11", "KLHL25",
  "LCA5L", "LGI3", "MMGT2", "MPP2", "PAQR8", "RGS9BP", "SERPINB8", "SLC16A6", "SORBS2",
  "TXNIP", "PTGS2", "NPAS4", "CYR61", "PHXR4", "KCNQ1OT1", "NFIL3", "TRIM29", "RGS16",
  "FOXD4", "4931406H21RIK", "CALHM2", "KMT2D", "NLRP1C-PS", "IL2", "GM438", "KCNK18",
  "HPCAL4", "OGN", "DLC1", "EGR2", "ZFP42", "ANKRD50", "TMEM132C", "APLNR", "TDRD1"
)

objects_dir <- "/iss-scratch/CoreBioinformatics/rk720/human_sucnr1/objects/seurat"

generate_initial_analysis_plots_full_dataset <- function(seurat_obj, genes_friends, dataset_name, output_dir, cell_type_col = "cell_type", microglia_col = "microglia") {
  dataset_output_dir <- paste0(output_dir, "/", dataset_name)
  print(dataset_output_dir)
  dir.create(dataset_output_dir, recursive = TRUE, showWarnings = FALSE)

  DimPlot(
    seurat_obj,
    reduction = "umap",
    group.by = cell_type_col,
    label = TRUE,
    repel = TRUE,
    pt.size = 0.2
  ) +
    ggtitle(paste0(dataset_name, " original full object UMAP"))

  ggsave(
    file.path(dataset_output_dir, "/original_full_object_cell_type_umap.png"),
    width = 8,
    height = 6,
    dpi = 300
  )

  # sucnr1_full_umap <- plot_gene_positive_full_umap(
  #   obj = seurat_obj,
  #   gene = "SUCNR1",
  #   reduction_name = "umap",
  #   dataset_name = dataset_name,
  #   output_file = file.path(
  #     dataset_output_dir,
  #     "SUCNR1_positive_full_object_umap.pdf"
  #   ),
  #   assay = "RNA",
  #   layer = "counts"
  # )

  sucnr1_full_umap <- plot_gene_positive_full_umap(
    obj = seurat_obj,
    gene = "SUCNR1",
    reduction_name = "umap",
    dataset_name = dataset_name,
    output_file = file.path(
      dataset_output_dir,
      "SUCNR1_positive_full_object_umap.png"
    ),
    assay = "RNA",
    layer = "counts"
  )

  sucnr1_celltype <- plot_sucnr1_proportions_by_celltype(
    obj = seurat_obj,
    assay = "RNA",
    layer = "counts",
    cell_type_col = cell_type_col,
    dataset_name = dataset_name,
    output_prefix = file.path(dataset_output_dir, "/SUCNR1_full_object")
  )

  sucnr1_friends_celltype <- plot_sucnr1_proportions_by_celltype(
    obj = seurat_obj,
    assay = "RNA",
    layer = "counts",
    genes = genes_friends,
    gene_set_name = "sucnr1_friends",
    cell_type_col = cell_type_col,
    dataset_name = dataset_name,
    output_prefix =  file.path(dataset_output_dir, "/SUCNR1_friends_full_object")
  )

  present_friends_genes <- genes_friends[genes_friends %in% rownames(seurat_obj)]

  mean_mat <- get_gene_celltype_mean_matrix(
    obj = seurat_obj,
    genes = present_friends_genes,
    assay = "SCT",
    layer = "data",
    cell_type_col = cell_type_col
  )

  other_cols <- setdiff(colnames(mean_mat), microglia_col)
  gene_scores <- mean_mat[, microglia_col] -
    apply(mean_mat[, other_cols, drop = FALSE], 1, max)

  ranked_genes <- sort(gene_scores, decreasing = TRUE)

  sucnr1_friends_celltype <- plot_sucnr1_proportions_by_celltype(
    obj = seurat_obj,
    genes = names(ranked_genes)[1:10],
    gene_set_name = "ranked_10_sucnr1_friends",
    cell_type_col = cell_type_col,
    dataset_name = dataset_name,
    output_prefix =  file.path(dataset_output_dir, "/SUCNR1_friends_full_object_ranked10")
  )
}

generate_microglia_analysis <- function(
  microglia_result,
  genes_friends,
  dataset_name,
  output_dir,
  object_output_dir,
  assay_use = "SCT",
  layer_use = "data",
  n_ma = 1500
) {
  dataset_output_dir <- file.path(output_dir, dataset_name)

  dir.create(
    dataset_output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  dir.create(
    object_output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  microglia_obj <- microglia_result$object
  reduction_name <- microglia_result$umap_name

  microglia_obj <- add_sucnr1_expression_categories(
    obj = microglia_obj
  )

  sucnr1_raw_plot <- plot_sucnr1_only_umap(
    obj = microglia_obj,
    reduction_name = reduction_name,
    dataset_name = dataset_name,
    output_file = file.path(
      dataset_output_dir,
      paste0("microglia_MA", n_ma, "_SUCNR1_raw_count_umap.png")
    ),
    sucnr1_col = "SUCNR1_raw_count",
    legend_title = "SUCNR1\nraw count"
  )

  sucnr1_sct_plot <- plot_sucnr1_only_umap(
    obj = microglia_obj,
    reduction_name = reduction_name,
    dataset_name = dataset_name,
    output_file = file.path(
      dataset_output_dir,
      paste0("microglia_MA", n_ma, "_SUCNR1_sct_expression_umap.png")
    ),
    sucnr1_col = "SUCNR1_sct_expr",
    legend_title = "SUCNR1\nSCT expression"
  )

  sucnr1_positive_summary <- summary(
    microglia_obj$SUCNR1_expr[
      microglia_obj$SUCNR1_expr > 0
    ]
  )

  print(sucnr1_positive_summary)

  average_expression_plot <- expression_avg_per_niche_plot(
    obj = microglia_obj,
    reduction_name = reduction_name,
    genes = genes_friends,
    dataset_name = paste(dataset_name, "microglia"),
    output_file = file.path(
      dataset_output_dir,
      "microglia_genes_friends_avg_expression_umap.png"
    ),
    assay = "RNA",
    layer = "counts"
  )

  active_gene_count_plot <- expression_active_gene_count_plot(
    obj = microglia_obj,
    reduction_name = reduction_name,
    genes = genes_friends,
    dataset_name = paste(dataset_name, "microglia"),
    output_file = file.path(
      dataset_output_dir,
      "microglia_genes_friends_active_gene_count_umap.png"
    ),
    assay = "RNA",
    layer = "counts"
  )

  object_file_label <- gsub(
    "[^A-Za-z0-9]+",
    "_",
    tolower(dataset_name)
  )

  saveRDS(
    microglia_obj,
    file = file.path(
      object_output_dir,
      paste0(object_file_label, "_microglia_obj.rds")
    )
  )

  return(list(
    object = microglia_obj,
    reduction_name = reduction_name,
    sucnr1_positive_summary = sucnr1_positive_summary,
    sucnr1_raw_plot = sucnr1_raw_plot,
    sucnr1_sct_plot = sucnr1_sct_plot,
    average_expression_plot = average_expression_plot,
    active_gene_count_plot = active_gene_count_plot
  ))
}


#### Absinta ####
absinta_celltypes <- readRDS(file.path(base_dir, "absinta_celltypes.rds"))

generate_initial_analysis_plots_full_dataset(absinta_celltypes, genes_friends, "Absinta", output_dir)

absinta_microglia <- run_microglia_ma_umap(
  obj = absinta_celltypes,
  dataset_name = "Absinta",
  assay_use = "SCT",
  n_ma = 1500,
  cell_type_col = "cell_type",
  microglia_label = "microglia",
  dims_use = 1:30
)

absinta_microglia_analysis <- generate_microglia_analysis(
  microglia_result = absinta_microglia,
  genes_friends = genes_friends,
  dataset_name = "Absinta",
  output_dir = output_dir,
  object_output_dir = objects_dir
)

absinta_microglia_obj <- absinta_microglia_analysis$object


#### Schirmer ####
schirmer_celltypes <- readRDS(file.path(base_dir, "schirmer_celltypes.rds"))

generate_initial_analysis_plots_full_dataset(schirmer_celltypes, genes_friends, "Schirmer", output_dir)

schirmer_microglia <- run_microglia_ma_umap(
  obj = schirmer_celltypes,
  dataset_name = "Schirmer",
  assay_use = "SCT",
  n_ma = 1500,
  cell_type_col = "cell_type",
  microglia_label = "microglia",
  dims_use = 1:30
)

schirmer_microglia_analysis <- generate_microglia_analysis(
  microglia_result = schirmer_microglia,
  genes_friends = genes_friends,
  dataset_name = "Schirmer",
  output_dir = output_dir,
  object_output_dir = objects_dir
)

schirmer_microglia_obj <- schirmer_microglia_analysis$object


#### MacNair ####
library(qs)

mcnair_celltypes <- qread("/servers/sutherland-scratch/andi/projects/0_2505_Luca_MacNair/objects/R/seurat/cleaned_aggr_wm_amsterdam_sct.qs", nthreads = 8)
table(mcnair_celltypes$type_broad)

# does not have PCA and UMAP precomputed.
mcnair_features <- VariableFeatures(
  mcnair_celltypes,
  assay = "SCT"
)

mcnair_celltypes <- RunPCA(
  mcnair_celltypes,
  assay = "SCT",
  features = mcnair_features,
  npcs = 50,
  rev.pca = FALSE,
  weight.by.var = TRUE,
  seed.use = 42,
  reduction.name = "pca",
  reduction.key = "PC_",
  verbose = TRUE
)

mcnair_celltypes <- RunUMAP(
  mcnair_celltypes,
  reduction = "pca",
  dims = 1:50,
  n.neighbors = 30,
  n.components = 2,
  metric = "cosine",
  learning.rate = 1,
  min.dist = 0.3,
  spread = 1,
  set.op.mix.ratio = 1,
  local.connectivity = 1,
  repulsion.strength = 1,
  negative.sample.rate = 5,
  umap.method = "uwot",
  seed.use = 42,
  reduction.name = "umap",
  reduction.key = "UMAP_",
  verbose = TRUE
)

generate_initial_analysis_plots_full_dataset(mcnair_celltypes, genes_friends, "MacNair", output_dir,  cell_type_col = "type_broad", microglia_col = "Microglia")

mcnair_microglia <- run_microglia_ma_umap(
  obj = mcnair_celltypes,
  dataset_name = "MacNair",
  assay_use = "SCT",
  n_ma = 2000,
  cell_type_col = "type_broad",
  microglia_label = "Microglia",
  dims_use = 1:30
)

mcnair_microglia_analysis <- generate_microglia_analysis(
  microglia_result = mcnair_microglia,
  genes_friends = genes_friends,
  dataset_name = "MacNair",
  output_dir = output_dir,
  object_output_dir = objects_dir,
  n_ma = 2000
)

