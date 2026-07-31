Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(1234)

library(BiocParallel)
BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)

if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
}


library(flowCore)
library(CATALYST)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(tidyr)
library(lsa) 
library(miloR)
library(scater)
library(ggrepel)

setwd("/Users/macbookpro/Desktop/Desktop_rafa/RA/cytof")
source("cytof_analysis/cytof_helper_functions.R")
dir.create("OUTPUTS_RAFA_v3/")

analysis_sets <- list(
  # full = character(0)
  # no_foxp3 = c("Foxp3"),
  no_foxp3_p2ry12 = c("Foxp3", "P2RY12")
)


run_basic_cytof_dataset <- function(
    folder_path,
    dataset_name,
    dataset_tag,
    annotation_map,
    cols_map,
    analysis_tag,
    exclude_markers = character(0)
) {
  
  pdf_infix_name <- paste(dataset_tag, analysis_tag, sep = "_")
  
  sce <- analyze_cytof_dataset(
    folder_path = folder_path,
    dataset_name = dataset_name,
    pdf_infix_name = pdf_infix_name,
    exclude_markers = exclude_markers,
    output_gene_expression=FALSE
  )
  
  save_plot_by_condition(
    sce,
    meta_key = NULL,
    facet_col = "condition",
    pdf_infix_name = pdf_infix_name
  )
  
  sce <- clustering_cytof_dataset(
    sce,
    pdf_infix_name = pdf_infix_name,
    maxK = 40
  )
  
  save_plot_by_condition(
    sce,
    meta_key = "meta22",
    facet_col = "condition",
    pdf_infix_name = pdf_infix_name
  )
  
  calculate_markers_distributions(sce, pdf_infix_name)
  
  # umaps gene expression for supplementary 
  save_gene_expression_umaps_supplementary(
    sce,
    pdf_infix_name = pdf_infix_name,
    point_size = 0.55
  )
  
  return(sce)
}


run_annotation_and_downstream <- function(
    sce,
    dataset_name,
    dataset_tag,
    annotation_map,
    cols_map,
    analysis_tag,
    exclude_markers = character(0),
    bcell_gate = NULL
) {
  
  pdf_infix_name <- paste(dataset_tag, analysis_tag, sep = "_")
  
  sce <- annotate_based_on_maps(
    sce,
    meta_key = "meta22",
    pdf_infix_name = pdf_infix_name,
    annotation_map = annotation_map,
    cols_map = cols_map
  )
  
  celltype_annotation_col <- "meta_cory_annotation"
  celltype_cols_map <- cols_map
  
  if (!is.null(bcell_gate)) {
    
    sce <- add_bcell_umap_gate_annotation(
      sce = sce,
      pdf_infix_name = pdf_infix_name,
      source_col = "meta_cory_annotation",
      output_col = "meta_cory_annotation_b_cells_included",
      gate_col = "manual_bcell_gate",
      gate_type = bcell_gate$gate_type,
      x_center = bcell_gate$x_center,
      y_center = bcell_gate$y_center,
      x_radius = bcell_gate$x_radius,
      y_radius = bcell_gate$y_radius,
      x_min = bcell_gate$x_min,
      x_max = bcell_gate$x_max,
      y_min = bcell_gate$y_min,
      y_max = bcell_gate$y_max,
      cols_map = cols_map,
      point_size = annotation_umap_point_size
    )
    
    celltype_annotation_col <- "meta_cory_annotation_b_cells_included"
    celltype_cols_map <- c(cols_map, "B cells" = all_cory_cols[["B cells"]])
  }
  
  marker_medians_by_condition <- plot_marker_heatmap_by_condition(
    sce,
    pdf_infix_name = pdf_infix_name,
    annotation_col = celltype_annotation_col,
    folder_save = "heatmaps_cell_types"
  )
  
  delta_marker_medians_by_condition <- plot_marker_delta_heatmap_by_condition(
    sce,
    pdf_infix_name = pdf_infix_name,
    annotation_col = celltype_annotation_col,
    folder_save = "delta_heatmaps_cell_types"
  )
  
  sce$meta22_cluster <- factor(
    as.character(cluster_ids(sce, k = "meta22")),
    levels = as.character(seq_len(22))
  )

  save_cluster_umap_by_condition(
    sce,
    pdf_infix_name = pdf_infix_name,
    cluster_col = "meta22_cluster",
    condition_col = "condition"
  )

  save_proportion_plots_by_condition(
    sce,
    pdf_infix_name = pdf_infix_name,
    cluster_col = "meta22_cluster",
    celltype_col = celltype_annotation_col,
    condition_col = "condition",
    celltype_cols = celltype_cols_map
  )
  
  marker_medians_by_condition_meta22 <- plot_marker_heatmap_by_condition(
    sce,
    pdf_infix_name = pdf_infix_name,
    annotation_col = "meta22_cluster",
    folder_save = "heatmaps_clusters"
  )
  
  delta_marker_medians_by_condition_meta22 <- plot_marker_delta_heatmap_by_condition(
    sce,
    pdf_infix_name = pdf_infix_name,
    annotation_col = "meta22_cluster",
    folder_save = "delta_heatmaps_clusters"
  )
  
  run_statistical_tests(
    sce,
    pdf_infix_name = pdf_infix_name,
    cluster_col = "meta22_cluster"
  )
  
  export_per_replicate_marker_expression(
    sce,
    pdf_infix_name = pdf_infix_name,
    annotation_col = celltype_annotation_col
  )
  
  export_per_replicate_marker_expression(
    sce,
    pdf_infix_name = pdf_infix_name,
    annotation_col = "meta22_cluster"
  )
  
  return(sce)
}

##### Initial basic analysis ####
all_results <- list()

for (analysis_tag in names(analysis_sets)) {

  exclude_markers <- analysis_sets[[analysis_tag]]

  all_results[[analysis_tag]] <- list(

    global = run_basic_cytof_dataset(
      folder_path = "global_sucnr1_ko_C-EAE_gated_exported",
      dataset_name = "Global KO",
      dataset_tag = "global",
      annotation_map = global_annotations,
      cols_map = global_cols,
      analysis_tag = analysis_tag,
      exclude_markers = exclude_markers
    ),

    chronic = run_basic_cytof_dataset(
      folder_path = "microglia_ko_Chronic_gated_exported",
      dataset_name = "Microglia KO Chronic",
      dataset_tag = "chronic",
      annotation_map = chronic_annotations,
      cols_map = chronic_cols,
      analysis_tag = analysis_tag,
      exclude_markers = exclude_markers
    ),

    acute = run_basic_cytof_dataset(
      folder_path = "microglia_ko_Acute_gated_exported",
      dataset_name = "Microglia KO Acute",
      dataset_tag = "acute",
      annotation_map = acute_annotations,
      cols_map = acute_cols,
      analysis_tag = analysis_tag,
      exclude_markers = exclude_markers
    )
  )
}


# save(
#   all_results,
#   file = "/Volumes/HD_rafael/Cytof/cytof_17_06_no_foxp3_p2ry12_deterministic.RData"
# )

#### Loading and downstream ####
load("/Volumes/HD_rafael/Cytof/cytof_17_06_no_foxp3_p2ry12_deterministic.RData")

sce_global <- all_results$no_foxp3_p2ry12$global
sce_chronic <- all_results$no_foxp3_p2ry12$chronic
sce_acute <- all_results$no_foxp3_p2ry12$acute

# b cells location for chronic and global (weren't on their own in the unsupervised clustering)
chronic_bcell_gate <- list(
  gate_type = "ellipse",
  x_center = -6.75,
  y_center = -5.55,
  x_radius = 0.75,
  y_radius = 0.45
)

global_bcell_gate <- list(
  gate_type = "ellipse",
  x_center = 6.15,
  y_center = -0.20,
  x_radius = 0.55,
  y_radius = 0.30
)


all_annotations_and_downstream_results <- list()

for (analysis_tag in names(analysis_sets)) {
  
  exclude_markers <- analysis_sets[[analysis_tag]]
  
  all_annotations_and_downstream_results[[analysis_tag]] <- list(
    
    global <- run_annotation_and_downstream(sce_global,
                                  dataset_name = "Global KO",
                                  dataset_tag = "global",
                                  annotation_map = global_annotations,
                                  cols_map = global_cols,
                                  analysis_tag = analysis_tag,
                                  exclude_markers = exclude_markers,
                                  bcell_gate=global_bcell_gate),
    
    chronic <- run_annotation_and_downstream(sce_chronic,
                                            dataset_name = "Microglia KO Chronic",
                                            dataset_tag = "chronic",
                                            annotation_map = chronic_annotations,
                                            cols_map = chronic_cols,
                                            analysis_tag = analysis_tag,
                                            exclude_markers = exclude_markers,
                                            bcell_gate=chronic_bcell_gate),
    
    acute <- run_annotation_and_downstream(sce_acute,
                                            dataset_name = "Microglia KO Acute",
                                            dataset_tag = "acute",
                                            annotation_map = acute_annotations,
                                            cols_map = acute_cols,
                                            analysis_tag = analysis_tag,
                                            exclude_markers = exclude_markers)
  )
}


#### Clusters of interest - expression boxplot ####
global_marker_cluster_map <- list(
  "IL-2"   = c(7, 20, 21, 22),
  "TNFa"   = c(7, 20, 21, 22),
  "iNOS"   = c(20, 21, 22),
  "IL-6"   = c(20, 21, 22),
  "IFNg"   = c(20, 21, 22),
  "IL-17A" = c(20, 21, 22)
)

acute_marker_cluster_map <- list(
  "TNFa" = c(3, 4, 16, 19),
  "iNOS" = c(3, 4, 10),
  "IFNg" = c(16),
  "IL-6" = c(2, 3, 4, 11),
  "IL-4" = c(1, 2, 7, 8, 9, 10)
)

chronic_marker_cluster_map <- list(
  "TNFa" = c(14, 19, 21),
  "iNOS" = c(13, 14, 19),
  "IFNg" = c(19),
  "IL-6" = c(19),
  "IL-4" = c(1, 2)
)

plot_marker_boxplot_by_selected_metaclusters(
  sce = sce_global,
  pdf_infix_name = "global_no_foxp3_p2ry12",
  marker_cluster_map = global_marker_cluster_map,
  meta_key = "meta22"
)

plot_marker_violin_by_selected_metaclusters(
  sce = sce_global,
  pdf_infix_name = "global_no_foxp3_p2ry12",
  marker_cluster_map = global_marker_cluster_map,
  meta_key = "meta22"
)

plot_marker_boxplot_by_selected_metaclusters(
  sce = sce_acute,
  pdf_infix_name = "acute_no_foxp3_p2ry12",
  marker_cluster_map = acute_marker_cluster_map,
  meta_key = "meta22"
)

plot_marker_violin_by_selected_metaclusters(
  sce = sce_acute,
  pdf_infix_name = "acute_no_foxp3_p2ry12",
  marker_cluster_map = acute_marker_cluster_map,
  meta_key = "meta22"
)

plot_marker_boxplot_by_selected_metaclusters(
  sce = sce_chronic,
  pdf_infix_name = "chronic_no_foxp3_p2ry12",
  marker_cluster_map = chronic_marker_cluster_map,
  meta_key = "meta22"
)

plot_marker_violin_by_selected_metaclusters(
  sce = sce_acute,
  pdf_infix_name = "acute_no_foxp3_p2ry12",
  marker_cluster_map = acute_marker_cluster_map,
  meta_key = "meta22"
)


#### diffcyt-DA-GLMM abundance tests ####
# load("/Volumes/HD_rafael/Cytof/cytof_20_07_no_foxp3_p2ry12_deterministic_downstream.RData")
# sce_global <- all_annotations_and_downstream_results$no_foxp3_p2ry12[[1]]
# sce_chronic <- all_annotations_and_downstream_results$no_foxp3_p2ry12[[2]]
# sce_acute <- all_annotations_and_downstream_results$no_foxp3_p2ry12[[3]]


sce_global$meta22_cluster <- factor(
  as.character(cluster_ids(sce_global, k = "meta22")),
  levels = as.character(seq_len(22))
)

global_cluster_abundance_diffcyt_glmm <- run_diffcyt_da_glmm(
  sce = sce_global,
  dataset_name = "global_no_foxp3_p2ry12",
  category_col = "meta22_cluster"
)

global_cluster_abundance_diffcyt_glmm
global_celltype_abundance_diffcyt_glmm <- run_diffcyt_da_glmm(
  sce = sce_global,
  dataset_name = "global_no_foxp3_p2ry12",
  category_col = "meta_cory_annotation_b_cells_included"
)
global_celltype_abundance_diffcyt_glmm


# dir.create(
#   "/Volumes/HD_rafael/Cytof/OUTPUTS_RAFA_v3/global_no_foxp3_p2ry12/stats",
#   showWarnings = FALSE,
#   recursive = TRUE
# )

write.csv(
  global_cluster_abundance_diffcyt_glmm,
  "/Volumes/HD_rafael/Cytof/OUTPUTS_RAFA_v3/global_no_foxp3_p2ry12/stats/global_no_foxp3_p2ry12_cluster_abundance_diffcyt_da_glmm_results.csv",
  row.names = FALSE
)

write.csv(
  global_celltype_abundance_diffcyt_glmm,
  "/Volumes/HD_rafael/Cytof/OUTPUTS_RAFA_v3/global_no_foxp3_p2ry12/stats/global_no_foxp3_p2ry12_celltype_abundance_diffcyt_da_glmm_results.csv",
  row.names = FALSE
)


sce_acute$meta22_cluster <- factor(
  as.character(cluster_ids(sce_acute, k = "meta22")),
  levels = as.character(seq_len(22))
)

acute_cluster_abundance_diffcyt_glmm <- run_diffcyt_da_glmm(
  sce = sce_acute,
  dataset_name = "acute_no_foxp3_p2ry12",
  category_col = "meta22_cluster"
)

acute_celltype_abundance_diffcyt_glmm <- run_diffcyt_da_glmm(
  sce = sce_acute,
  dataset_name = "acute_no_foxp3_p2ry12",
  category_col = "meta_cory_annotation"
)

acute_cluster_abundance_diffcyt_glmm
acute_celltype_abundance_diffcyt_glmm

dir.create(
  "/Volumes/HD_rafael/Cytof/OUTPUTS_RAFA_v3/acute_no_foxp3_p2ry12/stats",
  showWarnings = FALSE,
  recursive = TRUE
)

write.csv(
  acute_cluster_abundance_diffcyt_glmm,
  "/Volumes/HD_rafael/Cytof/OUTPUTS_RAFA_v3/acute_no_foxp3_p2ry12/stats/acute_no_foxp3_p2ry12_cluster_abundance_diffcyt_da_glmm_results.csv",
  row.names = FALSE
)

write.csv(
  acute_celltype_abundance_diffcyt_glmm,
  "/Volumes/HD_rafael/Cytof/OUTPUTS_RAFA_v3/acute_no_foxp3_p2ry12/stats/acute_no_foxp3_p2ry12_celltype_abundance_diffcyt_da_glmm_results.csv",
  row.names = FALSE
)


sce_chronic$meta22_cluster <- factor(
  as.character(cluster_ids(sce_chronic, k = "meta22")),
  levels = as.character(seq_len(22))
)

chronic_cluster_abundance_diffcyt_glmm <- run_diffcyt_da_glmm(
  sce = sce_chronic,
  dataset_name = "chronic_no_foxp3_p2ry12",
  category_col = "meta22_cluster"
)

chronic_celltype_abundance_diffcyt_glmm <- run_diffcyt_da_glmm(
  sce = sce_chronic,
  dataset_name = "chronic_no_foxp3_p2ry12",
  category_col = "meta_cory_annotation_b_cells_included"
)

chronic_cluster_abundance_diffcyt_glmm
chronic_celltype_abundance_diffcyt_glmm

dir.create(
  "/Volumes/HD_rafael/Cytof/OUTPUTS_RAFA_v3/chronic_no_foxp3_p2ry12/stats",
  showWarnings = FALSE,
  recursive = TRUE
)

write.csv(
  chronic_cluster_abundance_diffcyt_glmm,
  "/Volumes/HD_rafael/Cytof/OUTPUTS_RAFA_v3/chronic_no_foxp3_p2ry12/stats/chronic_no_foxp3_p2ry12_cluster_abundance_diffcyt_da_glmm_results.csv",
  row.names = FALSE
)

write.csv(
  chronic_celltype_abundance_diffcyt_glmm,
  "/Volumes/HD_rafael/Cytof/OUTPUTS_RAFA_v3/chronic_no_foxp3_p2ry12/stats/chronic_no_foxp3_p2ry12_celltype_abundance_diffcyt_da_glmm_results.csv",
  row.names = FALSE
)


#### Subset cluster delta plots ####
plot_cluster_delta_versions <- function(sce, dataset_tag, markers, versions) {
  sce$meta22_cluster <- factor(
    as.character(cluster_ids(sce, k = "meta22")),
    levels = as.character(seq_len(22))
  )

  results <- lapply(names(versions), function(version) {
    common_args <- list(
      sce = sce,
      pdf_infix_name = paste0(dataset_tag, "_no_foxp3_p2ry12"),
      annotation_col = "meta22_cluster",
      markers = markers,
      annotation_subset = versions[[version]],
      output_suffix = version
    )

    list(
      heatmap = do.call(
        plot_marker_delta_heatmap_by_condition,
        c(common_args, list(folder_save = "delta_heatmaps_clusters_subset"))
      ),
      bubbleplot = do.call(
        plot_marker_delta_bubble_by_condition,
        c(common_args, list(folder_save = "delta_bubbleplots_clusters_subset"))
      )
    )
  })

  setNames(results, names(versions))
}

inflammatory_markers <- c("iNOS", "TNFa", "IFNg", "IL-6", "IL-2", "IL-17A")

acute_delta_versions <- list(
  version_1 = c(3, 10, 15, 2, 4, 6, 7, 8, 9, 12),
  version_2 = c(3, 10, 15, 2, 4, 6, 7, 8, 9, 12, 5),
  version_3 = c(3, 10, 15, 2, 4, 6, 7, 9),
  version_4 = c(3, 10, 15, 2, 4, 6, 7, 9, 5)
)

chronic_delta_versions <- list(
  version_1 = c(19, 1, 14),
  version_2 = c(19, 1, 14, 11, 22)
)

global_delta_versions <- list(
  version_1 = c(1, 2, 3, 14, 9, 12, 4, 6, 7, 15, 18, 16, 17)
)

acute_delta_plots <- plot_cluster_delta_versions(
  sce_acute, "acute", inflammatory_markers, acute_delta_versions
)
chronic_delta_plots <- plot_cluster_delta_versions(
  sce_chronic, "chronic", inflammatory_markers, chronic_delta_versions
)
global_delta_plots <- plot_cluster_delta_versions(
  sce_global, "global", c("CX3CR1", "CCR2"), global_delta_versions
)
