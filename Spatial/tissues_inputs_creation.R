library(Matrix)
library(FNN)
library(ggplot2)
library(matrixStats)
library(dplyr)
library(Seurat)
library(qs2)
library(SeuratWrappers)

set.seed(42)

setwd("/Users/rafael/Desktop/RA/krzak_et_al_2025/Spatial")
source("spatial_helper_functions.R")

seu <- qs_read("/Users/rafael/Desktop/RA/MS/sharing_ms_st_for_flufftail/spatial_seurat_wm_combined.qs2")
tissues_alsema <- unique(seu$Image[seu$study_name == "alsema2024" & !is.na(seu$Image)])

preprocess_tissue <- function(seu, tissue, num_genes)  {
  seu.sub.tissue <- subset(seu, subset = Image == !!tissue)

  hvgs <- VariableFeatures(
    SCTransform(
      seu.sub.tissue,
      assay = "Spatial",
      return.only.var.genes = TRUE,
      verbose = FALSE
    )
  )

  # Get combined feature set
  feats <- union(hvgs[1:num_genes], "SUCNR1")

  seu.sub.tissue <- SCTransform(
    seu.sub.tissue,
    assay = "Spatial",
    residual.features = feats,
    verbose = T
  )

  seu.sub.tissue <- RunPCA(seu.sub.tissue, assay = "SCT", verbose  = FALSE)

  seu.sub.tissue <- RunUMAP(
    seu.sub.tissue,
    reduction = "pca",
    dims      = 1:30,
    assay     = "SCT",
    verbose   = FALSE
  )

  umap_plot <- DimPlot(seu.sub.tissue, group.by = "niches_detailed")
#   print(umap_plot)

  return(seu.sub.tissue)
}


# Preprocess tissues (full)
output_dir_preprocessed <- "objects/preprocessed_tissues"
dir.create(output_dir_preprocessed, recursive = TRUE, showWarnings = FALSE)

tissues_preprocessed <- list()

for (tissue in tissues_alsema) {
  message("Processing: ", tissue)
  seu.sub <- preprocess_tissue(seu, tissue, 500)
  tissues_preprocessed[[tissue]] <- seu.sub
  
  saveRDS(seu.sub, file = file.path(output_dir_preprocessed, paste0("preprocessed_", tissue, ".rds")))
  rm(seu.sub)
  gc()
}

# Creating minimum files needed for further analysis.
tissues_inputs <- list()

for (tissue in names(tissues_preprocessed)) {
  print(tissue)
  seu.sub <- tissues_preprocessed[[tissue]]
  expr_mat <- GetAssayData(seu.sub, assay = "SCT", slot = "data")
  coords <- get_coords(seu.sub)
  metadata <- seu.sub@meta.data
  tissues_inputs[[tissue]] <- list(expr_mat = expr_mat, coords = coords, metadata = metadata)
}
saveRDS(tissues_inputs, "objects/tissues_inputs.rds")
