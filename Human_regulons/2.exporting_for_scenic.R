library(Seurat)
library(SCopeLoomR)
library(Matrix)

project_dir <- "/iss-scratch/CoreBioinformatics/rk720/human_sucnr1"

source(file.path(project_dir, "helper_functions.R"))

objects_dir <- file.path(project_dir, "objects","seurat")

scenic_dir <- file.path(project_dir, "scenic")

# Absinta
absinta_microglia_obj <- readRDS(file.path(objects_dir,"absinta_microglia_obj.rds"))

absinta_umap_name <- grep(
  "^umap_ma",
  Reductions(absinta_microglia_obj),
  value = TRUE
)

absinta_export_result <- export_pyscenic_input(
  obj = absinta_microglia_obj,
  dataset_name = "Absinta",
  output_root = scenic_dir,
  reduction_name = absinta_umap_name,
  assay_use = "RNA",
  min_cell_fraction = 0.01,
  overwrite = TRUE
)

print(absinta_export_result$summary)


# Schirmer
schirmer_microglia_obj <- readRDS(file.path(objects_dir,"schirmer_microglia_obj.rds"))

schirmer_umap_name <- grep(
  "^umap_ma",
  Reductions(schirmer_microglia_obj),
  value = TRUE
)

schirmer_export_result <- export_pyscenic_input(
  obj = schirmer_microglia_obj,
  dataset_name = "Schirmer",
  output_root = scenic_dir,
  reduction_name = schirmer_umap_name,
  assay_use = "RNA",
  min_cell_fraction = 0.01,
  overwrite = TRUE
)

print(schirmer_export_result$summary)


# MacNair
macnair_microglia_obj <- readRDS(file.path(objects_dir,"macnair_microglia_obj.rds"))

macnair_umap_name <- grep(
  "^umap_ma",
  Reductions(macnair_microglia_obj),
  value = TRUE
)

macnair_export_result <- export_pyscenic_input(
  obj = macnair_microglia_obj,
  dataset_name = "MacNair",
  output_root = scenic_dir,
  reduction_name = macnair_umap_name,
  assay_use = "RNA",
  min_cell_fraction = 0.01,
  overwrite = TRUE
)

print(macnair_export_result$summary)