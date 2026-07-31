library(Seurat)
library(Matrix)
library(ggplot2)
library(tibble)
library(dplyr)

project_dir <- "/iss-scratch/CoreBioinformatics/rk720/human_sucnr1"

source(
  file.path(
    project_dir,
    "helper_functions.R"
  )
)

source(
  file.path(
    project_dir,
    "scenic_helper_functions.R"
  )
)

dataset_configs <- list(
  Absinta = list(
    object_file = file.path(
      project_dir,
      "objects",
      "seurat",
      "absinta_microglia_obj.rds"
    ),
    auc_file = file.path(
      project_dir,
      "scenic",
      "results",
      "Absinta",
      "Absinta_regulon_auc.csv"
    ),
    dataset_name = "Absinta",
    analysis_dir = file.path(
      project_dir,
      "scenic",
      "analysis_v3",
      "Absinta"
    ),
    output_object_file = file.path(
      project_dir,
      "objects",
      "seurat",
      "absinta_microglia_obj_with_scenic.rds"
    )
  ),

  Schirmer = list(
    object_file = file.path(
      project_dir,
      "objects",
      "seurat",
      "schirmer_microglia_obj.rds"
    ),
    auc_file = file.path(
      project_dir,
      "scenic",
      "results",
      "Schirmer",
      "Schirmer_regulon_auc.csv"
    ),
    dataset_name = "Schirmer",
    analysis_dir = file.path(
      project_dir,
      "scenic",
      "analysis_v3",
      "Schirmer"
    ),
    output_object_file = file.path(
      project_dir,
      "objects",
      "seurat",
      "schirmer_microglia_obj_with_scenic.rds"
    )
  ),

  MacNair = list(
    object_file = file.path(
      project_dir,
      "objects",
      "seurat",
      "macnair_microglia_obj.rds"
    ),
    auc_file = file.path(
      project_dir,
      "scenic",
      "results",
      "MacNair",
      "MacNair_regulon_auc.csv"
    ),
    dataset_name = "MacNair",
    analysis_dir = file.path(
      project_dir,
      "scenic",
      "analysis_v3",
      "MacNair"
    ),
    output_object_file = file.path(
      project_dir,
      "objects",
      "seurat",
      "macnair_microglia_obj_with_scenic.rds"
    )
  )
)

datasets_to_run <- names(dataset_configs)
# datasets_to_run <- c("Schirmer")
scenic_analysis_results <- lapply(
  datasets_to_run,
  function(dataset_name) {
    config <- dataset_configs[[dataset_name]]

    do.call(
      analyse_scenic_dataset,
      c(
        config,
        list(
          gene = "SUCNR1",
          scenic_assay = "SCENIC",
          top_n_mean_auc = 10
        )
      )
    )
  }
)

names(scenic_analysis_results) <- datasets_to_run


