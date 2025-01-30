library(Seurat)
library(dplyr)
library(ggplot2)
library(ClustAssess)

ncores <- 60
my_cluster <- parallel::makeCluster(
  ncores,
  type = "PSOCK"
)

RhpcBLASctl::blas_set_num_threads(1)
doParallel::registerDoParallel(cl = my_cluster)

options(warn=-1)
assay_name <- "SCT"

so_subset = readRDS('objects/so.rds')

n_abundant <- 2000
expr_matrix <- as.matrix(so_subset@assays[[assay_name]]@scale.data)
features <- dimnames(so_subset@assays[[assay_name]])[[1]]
var_features <- so_subset@assays[[assay_name]]@var.features

most_abundant_genes <- rownames(expr_matrix)[order(Matrix::rowSums(expr_matrix), decreasing = TRUE)][1:n_abundant]

gene_list = list(
  "Most_Abundant" = most_abundant_genes,
  "Highly_Variable" = so_subset@assays[[assay_name]]@var.features
)

steps_list = list(
  "Most_Abundant" = seq(from = 500, by = 500, to = n_abundant),
  "Highly_Variable" = seq(from = 500, by = 500, to = n_abundant)
)



test_automm = automatic_stability_assessment(expression_matrix = expr_matrix,
                                             n_repetitions = 20,
                                             temp_file = "clustassess_temp.rds",
                                             n_neigh_sequence = seq(from = 5, to = 50, by = 10),
                                             resolution_sequence = seq(from = 0.1, to = 1.5, by = 0.1),
                                             features_sets = gene_list,
                                             steps = steps_list,
                                             n_top_configs = 2,
                                             npcs = 30,
                                             post_processing = function(pca_emb) {
                                                 set.seed(0)
                                                 harmony::HarmonyMatrix(
                                                   data_mat = pca_emb,
                                                   do_pca = FALSE,
                                                   meta_data = so_subset$rep,
                                                   epsilon.cluster = 1e-07,
                                                   epsilon.harmony = 1e-07,
                                                   verbose = FALSE
                                                 )
                                             },
                                             umap_arguments = list(
                                                min_dist = 0.3,
                                                n_neighbors = 30,
                                                init = "spca",
                                                metric = "cosine"
                                             ))


saveRDS(test_automm, file.path('processed/test_automm.rds'))
on.exit(parallel::stopCluster(cl = my_cluster))

write_shiny_app(
  object = so_subset,
  assay_name = 'SCT',
  clustassess_object = test_automm,
  project_folder = 'processed/clustassess/'
)

