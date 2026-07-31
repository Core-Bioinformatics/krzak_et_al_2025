library(Seurat)
library(ClustAssess)
library(ggplot2)

seurat_obj <- readRDS(
  paste0(
    "/iss-corescratch/lp488/Luca_June2023/processed/",
    "so_subset_2025-10-sucnr1-noHB-PC2.rds"
  )
)

test_automm <- readRDS(
  paste0(
    "/iss-corescratch/lp488/Luca_June2023/processed/",
    "test_automm_2025-10-sucnr1-noHB-PC2.rds"
  )
)

# Add stable 20-cluster solution
chosen_clusters <- get_clusters_from_clustassess_object(
  test_automm,
  feature_type = "Most_Abundant",
  feature_size = 1950,
  clustering_method = "SLM",
  nclusters = 20
)

seurat_obj$stable_20_clusters <- factor(
  chosen_clusters[["20"]]$partitions[[1]]$mb,
  levels = 1:20
)

# Add UMAP used in the Shiny app
umap_mat <- test_automm[["Most_Abundant"]][["1950"]]$umap

stopifnot(
  nrow(umap_mat) == ncol(seurat_obj)
)

rownames(umap_mat) <- colnames(seurat_obj)
colnames(umap_mat) <- c("UMAP_1", "UMAP_2")

seurat_obj[["umap"]] <- CreateDimReducObject(
  embeddings = umap_mat,
  key = "UMAP_",
  assay = DefaultAssay(seurat_obj)
)

# Recreate marker-voting cell-type annotation
default_min_genes_ratio <- 0.3

gene_list <- list(
  Microglia_Homeostatic = list(
    genes = c(
      "Sparc", "P2ry12", "Siglech", "Tmem119",
      "Sall1", "Hexb", "Csf1r"
    ),
    min_n_genes = 3
  ),

  Microglia_DAM = list(
    genes = c(
      "C1qc", "C1qa", "Cst7", "Ly86",
      "Ccl12", "Spp1", "Lgals3"
    )
  ),

  Monocyte_derived_cells = list(
    genes = c(
      "Lyz2", "H2-Ab1", "Cd74",
      "Ly6c2", "Clec4e", "Ccr2"
    )
  ),

  Macrophages = list(
    genes = c(
      "Ly6c2", "Plac8", "Nos2",
      "Mgst1", "Arg1"
    )
  ),

  Neutrophils = list(
    genes = c(
      "Csf3r", "S100a8", "S100a9"
    )
  ),

  T_cells = list(
    genes = c(
      "Cd3e", "Bcl11b", "Il2ra",
      "Cd247", "Skap1", "Themis"
    )
  ),

  Neurons = list(
    genes = c(
      "Gnb4", "Syt1", "Snap25",
      "Rbfox1", "Sema6d", "Grik2"
    )
  ),

  Astrocytes = list(
    genes = c(
      "Aqp4", "Slc1a2", "Kcnj10", "Atp1a2"
    )
  ),

  Oligodendrocytes = list(
    genes = c(
      "Cnp", "Plp1", "Mbp", "Mobp"
    )
  ),

  Endothelial_Cells = list(
    genes = c(
      "Cldn5", "Ly6c1", "Pltp",
      "Bsg", "Itm2a", "Pecam1"
    )
  ),

  Ependymal_Neural_Stem_Cells = list(
    genes = c(
      "Tuba1a", "Foxj1", "Pifo", "Nnat",
      "Dbi", "Mt3", "Sox2", "Pax6"
    )
  ),

  Pericytes_Fibroblasts = list(
    genes = c(
      "Vtn", "Myl9", "Rgs5",
      "Acta2", "Tagln", "Mgp"
    )
  )
)

cell_type_thresholds <- c(
  Microglia_Homeostatic = 2.5,
  Microglia_DAM = 2.0,
  Monocyte_derived_cells = 2.5,
  Macrophages = 0.0,
  Neutrophils = 0.0,
  T_cells = 0.0,
  Neurons = 0.0,
  Astrocytes = 0.0,
  Oligodendrocytes = 2.5,
  Endothelial_Cells = 0.0,
  Ependymal_Neural_Stem_Cells = 0.0,
  Pericytes_Fibroblasts = 0.0
)

expr_mat <- LayerData(
  seurat_obj,
  assay = "SCT",
  layer = "data"
)

# Check that all marker genes are available.
marker_genes <- unique(
  unlist(
    lapply(
      gene_list,
      function(x) x$genes
    )
  )
)

missing_genes <- setdiff(
  marker_genes,
  rownames(expr_mat)
)

if (length(missing_genes) > 0) {
  stop(
    "Missing marker genes: ",
    paste(missing_genes, collapse = ", ")
  )
}

celltypes <- names(gene_list)
n_cells <- ncol(expr_mat)

marker_abundance <- matrix(
  NA_real_,
  nrow = n_cells,
  ncol = length(celltypes),
  dimnames = list(
    colnames(expr_mat),
    celltypes
  )
)

marker_counts <- matrix(
  0L,
  nrow = n_cells,
  ncol = length(celltypes),
  dimnames = list(
    colnames(expr_mat),
    celltypes
  )
)

for (celltype in celltypes) {

  genes <- gene_list[[celltype]]$genes

  expression_subset <- expr_mat[
    genes,
    ,
    drop = FALSE
  ]

  marker_abundance[, celltype] <-
    colSums(expression_subset)

  marker_counts[, celltype] <-
    colSums(expression_subset > 0)
}

assign_celltype_from_abundance <- function(
  marker_abundance,
  marker_counts,
  cell_type_thresholds,
  gene_list,
  default_min_ratio = 0.3
) {

  celltypes <- colnames(marker_abundance)

  assignments <- character(
    nrow(marker_abundance)
  )

  names(assignments) <-
    rownames(marker_abundance)

  for (i in seq_len(nrow(marker_abundance))) {

    abundance_values <-
      marker_abundance[i, ]

    count_values <-
      marker_counts[i, ]

    if (all(is.na(abundance_values))) {
      assignments[i] <- "Unassigned"
      next
    }

    best_index <- which.max(
      abundance_values
    )

    best_type <- celltypes[best_index]
    best_abundance <- abundance_values[best_index]
    best_count <- count_values[best_index]

    threshold <-
      cell_type_thresholds[[best_type]]

    if (
      !is.null(
        gene_list[[best_type]]$min_n_genes
      )
    ) {
      minimum_genes <-
        gene_list[[best_type]]$min_n_genes
    } else {
      minimum_genes <- max(
        1,
        ceiling(
          length(
            gene_list[[best_type]]$genes
          ) * default_min_ratio
        )
      )
    }

    if (
      is.na(best_abundance) ||
      best_abundance < threshold ||
      best_count < minimum_genes
    ) {
      assignments[i] <- "Unassigned"
    } else {
      assignments[i] <- best_type
    }
  }

  factor(
    assignments,
    levels = c(
      sort(celltypes),
      "Unassigned"
    )
  )
}

celltype_assignments <- assign_celltype_from_abundance(
  marker_abundance = marker_abundance,
  marker_counts = marker_counts,
  cell_type_thresholds = cell_type_thresholds,
  gene_list = gene_list,
  default_min_ratio = default_min_genes_ratio
)

# Explicitly align assignments to Seurat cell order.
seurat_obj$celltype_markers <-
  celltype_assignments[
    colnames(seurat_obj)
  ]

print(
  table(
    seurat_obj$celltype_markers,
    useNA = "ifany"
  )
)

# Add condition metadata
seurat_obj$eae_condition <- factor(
  paste0(
    seurat_obj$eae,
    "_",
    seurat_obj$condition
  )
)

# Keep stable clusters as identities.
Idents(seurat_obj) <- "stable_20_clusters"


output_file <- paste0(
  "/iss-scratch/CoreBioinformatics/rk720/",
  "nichenet_stuff/sucnr1_nichenet_object.rds"
)

saveRDS(
  seurat_obj,
  output_file
)

p_clusters <- DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "stable_20_clusters",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.2
) +
  ggtitle("Stable 20 clusters") +
  theme_classic()

ggsave(
  filename = paste0(
    "/iss-scratch/CoreBioinformatics/rk720/",
    "nichenet_stuff/full_dataset_umap_clusters.png"
  ),
  plot = p_clusters,
  width = 9,
  height = 7,
  dpi = 300,
  bg = "white"
)

p_celltypes <- DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "celltype_markers",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.2
) +
  ggtitle("Voting-based cell-type annotation") +
  theme_classic()

ggsave(
  filename = paste0(
    "/iss-scratch/CoreBioinformatics/rk720/",
    "nichenet_stuff/full_dataset_umap_celltype_markers.png"
  ),
  plot = p_celltypes,
  width = 11,
  height = 8,
  dpi = 300,
  bg = "white"
)

seurat_obj <- readRDS(output_file)

print(
  table(
    seurat_obj$eae_condition,
    useNA = "ifany"
  )
)

print(
  table(
    seurat_obj$celltype_markers,
    useNA = "ifany"
  )
)