library(ClustAssess)
library(Seurat)

# Callers must source 00_paths.R first (sets DATA_DIR, R_DIR, SEURAT_RDS, CLUSTASSESS_RDS).

load_annotated_object <- function(
    data_dir = DATA_DIR,
    r_dir = R_DIR,
    seurat_rds = SEURAT_RDS,
    clustassess_rds = CLUSTASSESS_RDS,
    feature_type = "Most_Abundant",
    feature_size = 1950,
    clustering_method = "SLM",
    n_clusters = 20,
    add_marker_celltypes = TRUE,
    default_min_ratio = 0.3
) {
    seurat_path <- file.path(data_dir, seurat_rds)
    automm_path <- file.path(data_dir, clustassess_rds)
    marker_script <- file.path(r_dir, "marker_genes.R")

    if (!file.exists(seurat_path)) {
        stop(
            "Missing Seurat object: ", seurat_path,
            "\nPlace the Seurat RDS in data/ as ", seurat_rds,
            " (or set SCRNA_DATA_DIR / SCRNA_SEURAT_RDS)."
        )
    }
    if (!file.exists(automm_path)) {
        stop(
            "Missing ClustAssess object: ", automm_path,
            "\nPlace the ClustAssess RDS in data/ as ", clustassess_rds,
            " (or set SCRNA_DATA_DIR / SCRNA_CLUSTASSESS_RDS)."
        )
    }
    if (add_marker_celltypes && !file.exists(marker_script)) {
        stop("Missing marker helper script: ", marker_script)
    }

    cat("Loading Seurat object from ", seurat_path, "...\n", sep = "")
    so <- readRDS(seurat_path)

    cat("Loading ClustAssess object from ", automm_path, "...\n", sep = "")
    test_automm <- readRDS(automm_path)

    chosen_clusters <- get_clusters_from_clustassess_object(
        test_automm,
        feature_type = feature_type,
        feature_size = feature_size,
        clustering_method = clustering_method,
        nclusters = n_clusters
    )

    cluster_col <- paste0("stable_", n_clusters, "_clusters")
    so[[cluster_col]] <- factor(
        chosen_clusters[[as.character(n_clusters)]]$partitions[[1]]$mb,
        levels = seq_len(n_clusters)
    )

    so$eae_condition <- ifelse(
        is.na(so$eae) | is.na(so$condition),
        NA_character_,
        paste0(so$eae, "_", so$condition)
    )
    so$eae_condition <- factor(so$eae_condition)

    if (add_marker_celltypes) {
        source(marker_script)
        so <- add_celltype_markers(
            so,
            assay = "SCT",
            layer = "data",
            default_min_ratio = default_min_ratio
        )
    }

    so
}
