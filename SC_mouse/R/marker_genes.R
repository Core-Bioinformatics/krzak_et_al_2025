add_celltype_markers <- function(
  so,
  assay = "SCT",
  layer = "data",
  default_min_ratio = 0.01,
  gene_list = NULL,
  cell_type_thresholds = NULL,
  stop_if_missing = TRUE
) {
  if (is.null(gene_list)) {
    gene_list <- list(
      "Microglia_Homeostatic" = list(
        genes = c("Sparc","P2ry12","Siglech","Tmem119","Sall1","Hexb","Csf1r"),
        min_n_genes = 3
      ),
      "Microglia_DAM" = list(genes = c("C1qc","C1qa","Cst7","Ly86","Ccl12","Spp1","Lgals3")),
      "Monocyte_derived_cells" = list(genes = c("Lyz2","H2-Ab1","Cd74","Ly6c2","Clec4e","Ccr2")),
      "Macrophages" = list(genes = c("Ly6c2","Plac8","Nos2","Mgst1","Arg1")),
      "Neutrophils" = list(genes = c("Csf3r","S100a8","S100a9")),
      "T_cells" = list(genes = c("Cd3e","Bcl11b","Il2ra","Cd247","Skap1","Themis")),
      "Neurons" = list(genes = c("Gnb4","Syt1","Snap25","Rbfox1","Sema6d","Grik2")),
      "Astrocytes" = list(genes = c("Aqp4","Slc1a2","Kcnj10","Atp1a2")),
      "Oligodendrocytes" = list(genes = c("Cnp","Plp1","Mbp","Mobp")),
      "Endothelial_Cells" = list(genes = c("Cldn5","Ly6c1","Pltp","Bsg","Itm2a","Pecam1")),
      "Ependymal_Neural_Stem_Cells" = list(genes = c("Tuba1a","Foxj1","Pifo","Nnat","Dbi","Mt3","Sox2","Pax6")),
      "Pericytes_Fibroblasts" = list(genes = c("Vtn","Myl9","Rgs5","Acta2","Tagln","Mgp"))
    )
  }

  if (is.null(cell_type_thresholds)) {
    cell_type_thresholds <- c(
      Microglia_Homeostatic        = 0.0,   # Try: 2.0, 3.0, 4.0
      Microglia_DAM                = 2.0,   # Try: 1.5, 2.5, 3.0
      Monocyte_derived_cells       = 0.0,   # Try: 2.0, 3.0
      Macrophages                  = 0.0,   # Try: 1.0, 2.0
      Neutrophils                  = 0.0,   # Try: 1.0, 2.0
      T_cells                      = 0.0,   # Try: 1.0, 2.0
      Neurons                      = 0.0,   # Try: 1.0, 2.0
      Astrocytes                   = 0.0,   # Try: 1.0, 2.0
      Oligodendrocytes             = 4.5,   # Try: 2.0, 3.0, 4.0
      Endothelial_Cells            = 0.0,   # Try: 1.0, 2.0
      Ependymal_Neural_Stem_Cells  = 0.0,   # Try: 1.0, 2.0
      Pericytes_Fibroblasts        = 0.0    # Try: 1.0, 2.0
    )
  }

  expr_mat <- LayerData(so, assay = assay, layer = layer)
  sct_features <- rownames(expr_mat)

  # missing genes check
  all_missing <- lapply(names(gene_list), function(ct) setdiff(gene_list[[ct]]$genes, sct_features))
  names(all_missing) <- names(gene_list)
  all_missing <- all_missing[lengths(all_missing) > 0]

  if (length(all_missing) > 0) {
    msg <- paste(
      "Missing genes in", paste0(assay, "/", layer), "layer:\n",
      paste(vapply(names(all_missing), function(ct) {
        paste0("  ", ct, ": ", paste(all_missing[[ct]], collapse = ", "))
      }, character(1)), collapse = "\n")
    )
    if (stop_if_missing) stop(msg) else message(msg)
  }

  celltypes <- names(gene_list)
  ncells <- ncol(expr_mat)

  marker_abundance <- matrix(
    NA_real_, nrow = ncells, ncol = length(celltypes),
    dimnames = list(colnames(expr_mat), celltypes)
  )
  marker_counts <- matrix(
    0L, nrow = ncells, ncol = length(celltypes),
    dimnames = list(colnames(expr_mat), celltypes)
  )

  for (ctype in celltypes) {
    g <- intersect(gene_list[[ctype]]$genes, sct_features)
    expr_sub <- expr_mat[g, , drop = FALSE]
    marker_abundance[, ctype] <- colSums(expr_sub)
    marker_counts[, ctype] <- colSums(expr_sub > 0)
  }

  assign_celltype_from_abundance <- function(marker_abundance, marker_counts, cell_type_thresholds, gene_list, default_min_ratio) {
    celltypes <- colnames(marker_abundance)
    ncell <- nrow(marker_abundance)
    assigned <- character(ncell)
    names(assigned) <- rownames(marker_abundance)

    for (i in seq_len(ncell)) {
      abund_i <- marker_abundance[i, ]
      counts_i <- marker_counts[i, ]

      if (all(is.na(abund_i))) {
        assigned[i] <- "Unassigned"
        next
      }

      best_idx <- which.max(abund_i)
      best_val <- abund_i[best_idx]
      best_type <- celltypes[best_idx]
      best_count <- counts_i[best_idx]

      threshold <- cell_type_thresholds[best_type]

      if (!is.null(gene_list[[best_type]]$min_n_genes)) {
        min_required <- gene_list[[best_type]]$min_n_genes
      } else {
        n_genes_total <- length(gene_list[[best_type]]$genes)
        min_required <- max(1, ceiling(n_genes_total * default_min_ratio))
      }

      if (is.na(best_val) || is.na(threshold) || best_val < threshold || best_count < min_required) {
        assigned[i] <- "Unassigned"
      } else {
        assigned[i] <- best_type
      }
    }

    factor(assigned, levels = c(sort(celltypes), "Unassigned"))
  }

  so$celltype_markers <- assign_celltype_from_abundance(
    marker_abundance,
    marker_counts,
    cell_type_thresholds,
    gene_list,
    default_min_ratio
  )


  # ------------------------------------------------------------
  # --- Add celltype palette ---
  # ------------------------------------------------------------
    celltype_palette <- c(
        Astrocytes                  = "#FFFF00",    # yellow
        Microglia_DAM               = "#79CDCD",    # turquoise
        Oligodendrocytes            = "#CDB5CD",    # lavender
        Microglia_Homeostatic       = "#4682B4",    # steel blue
        Endothelial_Cells           = "#CDAA7D",    # tan
        Neutrophils                 = "#FF00FF",    # magenta
        Ependymal_Neural_Stem_Cells = "#F4C27D",    # warm sand
        Macrophages                 = "#FF302B",    # bright red
        Neurons                     = "#FFA500",    # bright orange
        Monocyte_derived_cells      = "#123524",    # phthalo green
        Pericytes_Fibroblasts       = "#99C408",    # green
        T_cells                     = "#8B4513",    # brown
        Unassigned                  = "#A0A0A0"
    )

    # sanity: ensure all levels have a color (and no extra typos)
    lvl <- levels(so$celltype_markers)
    missing_cols <- setdiff(lvl, names(celltype_palette))
    if (length(missing_cols) > 0) {
    stop("No colors defined for: ", paste(missing_cols, collapse = ", "))
    }

    # store in Seurat object for downstream plotting
    so@misc$celltype_markers_palette <- celltype_palette
  # ------------------------------------------------------------


  so$celltype_markers_score <- apply(marker_abundance, 1, max, na.rm = TRUE)

  so$n_celltypes_best_tie <- apply(marker_abundance, 1, function(x) {
    if (all(is.na(x))) return(0L)
    max_val <- max(x, na.rm = TRUE)
    sum(x == max_val)
  })

  # optionally keep matrices for debugging
  attr(so, "celltype_markers_marker_abundance") <- marker_abundance
  attr(so, "celltype_markers_marker_counts") <- marker_counts

  so
}
