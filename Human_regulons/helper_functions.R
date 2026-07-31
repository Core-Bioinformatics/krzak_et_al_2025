# helper functions
get_assay_data_safe <- function(obj, assay, layer) {
  x <- tryCatch(
    GetAssayData(obj, assay = assay, layer = layer),
    error = function(e) NULL
  )

  if (is.null(x)) {
    x <- GetAssayData(obj, assay = assay, slot = layer)
  }

  x
}

run_microglia_ma_umap <- function(
  obj,
  dataset_name,
  assay_use = "SCT",
  n_ma = 1500,
  cell_type_col = "cell_type",
  microglia_label = "microglia",
  dims_use = 1:30,
  seed_use = 42
) {
  message("\n===== ", dataset_name, ": recalculating microglia-only MA", n_ma, " UMAP =====")

  stopifnot(cell_type_col %in% colnames(obj@meta.data))

  microglia_cells <- colnames(obj)[
    obj@meta.data[[cell_type_col]] == microglia_label
  ]

  message("Microglia cells: ", length(microglia_cells))

  obj_micro <- subset(obj, cells = microglia_cells)

  DefaultAssay(obj_micro) <- assay_use

  counts_mat <- get_assay_data_safe(
    obj_micro,
    assay = assay_use,
    layer = "counts"
  )

  most_abundant_genes <- rownames(counts_mat)[
    order(Matrix::rowSums(counts_mat), decreasing = TRUE)
  ]

  ma_genes <- most_abundant_genes[seq_len(min(n_ma, length(most_abundant_genes)))]
  ma_genes <- intersect(ma_genes, rownames(obj_micro))

  message("MA genes used: ", length(ma_genes))

  obj_micro <- ScaleData(
    obj_micro,
    assay = assay_use,
    features = ma_genes,
    verbose = TRUE
  )

  pca_name <- paste0("pca_ma", n_ma)
  umap_name <- paste0("umap_ma", n_ma)

  obj_micro <- RunPCA(
    obj_micro,
    assay = assay_use,
    features = ma_genes,
    npcs = max(dims_use),
    reduction.name = pca_name,
    reduction.key = paste0("MA", n_ma, "PC_"),
    verbose = TRUE
  )

  obj_micro <- RunUMAP(
    obj_micro,
    reduction = pca_name,
    dims = dims_use,
    reduction.name = umap_name,
    reduction.key = paste0("MA", n_ma, "UMAP_"),
    n.neighbors = 30,
    min.dist = 0.3,
    metric = "cosine",
    seed.use = seed_use,
    verbose = TRUE
  )

  list(
    object = obj_micro,
    ma_genes = ma_genes,
    pca_name = pca_name,
    umap_name = umap_name
  )
}

add_sucnr1_expression_categories <- function(obj) {
  raw_counts <- get_assay_data_safe(
    obj,
    assay = "RNA",
    layer = "counts"
  )

  sct_data <- get_assay_data_safe(
    obj,
    assay = "SCT",
    layer = "data"
  )

  if (!"SUCNR1" %in% rownames(raw_counts)) {
    stop("SUCNR1 not found in RNA counts.")
  }

  if (!"SUCNR1" %in% rownames(sct_data)) {
    stop("SUCNR1 not found in SCT data.")
  }

  obj$SUCNR1_raw_count <- as.numeric(
    raw_counts["SUCNR1", colnames(obj)]
  )

  obj$SUCNR1_sct_expr <- as.numeric(
    sct_data["SUCNR1", colnames(obj)]
  )

  obj$SUCNR1_expr <- obj$SUCNR1_raw_count
  obj$SUCNR1_positive <- obj$SUCNR1_raw_count > 0

  obj
}

make_umap_df <- function(obj, reduction_name) {
  emb <- Embeddings(obj, reduction = reduction_name) %>%
    as.data.frame()

  colnames(emb)[1:2] <- c("UMAP_1", "UMAP_2")

  emb %>%
    rownames_to_column("cell") %>%
    left_join(
      obj@meta.data %>%
        rownames_to_column("cell"),
      by = "cell"
    )
}


plot_sucnr1_only_umap <- function(
  obj,
  reduction_name,
  dataset_name,
  output_file,
  sucnr1_col = "SUCNR1_sct_expr",
  legend_title = "SUCNR1"
) {
  plot_df <- make_umap_df(obj, reduction_name)

  if (!sucnr1_col %in% colnames(plot_df)) {
    stop(sucnr1_col, " not found in metadata.")
  }

  plot_df <- plot_df %>%
    arrange(.data[[sucnr1_col]])

  p <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(
      aes(color = .data[[sucnr1_col]]),
      size = 0.5
    ) +
    scale_color_gradient(
      low = "grey85",
      high = "red",
      name = legend_title
    ) +
    theme_classic() +
    labs(
      title = paste0(dataset_name, " microglia: SUCNR1 expression"),
      x = "UMAP 1",
      y = "UMAP 2"
    )

  ggsave(
    filename = output_file,
    plot = p,
    width = 7,
    height = 6,
    dpi = 300
  )

  p
}

expression_avg_per_niche_plot <- function(
  obj,
  reduction_name,
  genes,
  dataset_name = "",
  output_file = NULL,
  assay = "SCT",
  layer = "data"
) {
  expr_mat <- get_assay_data_safe(obj, assay = assay, layer = layer)

  plot_df <- make_umap_df(obj, reduction_name)

  genes_present <- unique(genes[genes %in% rownames(expr_mat)])

  if (length(genes_present) == 0) {
    warning("No genes found in expression matrix for this tissue")
    return(NULL)
  }

  avg_expr <- Matrix::colMeans(expr_mat[genes_present, colnames(obj), drop = FALSE])

  plot_df$avg_expr <- as.numeric(avg_expr[plot_df$cell])

  plot_df <- plot_df %>%
    arrange(avg_expr)

  p <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(
      aes(color = avg_expr),
      size = 0.5
    ) +
   scale_color_viridis_c(
    option = "viridis",
    name = "Average\nexpression"
    ) +
    theme_classic() +
    labs(
      title = paste0(dataset_name, " average gene-set expression"),
      subtitle = paste0(length(genes_present), " / ", length(unique(genes)), " genes found"),
      x = "UMAP 1",
      y = "UMAP 2"
    )

  if (!is.null(output_file)) {
    ggsave(
      filename = output_file,
      plot = p,
      width = 7,
      height = 6,
      dpi = 300
    )
  }

  return(p)
}

expression_active_gene_count_plot <- function(
  obj,
  reduction_name,
  genes,
  dataset_name = "",
  output_file = NULL,
  assay = "SCT",
  layer = "data"
) {
  expr_mat <- get_assay_data_safe(obj, assay = assay, layer = layer)

  plot_df <- make_umap_df(obj, reduction_name)

  genes_present <- unique(genes[genes %in% rownames(expr_mat)])

  if (length(genes_present) == 0) {
    warning("No genes found in expression matrix for this tissue")
    return(NULL)
  }

  active_counts <- Matrix::colSums(
    expr_mat[genes_present, colnames(obj), drop = FALSE] > 0
  )

  plot_df$active_gene_count <- as.numeric(active_counts[plot_df$cell])

  plot_df <- plot_df %>%
    arrange(active_gene_count)

  p <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(
      aes(color = active_gene_count),
      size = 0.5
    ) +
    scale_color_viridis_c(
      option = "viridis",
      name = "Active\ngenes",
      breaks = seq(0, length(genes_present), by = 1)
    ) +
    theme_classic() +
    labs(
      title = paste0(dataset_name, " active gene count"),
      subtitle = paste0(
        length(genes_present), " / ", length(unique(genes)),
        " genes found; active = expression > 0"
      ),
      x = "UMAP 1",
      y = "UMAP 2"
    )

  if (!is.null(output_file)) {
    ggsave(
      filename = output_file,
      plot = p,
      width = 7,
      height = 6,
      dpi = 300
    )
  }

  return(p)
}

plot_sucnr1_proportions_by_celltype <- function(
  obj,
  genes = "SUCNR1",
  gene_set_name = NULL,
  assay = "SCT",
  layer = "data",
  cell_type_col = "cell_type",
  dataset_name = "",
  output_prefix = NULL,
  positive_threshold = 0,
  positive_only_density = FALSE
) {
  requested_genes <- unique(as.character(unlist(genes, use.names = FALSE)))
  requested_multi_gene <- length(requested_genes) > 1
  
  if (length(requested_genes) < 1) {
    stop("Please provide at least one gene.")
  }
  
  if (requested_multi_gene && is.null(gene_set_name)) {
    stop("Please provide gene_set_name when more than one gene is supplied.")
  }
  
  expr_mat <- get_assay_data_safe(obj, assay = assay, layer = layer)
  
  available_genes <- requested_genes[requested_genes %in% rownames(expr_mat)]
  missing_genes <- setdiff(requested_genes, available_genes)
  
  if (length(missing_genes) > 0) {
    warning(
      "These genes were not found and will be ignored: ",
      paste(missing_genes, collapse = ", ")
    )
  }
  
  if (length(available_genes) == 0) {
    stop(
      "None of the requested genes were found in ",
      assay,
      " ",
      layer,
      ": ",
      paste(requested_genes, collapse = ", ")
    )
  }
  
  if (!cell_type_col %in% colnames(obj@meta.data)) {
    stop(cell_type_col, " not found in obj@meta.data")
  }
  
  feature_label <- if (requested_multi_gene) {
    gene_set_name
  } else {
    available_genes[1]
  }
  
  file_label <- gsub("[^A-Za-z0-9_]+", "_", feature_label)
  
  expression_values <- if (requested_multi_gene) {
    colSums(expr_mat[available_genes, colnames(obj), drop = FALSE], na.rm = TRUE)
  } else {
    as.numeric(expr_mat[available_genes[1], colnames(obj)])
  }
  
  plot_df <- tibble(
    cell = colnames(obj),
    cell_type = obj@meta.data[colnames(obj), cell_type_col, drop = TRUE],
    expression = as.numeric(expression_values)
  ) %>%
    mutate(
      cell_type = ifelse(is.na(cell_type), "Unknown", as.character(cell_type)),
      positive = expression > positive_threshold
    )
  
  total_cells_all <- nrow(plot_df)
  total_positive_all <- sum(plot_df$positive, na.rm = TRUE)
  
  if (!requested_multi_gene && total_positive_all == 0) {
    stop(
      "No ",
      available_genes[1],
      "-positive cells found using threshold > ",
      positive_threshold
    )
  }
  
  summary_df <- plot_df %>%
    group_by(cell_type) %>%
    summarise(
      n_cells_in_celltype = n(),
      n_positive_in_celltype = sum(positive, na.rm = TRUE),
      mean_expr = mean(expression, na.rm = TRUE),
      median_expr = median(expression, na.rm = TRUE),
      pct_positive_within_celltype =
        n_positive_in_celltype / n_cells_in_celltype * 100,
      .groups = "drop"
    ) %>%
    mutate(
      total_cells_all = total_cells_all,
      total_positive_all = total_positive_all,
      pct_of_all_positive_cells = ifelse(
        total_positive_all > 0,
        n_positive_in_celltype / total_positive_all * 100,
        NA_real_
      ),
      requested_genes = paste(requested_genes, collapse = ";"),
      used_genes = paste(available_genes, collapse = ";"),
      missing_genes = paste(missing_genes, collapse = ";")
    ) %>%
    arrange(desc(mean_expr))
  
  plot_df$cell_type <- factor(plot_df$cell_type, levels = summary_df$cell_type)
  
  p_mean <- ggplot(
    summary_df,
    aes(x = reorder(cell_type, mean_expr), y = mean_expr)
  ) +
    geom_col() +
    geom_text(
      aes(label = round(mean_expr, 4)),
      hjust = -0.1,
      size = 3
    ) +
    coord_flip(clip = "off") +
    theme_classic() +
    labs(
      title = if (requested_multi_gene) {
        paste0(dataset_name, " ", feature_label, " summed expression by cell type")
      } else {
        paste0(dataset_name, " ", feature_label, " mean expression by cell type")
      },
      subtitle = if (requested_multi_gene) {
        paste0("Genes present: ", length(available_genes), " / ", length(requested_genes))
      } else {
        NULL
      },
      x = "Cell type",
      y = if (requested_multi_gene) {
        paste0("Mean summed ", assay, " ", layer, " expression")
      } else {
        paste0("Mean ", assay, " ", layer, " expression")
      }
    ) +
    scale_y_continuous(
      limits = c(0, max(summary_df$mean_expr, na.rm = TRUE) * 1.25)
    )
  
  p_pct_of_all_positive <- NULL
  p_pct_within_celltype <- NULL
  p_density <- NULL
  
  if (!requested_multi_gene) {
    gene <- available_genes[1]
    
    p_pct_of_all_positive <- ggplot(
      summary_df,
      aes(
        x = reorder(cell_type, pct_of_all_positive_cells),
        y = pct_of_all_positive_cells
      )
    ) +
      geom_col() +
      geom_text(
        aes(
          label = paste0(
            n_positive_in_celltype,
            " / ",
            total_positive_all,
            " cells, ",
            round(pct_of_all_positive_cells, 1),
            "%"
          )
        ),
        hjust = -0.1,
        size = 3
      ) +
      coord_flip(clip = "off") +
      theme_classic() +
      labs(
        title = paste0(dataset_name, " distribution of ", gene, "+ cells"),
        subtitle = paste0("Denominator = all ", gene, "+ cells in the whole object"),
        x = "Cell type",
        y = paste0("% of all ", gene, "+ cells")
      ) +
      scale_y_continuous(
        limits = c(
          0,
          max(summary_df$pct_of_all_positive_cells, na.rm = TRUE) * 1.25
        )
      )
    
    p_pct_within_celltype <- ggplot(
      summary_df,
      aes(
        x = reorder(cell_type, pct_positive_within_celltype),
        y = pct_positive_within_celltype
      )
    ) +
      geom_col() +
      geom_text(
        aes(
          label = paste0(
            n_positive_in_celltype,
            " / ",
            n_cells_in_celltype,
            " cells, ",
            round(pct_positive_within_celltype, 3),
            "%"
          )
        ),
        hjust = -0.1,
        size = 3
      ) +
      coord_flip(clip = "off") +
      theme_classic() +
      labs(
        title = paste0(dataset_name, " ", gene, "+ cells within each cell type"),
        subtitle = "Denominator = cells inside that cell type",
        x = "Cell type",
        y = paste0("% ", gene, "+ within cell type")
      ) +
      scale_y_continuous(
        limits = c(
          0,
          max(summary_df$pct_positive_within_celltype, na.rm = TRUE) * 1.25
        )
      )
    
    density_df <- plot_df
    
    if (positive_only_density) {
      density_df <- density_df %>% filter(positive)
    }
    
    p_density <- ggplot(density_df, aes(x = expression)) +
      geom_density(fill = "grey70", alpha = 0.6) +
      facet_wrap(~ cell_type, scales = "free_y") +
      theme_classic() +
      labs(
        title = paste0(dataset_name, " ", gene, " expression density by cell type"),
        subtitle = ifelse(
          positive_only_density,
          paste0("Showing only cells with ", gene, " expression > ", positive_threshold),
          "Showing all cells"
        ),
        x = paste0(gene, " ", assay, " ", layer, " expression"),
        y = "Density"
      )
  }
  
  if (!is.null(output_prefix)) {
    dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)
    
    write.csv(
      summary_df,
      paste0(output_prefix, "_", file_label, "_summary_by_celltype.csv"),
      row.names = FALSE
    )
    
    ggsave(
      paste0(output_prefix, "_", file_label, "_mean_by_celltype.png"),
      p_mean,
      width = 8,
      height = 6,
      dpi = 300
    )
    
    if (!requested_multi_gene) {
      gene <- available_genes[1]
      
      ggsave(
        paste0(output_prefix, "_", file_label, "_pct_of_all_positive_cells_by_celltype.png"),
        p_pct_of_all_positive,
        width = 9,
        height = 6,
        dpi = 300
      )
      
      ggsave(
        paste0(output_prefix, "_", file_label, "_pct_positive_within_celltype.png"),
        p_pct_within_celltype,
        width = 9,
        height = 6,
        dpi = 300
      )
      
      ggsave(
        paste0(output_prefix, "_", file_label, "_density_by_celltype.png"),
        p_density,
        width = 12,
        height = 8,
        dpi = 300
      )
    }
  }
  
  return(list(
    summary = summary_df,
    mean_plot = p_mean,
    pct_of_all_positive_cells_plot = p_pct_of_all_positive,
    pct_positive_within_celltype_plot = p_pct_within_celltype,
    density_plot = p_density,
    requested_genes = requested_genes,
    used_genes = available_genes,
    missing_genes = missing_genes
  ))
}

calc_gene_set_mean_summed_by_celltype <- function(
  obj,
  genes,
  assay = "SCT",
  layer = "data",
  cell_type_col = "cell_type",
  gene_set_name = NULL,
  sort_desc = TRUE
) {
  genes <- unique(as.character(unlist(genes, use.names = FALSE)))
  
  if (length(genes) == 0) {
    stop("Please provide at least one gene.")
  }
  
  expr_mat <- get_assay_data_safe(obj, assay = assay, layer = layer)
  
  available_genes <- genes[genes %in% rownames(expr_mat)]
  missing_genes <- setdiff(genes, available_genes)
  
  if (length(missing_genes) > 0) {
    warning(
      "These genes were not found and will be ignored: ",
      paste(missing_genes, collapse = ", ")
    )
  }
  
  if (length(available_genes) == 0) {
    stop("None of the requested genes were found.")
  }
  
  if (!cell_type_col %in% colnames(obj@meta.data)) {
    stop(cell_type_col, " not found in obj@meta.data")
  }
  
  cell_type <- obj@meta.data[colnames(expr_mat), cell_type_col, drop = TRUE]
  cell_type <- as.character(cell_type)
  cell_type[is.na(cell_type)] <- "Unknown"
  
  summed_expr <- if (length(available_genes) == 1) {
    as.numeric(expr_mat[available_genes, , drop = TRUE])
  } else if (inherits(expr_mat, "sparseMatrix")) {
    as.numeric(Matrix::colSums(expr_mat[available_genes, , drop = FALSE], na.rm = TRUE))
  } else {
    as.numeric(colSums(expr_mat[available_genes, , drop = FALSE], na.rm = TRUE))
  }
  
  agg <- rowsum(
    cbind(
      summed_expression_total = summed_expr,
      n_cells = 1
    ),
    group = cell_type,
    reorder = FALSE
  )
  
  out <- data.frame(
    cell_type = rownames(agg),
    n_cells = as.integer(agg[, "n_cells"]),
    mean_summed_expression =
      agg[, "summed_expression_total"] / agg[, "n_cells"],
    stringsAsFactors = FALSE
  )
  
  if (!is.null(gene_set_name)) {
    out$gene_set_name <- gene_set_name
  }
  
  out$n_used_genes <- length(available_genes)
  out$used_genes <- paste(available_genes, collapse = ";")
  
  if (length(missing_genes) > 0) {
    out$missing_genes <- paste(missing_genes, collapse = ";")
  }
  
  if (sort_desc) {
    out <- out[order(out$mean_summed_expression, decreasing = TRUE), ]
    rownames(out) <- NULL
  }
  
  return(out)
}

get_gene_celltype_mean_matrix <- function(
  obj,
  genes,
  assay = "SCT",
  layer = "data",
  cell_type_col = "cell_type"
) {
  expr_mat <- get_assay_data_safe(obj, assay = assay, layer = layer)
  
  genes <- unique(genes)
  genes <- genes[genes %in% rownames(expr_mat)]
  
  if (length(genes) == 0) {
    stop("No genes found in expression matrix.")
  }
  
  cell_type <- obj@meta.data[colnames(expr_mat), cell_type_col, drop = TRUE]
  cell_type <- as.character(cell_type)
  cell_type[is.na(cell_type)] <- "Unknown"
  
  cell_type <- factor(cell_type)
  
  design <- Matrix::sparse.model.matrix(~ 0 + cell_type)
  colnames(design) <- sub("^cell_type", "", colnames(design))
  
  sums <- expr_mat[genes, , drop = FALSE] %*% design
  counts <- as.numeric(table(cell_type))
  
  means <- sweep(as.matrix(sums), 2, counts, "/")
  
  return(means)
}

export_pyscenic_input <- function(
  obj,
  dataset_name,
  output_root,
  reduction_name,
  assay_use = "RNA",
  min_cell_fraction = 0.01,
  overwrite = FALSE
) {
  dataset_dir <- file.path(
    output_root,
    "input",
    dataset_name
  )

  dir.create(
    dataset_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  if (!assay_use %in% Assays(obj)) {
    stop(
      "Assay '", assay_use,
      "' is not present. Available assays: ",
      paste(Assays(obj), collapse = ", ")
    )
  }

  if (!reduction_name %in% Reductions(obj)) {
    stop(
      "Reduction '", reduction_name,
      "' is not present. Available reductions: ",
      paste(Reductions(obj), collapse = ", ")
    )
  }

  counts <- get_assay_data_safe(
    obj = obj,
    assay = assay_use,
    layer = "counts"
  )

  counts <- counts[, colnames(obj), drop = FALSE]

  if (anyDuplicated(rownames(counts))) {
    stop("Gene names are not unique.")
  }

  if (anyDuplicated(colnames(counts))) {
    stop("Cell identifiers are not unique.")
  }

  n_cells <- ncol(counts)

  min_cells_use <- max(
    3L,
    ceiling(min_cell_fraction * n_cells)
  )

  min_counts_use <- max(
    3L,
    ceiling(3 * min_cell_fraction * n_cells)
  )

  cells_per_gene <- Matrix::rowSums(counts > 0)
  counts_per_gene <- Matrix::rowSums(counts)

  keep_genes <- (
    cells_per_gene >= min_cells_use &
    counts_per_gene >= min_counts_use
  )

  gene_qc <- data.frame(
    gene = rownames(counts),
    total_counts = as.numeric(counts_per_gene),
    n_expressing_cells = as.numeric(cells_per_gene),
    kept = keep_genes,
    stringsAsFactors = FALSE
  )

  counts_filtered <- counts[
    keep_genes,
    ,
    drop = FALSE
  ]

  if (nrow(counts_filtered) == 0) {
    stop("No genes remained after filtering.")
  }

  embedding <- Embeddings(
    obj,
    reduction = reduction_name
  )

  embedding <- embedding[
    colnames(counts_filtered),
    1:2,
    drop = FALSE
  ]

  embedding <- as.data.frame(embedding)
  colnames(embedding) <- c("X", "Y")

  loom_file <- file.path(
    dataset_dir,
    paste0(dataset_name, "_microglia_scenic_input.loom")
  )

  if (file.exists(loom_file)) {
    if (!overwrite) {
      stop(
        loom_file,
        " already exists. Set overwrite = TRUE to replace it."
      )
    }

    unlink(loom_file)
  }

  message("Dataset: ", dataset_name)
  message("Cells: ", ncol(counts_filtered))
  message(
    "Genes retained: ",
    nrow(counts_filtered),
    " / ",
    nrow(counts)
  )
  message("Minimum expressing cells: ", min_cells_use)
  message("Minimum total counts: ", min_counts_use)

  SCopeLoomR::build_loom(
    file.name = loom_file,
    dgem = counts_filtered,
    title = paste(dataset_name, "microglia raw counts"),
    genome = "hg38",
    default.embedding = embedding,
    default.embedding.name = reduction_name,
    loom.spec.version = 3
  )

  metadata <- obj@meta.data[
    colnames(counts_filtered),
    ,
    drop = FALSE
  ]

  metadata <- data.frame(
    cell_id = rownames(metadata),
    metadata,
    row.names = NULL,
    check.names = FALSE
  )

  write.csv(
    metadata,
    file.path(
      dataset_dir,
      paste0(dataset_name, "_microglia_metadata.csv")
    ),
    row.names = FALSE
  )

  write.csv(
    embedding,
    file.path(
      dataset_dir,
      paste0(dataset_name, "_microglia_umap.csv")
    )
  )

  write.csv(
    gene_qc,
    file.path(
      dataset_dir,
      paste0(dataset_name, "_microglia_gene_filtering.csv")
    ),
    row.names = FALSE
  )

  summary_df <- data.frame(
    dataset = dataset_name,
    assay = assay_use,
    reduction = reduction_name,
    n_cells = ncol(counts_filtered),
    n_genes_before_filtering = nrow(counts),
    n_genes_after_filtering = nrow(counts_filtered),
    min_expressing_cells = min_cells_use,
    min_total_counts = min_counts_use,
    loom_file = loom_file
  )

  write.csv(
    summary_df,
    file.path(
      dataset_dir,
      paste0(dataset_name, "_microglia_export_summary.csv")
    ),
    row.names = FALSE
  )

  message("Created: ", loom_file)

  invisible(list(
    loom_file = loom_file,
    metadata_file = file.path(
      dataset_dir,
      paste0(dataset_name, "_microglia_metadata.csv")
    ),
    gene_qc = gene_qc,
    summary = summary_df
  ))
}

plot_gene_positive_full_umap <- function(
  obj,
  gene = "SUCNR1",
  reduction_name = "umap",
  dataset_name = "",
  output_file,
  assay = "RNA",
  layer = "counts"
) {
  expression_mat <- get_assay_data_safe(
    obj = obj,
    assay = assay,
    layer = layer
  )

  if (!gene %in% rownames(expression_mat)) {
    stop(gene, " not found in ", assay, "/", layer)
  }

  umap_df <- as.data.frame(
    Embeddings(
      obj,
      reduction = reduction_name
    )
  )

  colnames(umap_df)[1:2] <- c(
    "UMAP_1",
    "UMAP_2"
  )

  umap_df$cell <- rownames(umap_df)

  umap_df$expression <- as.numeric(
    expression_mat[
      gene,
      umap_df$cell
    ]
  )

  positive_df <- umap_df[
    umap_df$expression > 0,
    ,
    drop = FALSE
  ]

  max_expression <- max(
    umap_df$expression,
    na.rm = TRUE
  )

  colour_breaks <- seq(
    0,
    ceiling(max_expression),
    by = 1
  )

  p <- ggplot(
    umap_df,
    aes(
      x = UMAP_1,
      y = UMAP_2
    )
  ) +
    # All cells shown in grey
    geom_point(
      color = "grey85",
      size = 0.18
    ) +

    
    geom_point(
      data = positive_df,
      aes(
        color = expression
      ),
      size = 1.2
    ) +

    # geom_point( data = positive_df, shape = 21, size = 0.9, color = "black", stroke=0.3, alpha=0.2) + 

    # SUCNR1-positive cells coloured by raw count
    scale_color_gradientn(
      colours = c(
        "#E0E0E0", 
        "#FDBB84", 
        "#FC8D59", 
        "#E34A33", 
        "#B30000"  
      ),
      limits = c(
        0,
        max_expression
      ),
      breaks = colour_breaks,
      name = gene,
      guide = guide_colorbar(
        direction = "horizontal",
        title.position = "left",
        title.hjust = 0.5,
        barwidth = grid::unit(
          4.5,
          "cm"
        ),
        barheight = grid::unit(
          0.35,
          "cm"
        )
      )
    ) +

    theme_classic() +

    theme(
      legend.position = "bottom",
      legend.title = element_text(
        face = "italic",
        size = 11
      ),
      legend.text = element_text(
        size = 9
      ),
      axis.title = element_text(
        size = 11
      ),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(
        face = "italic",
        hjust = 0.5
      )
    ) +

    labs(
      title = dataset_name,
      x = "UMAP1",
      y = "UMAP2"
    )

  ggsave(
    filename = output_file,
    plot = p,
    width = 8,
    height = 6,
    dpi = 300,
    bg = "white"
  )

  p
}