# scenic_helper_functions
read_pyscenic_auc <- function(
  auc_file,
  seurat_cells
) {
  if (!file.exists(auc_file)) {
    stop("AUC file not found: ", auc_file)
  }

  auc <- read.csv(
    auc_file,
    row.names = 1,
    check.names = FALSE
  )

  auc <- as.matrix(auc)
  storage.mode(auc) <- "numeric"

  row_matches <- sum(
    rownames(auc) %in% seurat_cells
  )

  column_matches <- sum(
    colnames(auc) %in% seurat_cells
  )

  if (column_matches > row_matches) {
    message("Transposing AUC matrix so cells are rows.")
    auc <- t(auc)
  }

  if (anyDuplicated(rownames(auc))) {
    stop("Duplicated cell identifiers found in AUC matrix.")
  }

  original_regulon_names <- colnames(auc)
  colnames(auc) <- make.unique(
    original_regulon_names
  )

  common_cells <- intersect(
    seurat_cells,
    rownames(auc)
  )

  if (length(common_cells) == 0) {
    stop(
      "No matching cells between the Seurat object ",
      "and the pySCENIC AUC matrix."
    )
  }

  message(
    "Matched cells: ",
    length(common_cells),
    " / ",
    length(seurat_cells)
  )

  list(
    auc = auc,
    common_cells = common_cells,
    regulon_name_mapping = data.frame(
      original_name = original_regulon_names,
      seurat_name = colnames(auc),
      stringsAsFactors = FALSE
    )
  )
}


add_pyscenic_assay <- function(
  obj,
  auc_file,
  assay_name = "SCENIC"
) {
  auc_result <- read_pyscenic_auc(
    auc_file = auc_file,
    seurat_cells = colnames(obj)
  )

  if (
    length(auc_result$common_cells) <
      ncol(obj)
  ) {
    warning(
      ncol(obj) -
        length(auc_result$common_cells),
      " Seurat cells were absent from the AUC matrix."
    )

    obj <- subset(
      obj,
      cells = auc_result$common_cells
    )
  }

  auc <- auc_result$auc[
    colnames(obj),
    ,
    drop = FALSE
  ]

  obj[[assay_name]] <- CreateAssayObject(
    data = t(auc)
  )

  list(
    object = obj,
    auc = auc,
    regulon_name_mapping =
      auc_result$regulon_name_mapping
  )
}


validate_gene_detection_metadata <- function(
  obj,
  gene = "SUCNR1"
) {
  count_col <- paste0(
    gene,
    "_raw_count"
  )

  positive_col <- paste0(
    gene,
    "_positive"
  )

  sct_col <- paste0(
    gene,
    "_sct_expr"
  )

  required_columns <- c(
    count_col,
    positive_col
  )

  missing_columns <- setdiff(
    required_columns,
    colnames(obj@meta.data)
  )

  if (length(missing_columns) > 0) {
    stop(
      "Missing metadata columns: ",
      paste(
        missing_columns,
        collapse = ", "
      ),
      ". Rerun generate_microglia_analysis() ",
      "using the updated helper_functions.R first."
    )
  }

  if (
    !is.logical(
      obj@meta.data[[positive_col]]
    )
  ) {
    stop(
      positive_col,
      " must be a logical column."
    )
  }

  list(
    count_col = count_col,
    positive_col = positive_col,
    sct_col = if (
      sct_col %in%
        colnames(obj@meta.data)
    ) {
      sct_col
    } else {
      NULL
    }
  )
}


summarise_gene_detection <- function(
  obj,
  gene = "SUCNR1",
  count_col = paste0(
    gene,
    "_raw_count"
  ),
  positive_col = paste0(
    gene,
    "_positive"
  )
) {
  required_columns <- c(
    count_col,
    positive_col
  )

  if (
    !all(
      required_columns %in%
        colnames(obj@meta.data)
    )
  ) {
    stop(
      "Required gene-detection metadata columns are missing."
    )
  }

  data.frame(
    gene = gene,
    n_cells = ncol(obj),
    n_positive = sum(
      obj@meta.data[[positive_col]],
      na.rm = TRUE
    ),
    n_negative = sum(
      !obj@meta.data[[positive_col]],
      na.rm = TRUE
    ),
    percent_positive = mean(
      obj@meta.data[[positive_col]],
      na.rm = TRUE
    ) * 100,
    total_counts = sum(
      obj@meta.data[[count_col]],
      na.rm = TRUE
    ),
    stringsAsFactors = FALSE
  )
}


get_mean_regulon_auc <- function(
  obj,
  scenic_assay = "SCENIC",
  group_col = NULL,
  group_value = NULL
) {
  scenic_auc <- get_assay_data_safe(
    obj = obj,
    assay = scenic_assay,
    layer = "data"
  )

  cells_use <- colnames(scenic_auc)

  if (!is.null(group_col)) {
    if (
      !group_col %in%
        colnames(obj@meta.data)
    ) {
      stop(
        group_col,
        " not found in metadata."
      )
    }

    group <- obj@meta.data[
      cells_use,
      group_col,
      drop = TRUE
    ]

    cells_use <- cells_use[
      !is.na(group) &
        group == group_value
    ]
  }

  if (length(cells_use) == 0) {
    stop(
      "No cells selected for mean AUC calculation."
    )
  }

  mean_auc <- if (
    inherits(
      scenic_auc,
      "sparseMatrix"
    )
  ) {
    Matrix::rowMeans(
      scenic_auc[
        ,
        cells_use,
        drop = FALSE
      ]
    )
  } else {
    rowMeans(
      scenic_auc[
        ,
        cells_use,
        drop = FALSE
      ]
    )
  }

  out <- data.frame(
    regulon = rownames(scenic_auc),
    mean_auc = as.numeric(mean_auc),
    n_cells = length(cells_use),
    stringsAsFactors = FALSE
  )

  out <- out[
    order(
      out$mean_auc,
      decreasing = TRUE
    ),
    ,
    drop = FALSE
  ]

  rownames(out) <- NULL

  out
}


compare_regulon_activity <- function(
  obj,
  group_col = "SUCNR1_positive",
  scenic_assay = "SCENIC"
) {
  if (
    !group_col %in%
      colnames(obj@meta.data)
  ) {
    stop(
      group_col,
      " not found in metadata."
    )
  }

  scenic_auc <- get_assay_data_safe(
    obj = obj,
    assay = scenic_assay,
    layer = "data"
  )

  cell_ids <- colnames(scenic_auc)

  group <- as.logical(
    obj@meta.data[
      cell_ids,
      group_col,
      drop = TRUE
    ]
  )

  positive_cells <- cell_ids[
    !is.na(group) & group
  ]

  negative_cells <- cell_ids[
    !is.na(group) & !group
  ]

  message(
    "Positive cells: ",
    length(positive_cells)
  )

  message(
    "Negative cells: ",
    length(negative_cells)
  )

  if (length(positive_cells) < 3) {
    stop(
      "Too few positive cells for comparison."
    )
  }

  if (length(negative_cells) < 3) {
    stop(
      "Too few negative cells for comparison."
    )
  }

  results <- lapply(
    rownames(scenic_auc),
    function(regulon) {
      positive_values <- as.numeric(
        scenic_auc[
          regulon,
          positive_cells,
          drop = TRUE
        ]
      )

      negative_values <- as.numeric(
        scenic_auc[
          regulon,
          negative_cells,
          drop = TRUE
        ]
      )

      combined_values <- c(
        positive_values,
        negative_values
      )

      if (
        length(
          unique(combined_values)
        ) < 2
      ) {
        p_value <- 1
        rank_biserial <- 0
      } else {
        test <- wilcox.test(
          positive_values,
          negative_values,
          exact = FALSE
        )

        p_value <- test$p.value

        combined_ranks <- rank(
          combined_values,
          ties.method = "average"
        )

        n_positive <- length(
          positive_values
        )

        n_negative <- length(
          negative_values
        )

        positive_rank_sum <- sum(
          combined_ranks[
            seq_len(n_positive)
          ]
        )

        mann_whitney_u <- (
          positive_rank_sum -
            n_positive *
            (n_positive + 1) / 2
        )

        rank_biserial <- (
          2 * mann_whitney_u /
            (
              n_positive *
                n_negative
            )
        ) - 1
      }

      data.frame(
        regulon = regulon,
        n_positive =
          length(positive_values),
        n_negative =
          length(negative_values),
        mean_positive = mean(
          positive_values,
          na.rm = TRUE
        ),
        mean_negative = mean(
          negative_values,
          na.rm = TRUE
        ),
        median_positive = median(
          positive_values,
          na.rm = TRUE
        ),
        median_negative = median(
          negative_values,
          na.rm = TRUE
        ),
        mean_difference = (
          mean(
            positive_values,
            na.rm = TRUE
          ) -
            mean(
              negative_values,
              na.rm = TRUE
            )
        ),
        rank_biserial =
          rank_biserial,
        p_value = p_value,
        stringsAsFactors = FALSE
      )
    }
  )

  results <- do.call(
    rbind,
    results
  )

  results$adjusted_p_value <- p.adjust(
    results$p_value,
    method = "BH"
  )

  results$direction <- ifelse(
    results$rank_biserial > 0,
    "Higher in positive cells",
    ifelse(
      results$rank_biserial < 0,
      "Lower in positive cells",
      "No difference"
    )
  )

  results <- results[
    order(
      results$adjusted_p_value,
      -abs(results$rank_biserial)
    ),
    ,
    drop = FALSE
  ]

  rownames(results) <- NULL

  results
}


get_top_differential_regulons <- function(
  regulon_results,
  direction = c(
    "higher",
    "lower"
  ),
  top_n = 6,
  adjusted_p_value_cutoff = 0.05
) {
  direction <- match.arg(direction)

  if (direction == "higher") {
    keep <- (
      regulon_results$adjusted_p_value <
        adjusted_p_value_cutoff &
      regulon_results$rank_biserial > 0
    )
  } else {
    keep <- (
      regulon_results$adjusted_p_value <
        adjusted_p_value_cutoff &
      regulon_results$rank_biserial < 0
    )
  }

  out <- regulon_results[
    keep,
    ,
    drop = FALSE
  ]

  if (nrow(out) == 0) {
    return(character(0))
  }

  if (direction == "higher") {
    out <- out[
      order(
        out$adjusted_p_value,
        -out$rank_biserial
      ),
      ,
      drop = FALSE
    ]
  } else {
    out <- out[
      order(
        out$adjusted_p_value,
        out$rank_biserial
      ),
      ,
      drop = FALSE
    ]
  }

  head(
    out$regulon,
    top_n
  )
}


get_microglia_umap_name <- function(
  obj,
  reduction_name = NULL
) {
  if (!is.null(reduction_name)) {
    if (
      !reduction_name %in%
        Reductions(obj)
    ) {
      stop(
        reduction_name,
        " not found. Available reductions: ",
        paste(
          Reductions(obj),
          collapse = ", "
        )
      )
    }

    return(reduction_name)
  }

  reduction_name <- grep(
    "^umap_ma",
    Reductions(obj),
    value = TRUE
  )

  if (length(reduction_name) == 0) {
    reduction_name <- grep(
      "umap",
      Reductions(obj),
      value = TRUE
    )
  }

  if (length(reduction_name) == 0) {
    stop(
      "No UMAP reduction found."
    )
  }

  reduction_name[1]
}


plot_regulons_on_umap <- function(
  obj,
  regulons,
  reduction_name,
  dataset_name,
  output_file,
  plot_title,
  scenic_assay = "SCENIC",
  ncol = 2
) {
  regulons <- intersect(
    regulons,
    rownames(obj[[scenic_assay]])
  )

  if (length(regulons) == 0) {
    warning(
      "No regulons available for ",
      plot_title
    )

    return(NULL)
  }

  old_assay <- DefaultAssay(obj)

  DefaultAssay(obj) <-
    scenic_assay

  p <- FeaturePlot(
    obj,
    features = regulons,
    reduction = reduction_name,
    pt.size = 0.4,
    order = TRUE,
    combine = TRUE,
    ncol = ncol,
    keep.scale = "all"
  ) +
    patchwork::plot_annotation(
      title = paste0(
        dataset_name,
        ": ",
        plot_title
      )
    )

  DefaultAssay(obj) <-
    old_assay

  ggsave(
    filename = output_file,
    plot = p,
    width = 11,
    height = 4 *
      ceiling(
        length(regulons) / ncol
      ),
    dpi = 300
  )

  p
}


plot_mean_regulon_activity_umap <- function(
  obj,
  regulons,
  reduction_name,
  positive_col,
  dataset_name,
  output_file,
  plot_title = NULL,
  gene = "SUCNR1",
  scenic_assay = "SCENIC"
) {
  regulons <- intersect(
    regulons,
    rownames(obj[[scenic_assay]])
  )

  if (length(regulons) == 0) {
    warning(
      "No regulons available for mean activity UMAP."
    )

    return(NULL)
  }

  if (
    !positive_col %in%
      colnames(obj@meta.data)
  ) {
    stop(
      positive_col,
      " not found in metadata."
    )
  }

  scenic_auc <- get_assay_data_safe(
    obj = obj,
    assay = scenic_assay,
    layer = "data"
  )

  mean_activity <- if (
    inherits(
      scenic_auc,
      "sparseMatrix"
    )
  ) {
    Matrix::colMeans(
      scenic_auc[
        regulons,
        ,
        drop = FALSE
      ]
    )
  } else {
    colMeans(
      scenic_auc[
        regulons,
        ,
        drop = FALSE
      ]
    )
  }

  umap_df <- make_umap_df(
    obj = obj,
    reduction_name = reduction_name
  )

  umap_df$mean_regulon_activity <-
    as.numeric(
      mean_activity[
        umap_df$cell
      ]
    )

  umap_df$positive <- as.logical(
    umap_df[[positive_col]]
  )

  umap_df <- umap_df[
    order(
      umap_df$mean_regulon_activity
    ),
    ,
    drop = FALSE
  ]

  positive_df <- umap_df[
    !is.na(umap_df$positive) &
      umap_df$positive,
    ,
    drop = FALSE
  ]

  p <- ggplot(
    umap_df,
    aes(
      x = UMAP_1,
      y = UMAP_2
    )
  ) +
    geom_point(
      aes(
        color =
          mean_regulon_activity
      ),
      size = 0.5
    ) +
    geom_point(
      data = positive_df,
      shape = 1,
      color = "red",
      size = 1.5,
      stroke = 0.7
    ) +
    scale_color_gradient(
      low = "grey90",
      high = "steelblue4",
      name =
        "Mean regulon\nactivity"
    ) +
    theme_classic() +
    labs(
      title = if (is.null(plot_title)) {
        paste0(
          dataset_name,
          ": mean activity of ",
          gene,
          "+ regulons"
        )
      } else {
        paste0(
          dataset_name,
          ": ",
          plot_title
        )
      },
      subtitle = paste0(
        length(regulons),
        " regulons; red circles = ",
        gene,
        "+ cells"
      ),
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


analyse_scenic_dataset <- function(
  object_file,
  auc_file,
  dataset_name,
  analysis_dir,
  output_object_file,
  gene = "SUCNR1",
  scenic_assay = "SCENIC",
  reduction_name = NULL,
  top_n_mean_auc = 11,
  top_n_umap = 6,
  top_n_mean_activity = 14,
  adjusted_p_value_cutoff = 0.05
) {
  dir.create(
    analysis_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  message(
    "\n===== ",
    dataset_name,
    ": analysing pySCENIC results ====="
  )

  obj <- readRDS(object_file)

  metadata_columns <-
    validate_gene_detection_metadata(
      obj = obj,
      gene = gene
    )

  count_col <-
    metadata_columns$count_col

  positive_col <-
    metadata_columns$positive_col

  scenic_result <- add_pyscenic_assay(
    obj = obj,
    auc_file = auc_file,
    assay_name = scenic_assay
  )

  obj <- scenic_result$object

  gene_summary <- summarise_gene_detection(
    obj = obj,
    gene = gene,
    count_col = count_col,
    positive_col = positive_col
  )

  write.csv(
    gene_summary,
    file.path(
      analysis_dir,
      paste0(
        gene,
        "_summary.csv"
      )
    ),
    row.names = FALSE
  )

  mean_auc_all <- get_mean_regulon_auc(
    obj = obj,
    scenic_assay = scenic_assay
  )

  mean_auc_positive <- get_mean_regulon_auc(
    obj = obj,
    scenic_assay = scenic_assay,
    group_col = positive_col,
    group_value = TRUE
  )

  mean_auc_negative <- get_mean_regulon_auc(
    obj = obj,
    scenic_assay = scenic_assay,
    group_col = positive_col,
    group_value = FALSE
  )

  write.csv(
    mean_auc_all,
    file.path(
      analysis_dir,
      "mean_regulon_auc_all_microglia.csv"
    ),
    row.names = FALSE
  )

  write.csv(
    mean_auc_positive,
    file.path(
      analysis_dir,
      paste0(
        "mean_regulon_auc_",
        gene,
        "_positive.csv"
      )
    ),
    row.names = FALSE
  )

  write.csv(
    mean_auc_negative,
    file.path(
      analysis_dir,
      paste0(
        "mean_regulon_auc_",
        gene,
        "_negative.csv"
      )
    ),
    row.names = FALSE
  )

  write.csv(
    head(
      mean_auc_all,
      top_n_mean_auc
    ),
    file.path(
      analysis_dir,
      paste0(
        "top_",
        top_n_mean_auc,
        "_mean_auc_all_microglia.csv"
      )
    ),
    row.names = FALSE
  )

  write.csv(
    head(
      mean_auc_positive,
      top_n_mean_auc
    ),
    file.path(
      analysis_dir,
      paste0(
        "top_",
        top_n_mean_auc,
        "_mean_auc_",
        gene,
        "_positive.csv"
      )
    ),
    row.names = FALSE
  )

  regulon_results <- compare_regulon_activity(
    obj = obj,
    group_col = positive_col,
    scenic_assay = scenic_assay
  )

  write.csv(
    regulon_results,
    file.path(
      analysis_dir,
      paste0(
        gene,
        "_positive_vs_negative_regulons.csv"
      )
    ),
    row.names = FALSE
  )

  higher_results <- regulon_results[
    regulon_results$adjusted_p_value <
      adjusted_p_value_cutoff &
      regulon_results$rank_biserial > 0,
    ,
    drop = FALSE
  ]

  lower_results <- regulon_results[
    regulon_results$adjusted_p_value <
      adjusted_p_value_cutoff &
      regulon_results$rank_biserial < 0,
    ,
    drop = FALSE
  ]

  write.csv(
    higher_results,
    file.path(
      analysis_dir,
      paste0(
        "regulons_higher_in_",
        gene,
        "_positive.csv"
      )
    ),
    row.names = FALSE
  )

  write.csv(
    lower_results,
    file.path(
      analysis_dir,
      paste0(
        "regulons_lower_in_",
        gene,
        "_positive.csv"
      )
    ),
    row.names = FALSE
  )

  all_significant_higher_regulons <-
    higher_results$regulon

  # regulon_results is ordered by:
  #   1. adjusted p-value, ascending
  #   2. absolute rank-biserial effect size, descending
  top_higher_regulon_statistics <-
    regulon_results[
      regulon_results$rank_biserial > 0,
      ,
      drop = FALSE
    ]

  top_higher_regulon_statistics <-
    head(
      top_higher_regulon_statistics,
      top_n_mean_activity
    )

  top_higher_regulon_statistics$rank <-
    seq_len(
      nrow(top_higher_regulon_statistics)
    )

  # Put rank first in the output table.
  top_higher_regulon_statistics <-
    top_higher_regulon_statistics[
      ,
      c(
        "rank",
        setdiff(
          colnames(
            top_higher_regulon_statistics
          ),
          "rank"
        )
      ),
      drop = FALSE
    ]

  top_higher_regulons_for_mean_activity <-
    top_higher_regulon_statistics$regulon

  write.csv(
    data.frame(
      regulon =
        all_significant_higher_regulons,
      stringsAsFactors = FALSE
    ),
    file.path(
      analysis_dir,
      "all_significant_higher_regulons_used_for_mean_activity_umap.csv"
    ),
    row.names = FALSE
  )

  # Names-only file retained for convenience.
  write.csv(
    data.frame(
      regulon =
        top_higher_regulons_for_mean_activity,
      stringsAsFactors = FALSE
    ),
    file.path(
      analysis_dir,
      paste0(
        "top_",
        top_n_mean_activity,
        "_higher_regulons_used_for_mean_activity_umap.csv"
      )
    ),
    row.names = FALSE
  )

  # Full statistical table for the exact regulons used
  # in the top-N mean-activity UMAP.
  write.csv(
    top_higher_regulon_statistics,
    file.path(
      analysis_dir,
      paste0(
        "top_",
        top_n_mean_activity,
        "_higher_regulons_statistics.csv"
      )
    ),
    row.names = FALSE
  )

  top_higher_regulons <-
    get_top_differential_regulons(
      regulon_results =
        regulon_results,
      direction = "higher",
      top_n = top_n_umap,
      adjusted_p_value_cutoff =
        adjusted_p_value_cutoff
    )

  top_lower_regulons <-
    get_top_differential_regulons(
      regulon_results =
        regulon_results,
      direction = "lower",
      top_n = top_n_umap,
      adjusted_p_value_cutoff =
        adjusted_p_value_cutoff
    )

  write.csv(
    data.frame(
      regulon =
        top_higher_regulons,
      stringsAsFactors = FALSE
    ),
    file.path(
      analysis_dir,
      "top_higher_regulons_used_for_umap.csv"
    ),
    row.names = FALSE
  )

  write.csv(
    data.frame(
      regulon =
        top_lower_regulons,
      stringsAsFactors = FALSE
    ),
    file.path(
      analysis_dir,
      "top_lower_regulons_used_for_umap.csv"
    ),
    row.names = FALSE
  )

  reduction_name <- get_microglia_umap_name(
    obj = obj,
    reduction_name = reduction_name
  )

  # Individual FDR-significant regulons higher
  # in SUCNR1-positive cells.
  plot_regulons_on_umap(
    obj = obj,
    regulons = top_higher_regulons,
    reduction_name = reduction_name,
    dataset_name = dataset_name,
    output_file = file.path(
      analysis_dir,
      "top_higher_regulons_umap.png"
    ),
    plot_title = paste0(
      "top regulons higher in ",
      gene,
      "+ microglia"
    ),
    scenic_assay = scenic_assay
  )

  # Individual FDR-significant regulons lower
  # in SUCNR1-positive cells.
  plot_regulons_on_umap(
    obj = obj,
    regulons = top_lower_regulons,
    reduction_name = reduction_name,
    dataset_name = dataset_name,
    output_file = file.path(
      analysis_dir,
      "top_lower_regulons_umap.png"
    ),
    plot_title = paste0(
      "top regulons lower in ",
      gene,
      "+ microglia"
    ),
    scenic_assay = scenic_assay
  )

  # Mean activity of all FDR-significant regulons
  # higher in SUCNR1-positive cells.
  plot_mean_regulon_activity_umap(
    obj = obj,
    regulons =
      all_significant_higher_regulons,
    reduction_name =
      reduction_name,
    positive_col =
      positive_col,
    dataset_name =
      dataset_name,
    output_file = file.path(
      analysis_dir,
      "all_significant_higher_regulons_mean_activity_umap.png"
    ),
    plot_title = paste0(
      "mean activity of all significant regulons higher in ",
      gene,
      "+ microglia"
    ),
    gene = gene,
    scenic_assay = scenic_assay
  )

  # Mean activity of the top positive-direction
  # regulons. These are not necessarily FDR-significant.
  plot_mean_regulon_activity_umap(
    obj = obj,
    regulons =
      top_higher_regulons_for_mean_activity,
    reduction_name =
      reduction_name,
    positive_col =
      positive_col,
    dataset_name =
      dataset_name,
    output_file = file.path(
      analysis_dir,
      paste0(
        "top_",
        top_n_mean_activity,
        "_higher_regulons_mean_activity_umap.png"
      )
    ),
    plot_title = paste0(
      "mean activity of top ",
      top_n_mean_activity,
      " regulons higher in ",
      gene,
      "+ microglia"
    ),
    gene = gene,
    scenic_assay = scenic_assay
  )

  if ("SCT" %in% Assays(obj)) {
    DefaultAssay(obj) <- "SCT"
  }

  saveRDS(
    obj,
    output_object_file
  )

  message(
    "Saved updated object: ",
    output_object_file
  )

  message(
    "Saved analysis outputs: ",
    analysis_dir
  )

  invisible(
    list(
      object = obj,
      gene_summary =
        gene_summary,
      mean_auc_all =
        mean_auc_all,
      mean_auc_positive =
        mean_auc_positive,
      mean_auc_negative =
        mean_auc_negative,
      regulon_results =
        regulon_results,
      top_higher_regulons =
        top_higher_regulons,
      top_lower_regulons =
        top_lower_regulons,
      all_significant_higher_regulons =
        all_significant_higher_regulons,
      top_higher_regulons_for_mean_activity =
        top_higher_regulons_for_mean_activity,
      top_higher_regulon_statistics =
        top_higher_regulon_statistics
    )
  )
}