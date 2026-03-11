library(Seurat)
library(ggplot2)

# Validate and extract the tissue input object
.get_ti <- function(ti) {
  stopifnot(is.list(ti), !is.null(ti$expr_mat), !is.null(ti$coords))
  stopifnot(!is.null(ti$coords$df), !is.null(ti$coords$mat))
  ti
}

# Expression getter
.get_expr <- function(ti) {
  .get_ti(ti)$expr_mat
}

# Coords getters
.get_coords_df <- function(ti) {
  .get_ti(ti)$coords$df
}
.get_xy <- function(ti) {
  .get_ti(ti)$coords$mat
}

get_coords <- function(seu) {
  coords <- GetTissueCoordinates(seu, image = NULL)
  xy <- as.matrix(coords[, c("imagecol", "imagerow")])
  rownames(xy) <- rownames(coords)
  return(list(df = coords, mat = xy))
}

# ----------------- compute neighbors and activation --------------

# For each spot, count neighbors with non-zero expression for a gene
compute_neighbors_nonzero <- function(ti,
                                      gene,
                                      k = 10) {
  ti <- .get_ti(ti)
  coords_df <- .get_coords_df(ti)
  xy <- .get_xy(ti)                       # numeric matrix of spot coordinates (rows = spots)
  expr_mat <- .get_expr(ti)

  if (!gene %in% rownames(expr_mat)) {
    stop("Gene not found in expression matrix: ", gene)
  }

  # expression values for the gene; ensure alignment to xy rownames
  expr_vals <- expr_mat[gene, , drop = TRUE]
  spot_names <- rownames(xy)
  expr_vals <- expr_vals[colnames(expr_mat) %in% spot_names]        # keep only intersect
  # Reindex to spot order; missing -> NA (then treated as zero)
  aligned <- rep(NA_real_, length(spot_names))
  names(aligned) <- spot_names
  aligned[names(expr_vals)] <- expr_vals
  aligned[is.na(aligned)] <- 0

  active_spots <- names(aligned)[aligned != 0]

  # kNN on spatial coords (exclude self by using k+1 and dropping the first)
  knn_res <- get.knn(xy, k = k + 1)
  neighbor_idx <- knn_res$nn.index

  counts <- vapply(seq_len(nrow(neighbor_idx)), function(i) {
    neigh <- neighbor_idx[i, -1]                    # drop self
    neigh_spots <- rownames(xy)[neigh]
    sum(neigh_spots %in% active_spots)
  }, integer(1))

  df <- data.frame(
    spot = rownames(xy),
    neighbors_nonzero = counts,
    is_spot_itself_active = rownames(xy) %in% active_spots,
    stringsAsFactors = FALSE
  )

  # attach image coordinates
  df$imagecol <- coords_df[df$spot, "imagecol"]
  df$imagerow <- coords_df[df$spot, "imagerow"]

  rownames(df) <- df$spot
  df
}

# Activation summary + hotspots using neighbors criterion
compute_activation_percent_neighbors <- function(ti,
                                                 genes,
                                                 k,
                                                 pct_thresh) {
  ti <- .get_ti(ti)
  expr_mat <- .get_expr(ti)

  present_genes <- genes[genes %in% rownames(expr_mat)]
  if (length(present_genes) == 0) {
    return(list(summary = data.frame(), hotspots = list()))
  }

  summary_list <- vector("list", length(present_genes))
  hotspots_list <- vector("list", length(present_genes))
  names(summary_list) <- names(hotspots_list) <- present_genes

  for (gene in present_genes) {
    neighbors_current_gene <- compute_neighbors_nonzero(ti, gene, k = k)
    neighbors_current_gene$pct_active_neighbors <- neighbors_current_gene$neighbors_nonzero / k

    # mark hotspots
    neighbors_current_gene$hotspot <- neighbors_current_gene$pct_active_neighbors >= pct_thresh

    # store hotspot IDs
    hotspots_list[[gene]] <- rownames(neighbors_current_gene)[neighbors_current_gene$hotspot]

    n_spots <- nrow(neighbors_current_gene)
    n_hotspots <- sum(neighbors_current_gene$hotspot, na.rm = TRUE)
    pct_hotspots <- n_hotspots / n_spots

    summary_list[[gene]] <- data.frame(
      gene = gene,
      n_spots = n_spots,
      n_hotspots = n_hotspots,
      pct_hotspots = pct_hotspots,
      stringsAsFactors = FALSE
    )
  }

  summary_df <- do.call(rbind, summary_list)

  list(summary = summary_df, hotspots = hotspots_list)
}

# --------------------- voting/aggregation ------------------------
voting_multi_genes_plots_spatial <- function(ti,
                                             genes,
                                             tissue,
                                             vote_fraction_threshold = 0.5,
                                             vote_count_threshold = NULL,
                                             use_fraction = TRUE,
                                             aggregation = c("mean", "sum", "median")) {
  aggregation <- match.arg(aggregation)
  ti <- .get_ti(ti)
  expr_mat <- .get_expr(ti)
  coords_df <- .get_coords_df(ti)

  message("  Running voting/aggregate on ", length(genes), " genes")
  genes_present <- genes[genes %in% rownames(expr_mat)]
  if (length(genes_present) == 0) {
    stop("None of the provided genes are present in expr_mat.")
  }

  expr_sub <- expr_mat[genes_present, , drop = FALSE]

  # compute per-spot aggregates (always compute mean, sum and median so we can return them)
  avg_expr <- as.numeric(Matrix::colMeans(expr_sub))
  sum_expr <- as.numeric(Matrix::colSums(expr_sub))
  expr_sub_dense <- as.matrix(expr_sub)   # we compute dense just once

  med_expr <- as.numeric(matrixStats::colMedians(expr_sub_dense))

  # choose the aggregation to plot
  if (aggregation == "mean") {
    agg_expr <- avg_expr
  } else if (aggregation == "sum") {
    agg_expr <- sum_expr
  } else if (aggregation == "median") {
    agg_expr <- med_expr
  } else {
    stop("Unsupported aggregation: ", aggregation)
  }

  # gene counts active (non-zero genes per spot)
  gene_counts_active <- as.integer(colSums(expr_sub != 0))

  if (use_fraction) {
    threshold_needed <- ceiling(length(genes_present) * vote_fraction_threshold)
  } else {
    threshold_needed <- if (is.null(vote_count_threshold)) 1 else vote_count_threshold
  }
  vote_pass <- gene_counts_active >= threshold_needed
  vote_fraction <- gene_counts_active / length(genes_present)

  df <- data.frame(
    spot = colnames(expr_sub),
    avg_expr = avg_expr,
    sum_expr = sum_expr,
    med_expr = med_expr,
    agg_expr = as.numeric(agg_expr),
    gene_count = gene_counts_active,
    vote_pass = as.logical(vote_pass),
    vote_fraction = as.numeric(vote_fraction),
    stringsAsFactors = FALSE
  )
  rownames(df) <- df$spot

  # attach coords (preserve order)
  # coords_df must have rownames matching spot names
  df$imagecol <- coords_df[df$spot, "imagecol"]
  df$imagerow <- coords_df[df$spot, "imagerow"]

  # label for color scale
  agg_label <- switch(aggregation,
                      "mean" = "Mean expr",
                      "sum" = "Sum expr",
                      "median" = "Median expr")

  # Plot: aggregated chosen metric
  plot_agg <- ggplot(df, aes(x = imagecol, y = imagerow, color = agg_expr)) +
    geom_point(size = 2.5) +
    scale_color_viridis_c(option = "D", name = agg_label) +
    scale_y_reverse() +
    theme_void() +
    ggtitle(paste0("Aggregated (", aggregation, ") expression - ", tissue))

  # Plot: median (kept for comparison)
  plot_med <- ggplot(df, aes(x = imagecol, y = imagerow, color = med_expr)) +
    geom_point(size = 2.5) +
    scale_color_viridis_c(option = "D", name = "Median expr") +
    scale_y_reverse() +
    theme_void() +
    ggtitle(paste0("Aggregated median expression - ", tissue))

  # Plot: vote binary
  plot_vote <- ggplot(df, aes(x = imagecol, y = imagerow, color = vote_pass)) +
    geom_point(size = 2.5) +
    scale_color_manual(name = paste0("Pass (≥", threshold_needed, " genes)"),
                       values = c("TRUE" = "red", "FALSE" = "grey80")) +
    scale_y_reverse() +
    theme_void() +
    ggtitle(paste0("Voting pass (binary) - ", tissue))

  # gene count discrete plot
  df$gene_count_discrete <- factor(df$gene_count, levels = sort(unique(df$gene_count)))
  plot_gene_count <- ggplot(df, aes(x = imagecol, y = imagerow, color = gene_count_discrete)) +
    geom_point(size = 2.5) +
    scale_color_viridis_d(name = "Active gene count", option = "D") +
    scale_y_reverse() +
    theme_void() +
    ggtitle(paste0("Active gene count per spot (discrete) - ", tissue))

  list(
    plot_agg = plot_agg,
    plot_med = plot_med,
    plot_vote = plot_vote,
    plot_gene_count = plot_gene_count,
    summary_df = df,
    aggregation = aggregation
  )
}


# -------------------------- plotting helpers ----------------------------------
plot_tissue_niches <- function(ti, feature = "niches_detailed", color_scheme_option = "D") {
  coords <- ti$coords$df
  metadata <- ti$metadata

  if (!feature %in% colnames(metadata)) {
    stop(paste("Feature", feature, "not found in metadata"))
  }

  coords$feature <- metadata[[feature]]

  if (is.factor(coords$feature) || is.character(coords$feature)) {
    ggplot(coords, aes(x = imagecol, y = imagerow, color = feature)) +
      geom_point(size = 2.5) +
      scale_color_viridis_d(name = feature, option = color_scheme_option) +
      scale_y_reverse() +
      theme_void() +
      ggtitle(feature)
  } else {
    ggplot(coords, aes(x = imagecol, y = imagerow, color = feature)) +
      geom_point(size = 2.5) +
      scale_color_viridis_c(name = feature, option = color_scheme_option) +
      scale_y_reverse() +
      theme_void() +
      ggtitle(feature)
  }
}

plot_gene_with_hotspots <- function(ti, gene, hotspot_spots) {
  ti <- .get_ti(ti)
  coords_df <- .get_coords_df(ti)
  expr_mat <- .get_expr(ti)

  if (!gene %in% rownames(expr_mat)) return(invisible(NULL))

  expr_vals <- expr_mat[gene, , drop = TRUE]
  coords_df$expression <- expr_vals[rownames(coords_df)]
  coords_df$expression[is.na(coords_df$expression)] <- 0

  coords_df$hotspot <- rownames(coords_df) %in% hotspot_spots

  p <- ggplot(coords_df, aes(x = imagecol, y = imagerow)) +
    geom_point(aes(color = expression), size = 2.5) +
    geom_point(data = subset(coords_df, hotspot),
               color = "red", size = 3, shape = 1, stroke = 1) +
    scale_color_viridis_c(name = gene) +
    scale_y_reverse() +
    theme_void() +
    ggtitle(gene)

  return(p)
}

#  Quick single-gene spatial plot from tissues_inputs_app entry
get_spatial_ggplot_gene <- function(ti, gene, gene_label = gene, color_scheme_option = "D") {
  ti <- .get_ti(ti)
  coords_df <- .get_coords_df(ti)
  expr_mat <- .get_expr(ti)

  if (!gene %in% rownames(expr_mat)) stop("Gene not found: ", gene)

  expr_vals <- expr_mat[gene, , drop = TRUE]
  coords_df$expression <- expr_vals[rownames(coords_df)]
  coords_df$expression[is.na(coords_df$expression)] <- 0

  p <- ggplot(coords_df, aes(x = imagecol, y = imagerow, color = expression)) +
    geom_point(size = 2.7) +
    scale_color_viridis_c(name = gene_label, option = color_scheme_option) +
    scale_y_reverse() +
    theme_void() +
    ggtitle(paste0("GENE: ", gene))

  return(p)
}

get_spatial_ggplot_tissue <- function(seu, coords, feature, feature_name, categorical = F, color_scheme_option="D") {

  coords$feature <- seu@meta.data[[feature]]

  if (categorical) {
    ggplot(coords, aes(x = imagecol, y = imagerow, color = feature)) +
      geom_point(size = 2.5) +
      scale_color_viridis_d(name = feature_name, option=color_scheme_option) +
      scale_y_reverse() +
      theme_void()
  } else {
    ggplot(coords, aes(x = imagecol, y = imagerow, color = feature)) +
      geom_point(size = 2.5) +
      scale_color_viridis_c(name = feature_name, option=color_scheme_option) +
      scale_y_reverse() +
      theme_void()
  }
}

expression_avg_per_niche_plot <- function(ti, genes, tissue) {
  expr_mat <- ti$expr_mat
  metadata <- ti$metadata

  genes_present <- genes[genes %in% rownames(expr_mat)]
  if (length(genes_present) == 0) {
    warning("No genes found in expression matrix for this tissue")
    return(NULL)
  }

  avg_expr <- Matrix::colMeans(expr_mat[genes_present, , drop = FALSE])

  df <- data.frame(
    cell = names(avg_expr),
    avg_expr = avg_expr
  ) %>%
    left_join(
      metadata %>% tibble::rownames_to_column("cell"),
      by = "cell"
    )

  p <- ggplot(df, aes(x = niches_detailed, y = avg_expr, fill = niches_detailed)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, size = 0.5, alpha = 0.3) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      x = "Niche",
      y = "Average expression (geneset)",
      title = paste0("Average expression per niche (", length(genes_present), " genes) - Tissue ", tissue)
    )

  return(p)
}

hotspot_counts_per_niche <- function(ti, hotspots_list, tissue, pathway_name) {
  metadata <- ti$metadata

  df <- do.call(rbind, lapply(names(hotspots_list), function(g) {
    if (length(hotspots_list[[g]]) == 0) return(NULL)
    data.frame(
      spot = hotspots_list[[g]],
      gene = g,
      stringsAsFactors = FALSE
    )
  }))

  if (is.null(df)) {
    return(NULL)
  }

  df <- df %>%
    left_join(
      metadata %>% tibble::rownames_to_column("spot") %>%
        select(spot, niches_detailed),
      by = "spot"
    )

  counts <- df %>%
    group_by(niches_detailed) %>%
    summarise(n_hotspots = n(), .groups = "drop")

  p <- ggplot(counts, aes(x = niches_detailed, y = n_hotspots, fill = niches_detailed)) +
    geom_col() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      x = "Niche",
      y = "Number of hotspots",
      title = paste0("Hotspot counts per niche - ", pathway_name, " (", tissue, ")")
    )

  return(p)
}

# Robust niche-percentage plotter: coerces list-columns to atomic vector safely
plot_niche_percentages <- function(ti, tissue) {
  if (is.null(ti) || is.null(ti$metadata)) {
    warning("plot_niche_percentages: missing tissue or metadata")
    return(NULL)
  }

  metadata <- ti$metadata

  # defensive: handle the case where niches_detailed is a list-column
  if (!"niches_detailed" %in% colnames(metadata)) {
    warning("plot_niche_percentages: metadata does not contain 'niches_detailed'")
    return(NULL)
  }

  niche_col <- metadata$niches_detailed

  # If it's a list (e.g. nested vectors), coerce to a single character per row:
  if (is.list(niche_col)) {
    # collapse multi-element entries, keep first element if multiple, convert NULL/length-0 to NA
    niche_col <- vapply(niche_col, FUN.VALUE = character(1), FUN = function(x) {
      if (is.null(x)) return(NA_character_)
      if (length(x) == 0) return(NA_character_)
      # if named or length>1, take the first element and coerce to character
      as.character(x[[1]])
    })
  } else {
    # not a list — coerce factors to character, keep as-is if already character
    if (is.factor(niche_col)) niche_col <- as.character(niche_col)
    if (!is.character(niche_col)) niche_col <- as.character(niche_col)
  }

  # attach cleaned column back to a temporary df so dplyr verbs work predictably
  tmp <- tibble::tibble(niches_detailed = niche_col)

  # drop NA rows (or keep them flagged as "unknown")
  if (all(is.na(tmp$niches_detailed))) {
    # nothing useful to plot
    return(ggplot() + theme_void() + ggtitle(paste0("No niche information available for ", tissue)))
  }

  counts <- tmp %>%
    dplyr::count(niches_detailed) %>%
    dplyr::mutate(pct = n / sum(n) * 100)

  ggplot(counts, aes(x = niches_detailed, y = pct, fill = niches_detailed)) +
    geom_col() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Niche", y = "Percentage of spots", title = paste0("Niche distribution in ", tissue))
}


# --------------- top-level runner across tissues -----------------
run_activation_neighborhood_approach <- function(tissues,
                                                 pathway_genes_list,
                                                 k = 10,
                                                 pct_thresh_active_neighbors = 0.6,
                                                 min_hotspots = 3,
                                                 tissue_inputs_map = tissues_inputs_app,
                                                 output_dir = NULL,
                                                 save_pdf = FALSE) {
  activation_df_list <- list()

  for (tissue in tissues) {
    message("=== Tissue: ", tissue, " ===")
    ti <- tissue_inputs_map[[tissue]]
    if (is.null(ti)) {
      warning("Missing tissue in tissue_inputs_map: ", tissue)
      next
    }
    ti <- .get_ti(ti)

    activation_df_list[[tissue]] <- list()

    for (pathway_name in names(pathway_genes_list)) {
      genes <- pathway_genes_list[[pathway_name]]
      activation_res <- compute_activation_percent_neighbors(
        ti,
        genes,
        k = k,
        pct_thresh = pct_thresh_active_neighbors
      )

      activation_df <- activation_res$summary
      hotspots_all  <- activation_res$hotspots

      if (nrow(activation_df) == 0) {
        activation_df_list[[tissue]][[pathway_name]] <- list(
          summary  = activation_df,
          hotspots = list(),
          plots    = list()
        )
        next
      }

      filtered_activation_df <- activation_df %>% filter(n_hotspots > min_hotspots)

      plots_list <- list()
      filtered_hotspots <- list()
      if (nrow(filtered_activation_df) > 0) {
        for (g in filtered_activation_df$gene) {
          filtered_hotspots[[g]] <- hotspots_all[[g]]
          plots_list[[g]] <- plot_gene_with_hotspots(ti, g, hotspots_all[[g]])
        }
      }

      activation_df_list[[tissue]][[pathway_name]] <- list(
        summary  = filtered_activation_df,
        hotspots = filtered_hotspots,
        plots    = plots_list
      )

      if (save_pdf && !is.null(output_dir) && nrow(filtered_activation_df) > 0) {
        out_dir_path <- file.path(output_dir, tissue, pathway_name)
        dir.create(out_dir_path, recursive = TRUE, showWarnings = FALSE)
        pdf(file.path(out_dir_path,
                      paste0("active_genes_circular_neighbors_", tissue, "_", pathway_name, ".pdf")))
        for (p in plots_list) print(p)
        dev.off()
      }
    }
  }

  activation_df_list
}

plot_compare_avg_niche <- function(
    tissues_selected,
    res_all,
    pathway_name,
    tissue_inputs_map,
    aggregation = c("mean", "sum", "median"),
    restrict_to_common_niches = FALSE,
    drop_na_niche = TRUE
) {
  aggregation <- match.arg(aggregation)

  # pick the summary function for the white diamond
  summary_fun <- switch(aggregation,
                        mean = mean,
                        sum = mean,    # show mean of sums
                        median = median)

  build_one <- function(t) {
    # genes that PASSED for this tissue & pathway
    genes_passed <- character(0)
    if (!is.null(res_all) &&
        !is.null(res_all[[t]]) &&
        !is.null(res_all[[t]][[pathway_name]]) &&
        !is.null(res_all[[t]][[pathway_name]]$summary)) {
      genes_passed <- res_all[[t]][[pathway_name]]$summary$gene
    }
    if (length(genes_passed) == 0) return(NULL)  # nothing passed in this tissue

    ti <- tissue_inputs_map[[t]]
    if (is.null(ti) || is.null(ti$expr_mat) || is.null(ti$metadata)) return(NULL)

    expr_mat <- ti$expr_mat
    md <- ti$metadata

    genes_present <- intersect(genes_passed, rownames(expr_mat))
    if (length(genes_present) == 0) return(NULL)

    expr_sub <- expr_mat[genes_present, , drop = FALSE]
    if (ncol(expr_sub) == 0) return(NULL)

    agg_value <- switch(
      aggregation,
      mean   = as.numeric(Matrix::colMeans(expr_sub)),
      sum    = as.numeric(Matrix::colSums(expr_sub)),
      median = as.numeric(matrixStats::colMedians(as.matrix(expr_sub)))
    )

    df <- data.frame(
      cell = colnames(expr_sub),
      agg_value = agg_value,
      tissue = t,
      stringsAsFactors = FALSE
    )

    # Align niches without joins (keeps length == nrow(df))
    niche_raw <- md[match(df$cell, rownames(md)), "niches_detailed"]

    niche <- if (is.list(niche_raw)) {
      vapply(niche_raw, function(x) {
        if (is.null(x) || length(x) == 0) NA_character_ else as.character(x[[1]])
      }, character(1))
    } else {
      nr <- niche_raw
      if (is.factor(nr)) nr <- as.character(nr)
      if (!is.character(nr)) nr <- as.character(nr)
      nr
    }

    df$niche <- niche
    if (drop_na_niche) {
      df <- df[!is.na(df$niche), , drop = FALSE]
      if (nrow(df) == 0) return(NULL)
    }
    df
  }

  parts <- lapply(tissues_selected, build_one)
  parts <- Filter(Negate(is.null), parts)
  if (!length(parts)) {
    return(ggplot() + theme_void() + ggtitle("No data: no genes passed or missing metadata"))
  }

  df_all <- dplyr::bind_rows(parts)

  if (restrict_to_common_niches) {
    by_tissue <- split(unique(df_all$niche), df_all$tissue)
    common <- Reduce(intersect, by_tissue)
    df_all <- df_all[df_all$niche %in% common, , drop = FALSE]
    if (nrow(df_all) == 0) {
      return(ggplot() + theme_void() + ggtitle("No common niches across selected tissues"))
    }
  }

  df_all$tissue <- factor(df_all$tissue, levels = tissues_selected)

  ggplot(df_all, aes(x = niche, y = agg_value, fill = niche)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    # geom_jitter(width = 0.2, size = 0.5, alpha = 0.3) +
    stat_summary(fun = summary_fun, geom = "point", shape = 23, size = 3, fill = "white") +
    facet_wrap(~ tissue, scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      x = "Niche",
      y = paste0(tools::toTitleCase(aggregation), " expression (passed genes)"),
      title = paste0(tools::toTitleCase(aggregation), " expression per niche across tissues (passed genes)")
    )
}

plot_compare_avg_niche_by_sample <- function(
    tissues_selected,
    res_all,
    pathway_name,
    tissue_inputs_map,
    aggregation = c("mean", "sum", "median"),
    restrict_to_common_niches = FALSE,
    drop_na_niche = TRUE,
    harmonize_genes = c("per_tissue", "common", "union")
) {
  aggregation <- match.arg(aggregation)
  harmonize_genes <- match.arg(harmonize_genes)

  # summary function used for the white diamond
  summary_fun <- switch(aggregation,
                        mean = mean,
                        sum  = mean,   # show mean of sums
                        median = median)

  # gather passed genes per tissue
  genes_passed_by_tissue <- lapply(tissues_selected, function(t) {
    print(t)
    gp <- character(0)

    if (!is.null(res_all) && !is.null(res_all[[t]]) &&
        !is.null(res_all[[t]][[pathway_name]]) &&
        !is.null(res_all[[t]][[pathway_name]]$summary)) {
      gp <- res_all[[t]][[pathway_name]]$summary$gene
    }
    print(gp)
    gp

  })
  names(genes_passed_by_tissue) <- tissues_selected

  # decide gene set to use per harmonization strategy
  gene_set_common <- NULL
  gene_set_union <- NULL
  if (harmonize_genes == "common") {
    # intersection across tissues
    gene_set_common <- Reduce(intersect, genes_passed_by_tissue)
  } else if (harmonize_genes == "union") {
    gene_set_union <- Reduce(union, genes_passed_by_tissue)
  }

  # build sample-level (tissue-level) per-niche summaries
  build_sample_niche <- function(t) {
    gp <- genes_passed_by_tissue[[t]]
    ti <- tissue_inputs_map[[t]]
    if (is.null(ti) || is.null(ti$expr_mat) || is.null(ti$metadata)) return(NULL)

    expr_mat <- ti$expr_mat
    md <- ti$metadata

    # choose genes to use for this tissue according to harmonize_genes
    genes_use <- switch(harmonize_genes,
                        per_tissue = gp,
                        common = gene_set_common,
                        union = gene_set_union)

    # if union and some genes are not present in this expr_mat, ignore them
    genes_use <- genes_use[genes_use %in% rownames(expr_mat)]
    if (length(genes_use) == 0) return(NULL)

    expr_sub <- expr_mat[genes_use, , drop = FALSE]
    if (ncol(expr_sub) == 0) return(NULL)

    # compute per-cell aggregate
    agg_value_per_cell <- switch(
      aggregation,
      mean   = as.numeric(Matrix::colMeans(expr_sub)),
      sum    = as.numeric(Matrix::colSums(expr_sub)),
      median = as.numeric(matrixStats::colMedians(as.matrix(expr_sub)))
    )

    df_cells <- data.frame(
      cell = colnames(expr_sub),
      agg_value = agg_value_per_cell,
      tissue = t,
      stringsAsFactors = FALSE
    )

    # attach niche, preserving order, handling list-columns
    niche_raw <- md[match(df_cells$cell, rownames(md)), "niches_detailed"]
    niche <- if (is.list(niche_raw)) {
      vapply(niche_raw, function(x) {
        if (is.null(x) || length(x) == 0) NA_character_ else as.character(x[[1]])
      }, character(1))
    } else {
      nr <- niche_raw
      if (is.factor(nr)) nr <- as.character(nr)
      if (!is.character(nr)) nr <- as.character(nr)
      nr
    }
    df_cells$niche <- niche

    if (drop_na_niche) df_cells <- df_cells[!is.na(df_cells$niche), , drop = FALSE]
    if (nrow(df_cells) == 0) return(NULL)

    # summarise across cells to get one sample-level value per niche
    df_sample_niche <- df_cells %>%
      dplyr::group_by(niche) %>%
      dplyr::summarise(
        sample_agg = summary_fun(agg_value, na.rm = TRUE),
        n_cells = dplyr::n(),
        .groups = "drop"
      ) %>%
      dplyr::mutate(tissue = t) %>%
      dplyr::select(tissue, niche, sample_agg, n_cells)

    # rename sample_agg to agg_value to match plotting expectations
    colnames(df_sample_niche)[colnames(df_sample_niche) == "sample_agg"] <- "agg_value"
    df_sample_niche
  }

  parts <- lapply(tissues_selected, build_sample_niche)
  parts <- Filter(Negate(is.null), parts)
  if (length(parts) == 0) {
    return(ggplot() + theme_void() + ggtitle("No data: no genes passed or missing metadata"))
  }

  df_samples <- dplyr::bind_rows(parts)

  # optionally restrict to niches that appear in all tissues
  if (restrict_to_common_niches) {
    by_tissue <- split(unique(df_samples$niche), df_samples$tissue)
    common_niches <- Reduce(intersect, by_tissue)
    df_samples <- df_samples[df_samples$niche %in% common_niches, , drop = FALSE]
    if (nrow(df_samples) == 0) {
      return(ggplot() + theme_void() + ggtitle("No common niches across selected tissues"))
    }
  }

  # make niche an ordered factor by frequency to stabilise x ordering
  niche_levels <- df_samples %>%
    dplyr::count(niche) %>%
    dplyr::arrange(desc(n)) %>%
    dplyr::pull(niche)
  df_samples$niche <- factor(df_samples$niche, levels = niche_levels)
  df_samples <- df_samples %>% filter(niche != "GM")

  # plot: boxplot of sample-level values per niche, jittered points coloured by tissue
  p <- ggplot(df_samples, aes(x = niche, y = agg_value, fill = niche)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    stat_summary(fun = summary_fun, geom = "point", shape = 23, size = 3, fill = "white") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      x = "Niche",
      y = paste0(tools::toTitleCase(aggregation), " expression"),
      title = paste0(tools::toTitleCase(aggregation), " expression per niche across samples (passed genes)")
    )

  return(p)
}

plot_compare_avg_niche_by_niche <- function(
    tissues_selected,
    res_all,
    pathway_name,
    tissue_inputs_map,
    aggregation = c("mean", "sum", "median"),
    restrict_to_common_niches = FALSE,
    drop_na_niche = TRUE,
    pad_y = 1.05   # multiply max by this for small headroom; set to 1 for exact max
) {
  aggregation <- match.arg(aggregation)

  summary_fun <- switch(aggregation,
                        mean = mean,
                        sum  = mean,
                        median = median)

  build_one <- function(t) {
    genes_passed <- character(0)
    if (!is.null(res_all) && !is.null(res_all[[t]]) &&
        !is.null(res_all[[t]][[pathway_name]]) &&
        !is.null(res_all[[t]][[pathway_name]]$summary)) {
      genes_passed <- res_all[[t]][[pathway_name]]$summary$gene
    }
    if (length(genes_passed) == 0) return(NULL)

    ti <- tissue_inputs_map[[t]]
    if (is.null(ti) || is.null(ti$expr_mat) || is.null(ti$metadata)) return(NULL)

    expr_mat <- ti$expr_mat
    md <- ti$metadata

    genes_present <- intersect(genes_passed, rownames(expr_mat))
    if (length(genes_present) == 0) return(NULL)

    expr_sub <- expr_mat[genes_present, , drop = FALSE]
    if (ncol(expr_sub) == 0) return(NULL)

    agg_value <- switch(
      aggregation,
      mean   = as.numeric(Matrix::colMeans(expr_sub)),
      sum    = as.numeric(Matrix::colSums(expr_sub)),
      median = as.numeric(matrixStats::colMedians(as.matrix(expr_sub)))
    )

    df <- data.frame(
      cell = colnames(expr_sub),
      agg_value = agg_value,
      tissue = t,
      stringsAsFactors = FALSE
    )

    niche_raw <- md[match(df$cell, rownames(md)), "niches_detailed"]
    niche <- if (is.list(niche_raw)) {
      vapply(niche_raw, function(x) {
        if (is.null(x) || length(x) == 0) NA_character_ else as.character(x[[1]])
      }, character(1))
    } else {
      nr <- niche_raw
      if (is.factor(nr)) nr <- as.character(nr)
      if (!is.character(nr)) nr <- as.character(nr)
      nr
    }

    df$niche <- niche
    if (drop_na_niche) {
      df <- df[!is.na(df$niche), , drop = FALSE]
      if (nrow(df) == 0) return(NULL)
    }
    df
  }

  parts <- lapply(tissues_selected, build_one)
  parts <- Filter(Negate(is.null), parts)
  if (!length(parts)) return(ggplot() + theme_void() + ggtitle("No data: no genes passed or missing metadata"))

  df_all <- dplyr::bind_rows(parts)

  if (restrict_to_common_niches) {
    by_tissue <- split(unique(df_all$niche), df_all$tissue)
    common <- Reduce(intersect, by_tissue)
    df_all <- df_all[df_all$niche %in% common, , drop = FALSE]
    if (nrow(df_all) == 0) return(ggplot() + theme_void() + ggtitle("No common niches across selected tissues"))
  }

  # exclude GM and drop NA niches if requested
  df_all <- df_all %>% dplyr::filter(niche != "GM")
  if (drop_na_niche) df_all <- df_all %>% dplyr::filter(!is.na(niche))
  if (nrow(df_all) == 0) return(ggplot() + theme_void() + ggtitle("No data after filtering (GM excluded)"))

  # preserve overall tissue order for consistency, but note facets will only show tissues present
  df_all$tissue <- factor(df_all$tissue, levels = tissues_selected)

  # order niches by frequency for consistent facet order
  niche_levels <- df_all %>% dplyr::count(niche) %>% dplyr::arrange(desc(n)) %>% dplyr::pull(niche)
  df_all$niche <- factor(df_all$niche, levels = niche_levels)

  # compute global max for y axis, add small padding
  max_val <- max(df_all$agg_value, na.rm = TRUE)
  if (is.na(max_val) || max_val <= 0) max_val <- 1
  y_upper <- max_val * pad_y

  ggplot(df_all, aes(x = tissue, y = agg_value, fill = tissue)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.75) +
    stat_summary(fun = summary_fun, geom = "point", shape = 23, size = 3, fill = "white") +
    facet_wrap(~ niche, scales = "free_x", nrow = 1) +   # <-- each facet shows only tissues present there
    scale_y_continuous(limits = c(0, y_upper)) +         # <-- shared y axis across facets
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(face = "bold")) +
    labs(
      x = "Tissue",
      y = paste0(tools::toTitleCase(aggregation), " expression (passed genes)"),
      title = paste0(tools::toTitleCase(aggregation), " expression per tissue (passed genes)")
    )
}