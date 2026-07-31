library(ComplexHeatmap)
library(circlize)
library(grid)
library(RColorBrewer)
library(matrixStats)

panel <- data.frame(
  fcs_colname = c(
    "Sm152Di", "Nd148Di", "Er168Di", "Nd145Di", "Y89Di", "Dy164Di", "Bi209Di", "Tm169Di", 
    "Sm149Di", "Nd150Di", "Eu151Di", "Er170Di", "Yb173Di", "Tb159Di", "Sm154Di", "Lu175Di", 
    "Yb171Di", "Yb172Di",
    "Ho165Di", "Nd144Di", "Yb174Di", "Er167Di", "Pr141Di", "Gd155Di", "Yb176Di", "Dy162Di", 
    "Gd156Di", "Dy161Di"
  ),
  antigen = c(
    "CD3e", "CD11b", "CD8a", "CD4", "CD45", "CX3CR1", "CD11c", "CD206", 
    "CD19", "Ly6C", "Ly6G", "CD68", "TREM2", "P2RY12", "Dectin-1", "CCR2", 
    "CD115", "CD86", "IFNg", "IL-2", "IL-17A", "IL-6", "TNFa", "IL-4", "Foxp3", "Ki67", 
    "p-p38", "iNOS"
  ),
  marker_class = c(rep("type", 18), rep("state", 10))
)


plotDR_small <- function(..., point_size = 0.2, stroke = 0) {
  p <- plotDR(...)
  p$layers[[1]]$aes_params$size <- point_size
  p$layers[[1]]$aes_params$stroke <- stroke
  p$layers[[1]] <- ggrastr::rasterise(p$layers[[1]], dpi = 300, dev = "ragg")
  p
}

annotation_umap_point_size <- 0.4

analyze_cytof_dataset <- function(
    folder_path,
    dataset_name,
    pdf_infix_name,
    exclude_markers = character(0),
    output_gene_expression=TRUE
) {
  cat("\n======================================================\n")
  cat("Processing Dataset:", dataset_name, "\n")
  cat("Excluded markers:", paste(exclude_markers, collapse = ", "), "\n")
  cat("======================================================\n")
  
  exclude_markers <- unique(exclude_markers)
  
  panel_use <- panel
  
  if (length(exclude_markers) > 0) {
    panel_use <- panel_use[
      !tolower(panel_use$antigen) %in% tolower(exclude_markers),
      ,
      drop = FALSE
    ]
  }
  
  fcs_files <- list.files(path = folder_path, pattern = "\\.fcs$", full.names = TRUE)
  
  c_num <- sub(".*(c[0-9]+).*", "\\1", basename(fcs_files))
  
  pop <- ifelse(grepl("myeloid", fcs_files, ignore.case = TRUE), "Myeloid", "Lymphocyte")
  
  md <- data.frame(
    file_name = fcs_files,
    sample_id = paste(c_num, pop, sep = "_"),
    population = pop,
    stringsAsFactors = FALSE
  )
  
  md$condition <- ifelse(
    grepl("wt", md$file_name, ignore.case = TRUE), "WT",
    ifelse(
      grepl("ko", md$file_name, ignore.case = TRUE), "KO",
      ifelse(
        grepl("c11|c12|c13|c14", md$file_name), "WT",
        ifelse(grepl("c15|c16|c17|c18", md$file_name), "KO", "Unknown")
      )
    )
  )
  
  print(md[, c("sample_id", "condition", "population")])
  
  fs <- read.flowSet(fcs_files, transformation = FALSE, truncate_max_range = FALSE)
  md$file_name <- sampleNames(fs)
  md$patient_id <- md$sample_id 
  
  sce <- prepData(fs, panel_use, md, transform = TRUE, cofactor = 5)
  
  metadata(sce)$excluded_markers <- exclude_markers
  metadata(sce)$panel_used <- panel_use
  
  cells_per_condition <- table(sce$condition)
  min_cells <- min(cells_per_condition)
  cat("Downsampling to exactly", min_cells, "cells per WT/KO group...\n")
  
  set.seed(1234)
  wt_indices <- sample(which(sce$condition == "WT"), min_cells)
  ko_indices <- sample(which(sce$condition == "KO"), min_cells)
  sce_down <- sce[, c(wt_indices, ko_indices)]
  
  type_markers <- rownames(sce_down)[rowData(sce_down)$marker_class == "type"]
  
  cat("Calculating UMAP coordinates...\n")
  sce_down <- runDR(sce_down, dr = "UMAP", features = type_markers)
  
  dir.create(paste0("/Volumes/HD_rafael/Cytof/OUTPUTS_RAFA_v3/", pdf_infix_name), showWarnings = FALSE, recursive = TRUE)
  
  pdf(
    paste0("/Volumes/HD_rafael/Cytof/OUTPUTS_RAFA_v3/", pdf_infix_name, "/UMAP_condition_single_", pdf_infix_name, ".pdf"),
    width = 8, height = 6
  )
  print(plotDR_small(sce_down, "UMAP", color_by = "condition"))
  dev.off()
  
  pdf(
    paste0("/Volumes/HD_rafael/Cytof/OUTPUTS_RAFA_v3/", pdf_infix_name, "/UMAP_condition_facet_", pdf_infix_name, ".pdf"),
    width = 8, height = 6
  )
  print(plotDR_small(sce_down, "UMAP", color_by = "condition", facet_by = "condition"))
  dev.off()
  
  if (output_gene_expression) {
    gene_out_dir <- paste0("/Volumes/HD_rafael/Cytof/OUTPUTS_RAFA_v3/", pdf_infix_name, "/per_gene_expression")
    dir.create(gene_out_dir, showWarnings = FALSE, recursive = TRUE)
    
    for (gene in rownames(assay(sce_down, "exprs"))) {
      cat("Plotting:", gene, "\n")
      
      pdf(
        paste0(gene_out_dir, "/UMAP_", gene, "_", pdf_infix_name, ".pdf"),
        width = 8, height = 6
      )
      
      print(plotDR_small(sce_down, "UMAP", color_by = gene) + theme_classic())
      dev.off()
    }
  }
  
  return(sce_down)
}

save_gene_expression_umaps_supplementary <- function(
    sce,
    pdf_infix_name,
    output_dir = NULL,
    genes = rownames(assay(sce, "exprs")),
    dr_name = "UMAP",
    assay_name = "exprs",
    point_size = 0.45,
    alpha_val = 0.85,
    width = 5,
    height = 4,
    seed = 1234
) {
  if (is.null(output_dir)) {
    output_dir <- paste0(
      "/Volumes/HD_rafael/Cytof/OUTPUTS_RAFA_v3/",
      pdf_infix_name,
      "/per_gene_expression_single"
    )
  }
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  umap <- as.data.frame(reducedDim(sce, dr_name))
  colnames(umap)[1:2] <- c("UMAP_1", "UMAP_2")
  
  expr_mat <- assay(sce, assay_name)
  genes <- genes[genes %in% rownames(expr_mat)]
  
  for (gene in genes) {
    cat("Plotting:", gene, "\n")
    
    df <- umap
    df$expression_raw <- expr_mat[gene, ]
    
    # Scale each marker from 0 to 1, like the example-style marker panels
    rng <- range(df$expression_raw, na.rm = TRUE)
    df$expression_scaled <- if (diff(rng) == 0) {
      0
    } else {
      (df$expression_raw - rng[1]) / diff(rng)
    }
    
    # Randomise, then draw high-expression cells last so signal is visible
    set.seed(seed)
    df <- df[sample(seq_len(nrow(df))), , drop = FALSE]
    df <- df[order(df$expression_scaled), , drop = FALSE]
    
    p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, colour = expression_scaled)) +
      geom_point(
        size = point_size,
        alpha = alpha_val,
        shape = 16,
        stroke = 0
      ) +
      scale_colour_gradientn(
        colours = c("black", "#3b3b00", "#8c8c00", "#d9c900", "#ffb000", "#fff176"),
        limits = c(0, 1),
        name = "scaled\nexpression"
      ) +
      coord_equal() +
      theme_void() +
      ggtitle(gene) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14),
        legend.position = "none"
      )
    
    pdf(
      file = paste0(output_dir, "/UMAP_", gene, "_", pdf_infix_name, ".pdf"),
      width = width,
      height = height
    )
    print(p)
    dev.off()
  }
}

plotExprHeatmap_log <- function(x, features = NULL, by = c("sample_id", "cluster_id", "both"), 
                                k = "meta20", m = NULL, assay = "exprs", fun = c("median", "mean", "sum"), 
                                scale = c("first", "last", "never"), q = 0.01, 
                                row_anno = TRUE, col_anno = TRUE, row_clust = TRUE, col_clust = TRUE, 
                                row_dend = TRUE, col_dend = TRUE, bars = FALSE, perc = FALSE, 
                                bin_anno = FALSE, hm_pal = rev(brewer.pal(11, "RdYlBu")), 
                                k_pal = CATALYST:::.cluster_cols, m_pal = k_pal, distance = c("euclidean", 
                                                                                              "maximum", "manhattan", "canberra", "binary", "minkowski"), 
                                linkage = c("average", "ward.D", "single", "complete", "mcquitty", 
                                            "median", "centroid", "ward.D2")) 
{
  args <- as.list(environment())
  CATALYST:::.check_args_plotExprHeatmap(args)
  distance <- match.arg(distance)
  linkage <- match.arg(linkage)
  scale <- match.arg(scale)
  fun <- match.arg(fun)
  by <- match.arg(by)
  
  x <- x[unique(CATALYST:::.get_features(x, features)), ]
  
  if (by != "sample_id") {
    CATALYST:::.check_k(x, k)
    x$cluster_id <- cluster_ids(x, k)
  }
  if (by == "both") 
    by <- c("cluster_id", "sample_id")
  
  .do_agg <- function() {
    z <- CATALYST:::.agg(x, by, fun, assay)
    if (length(by) > 1) {
      z <- do.call("rbind", z)
      rownames(z) <- levels(x$cluster_id)
    }
    return(z)
  }
  
  .do_scale <- function() {
    if (scale == "first") {
      z <- assay(x, assay)
      z <- CATALYST:::.scale_exprs(z, 1, q)
      assay(x, assay, FALSE) <- z
      return(x)
    }
    else CATALYST:::.scale_exprs(z, 1, q)
  }
  
  z <- switch(scale, first = {
    x <- .do_scale()
    .do_agg()
  }, last = {
    z <- .do_agg()
    .do_scale()
  }, never = {
    .do_agg()
  })
  
  if (length(by) == 1) 
    z <- t(z)
  
  if (scale != "never" && !(assay == "counts" && fun == "sum")) {
    qs <- round(quantile(z, c(0.01, 0.99)) * 5)/5
    lgd_aes <- list(at = seq(qs[1], qs[2], 0.2))
  } else {
    lgd_aes <- list()
  }
  
  lgd_aes$title_gp <- gpar(fontsize = 10, fontface = "bold", lineheight = 0.8)
  sids <- levels(droplevels(factor(x$sample_id)))
  
  if (!isFALSE(row_anno)) {
    left_anno <- switch(by[1], sample_id = CATALYST:::.anno_factors(x, 
                                                                    sids, row_anno, "row"), CATALYST:::.anno_clusters(x, k, m, 
                                                                                                                      k_pal, m_pal))
  } else {
    left_anno <- NULL
  }
  
  if (!isFALSE(col_anno) && length(by) == 2) {
    top_anno <- CATALYST:::.anno_factors(x, sids, col_anno, "colum")
  } else {
    top_anno <- NULL
  }
  
  # ==========================================
  # CUSTOM LOG-SCALE BARS IMPLEMENTATION START
  # ==========================================
  if (bars) {
    # Calculate cell counts for each cluster
    nk <- table(x[[by[1]]])
    counts <- as.numeric(nk)
    
    # Calculate percentage labels
    p <- counts / sum(counts)
    txt <- sprintf("%.2f%%", p * 100)
    if (perc) {
      # Adds the cluster number (e.g. "6.7%(1)")
      txt <- paste0(txt, "(", names(nk), ")")
    }
    
    # Define log-transformed bar annotation using log10(counts + 1)
    right_anno <- rowAnnotation(
      "n_cells" = anno_barplot(
        log10(counts + 1), 
        gp = gpar(fill = "grey"),
        axis_param = list(
          side = "top",
          at = c(0, 1, 2, 3, 4, 5), # log10 values for ticks
          labels = c("1", "10", "100", "1k", "10k", "100k") # visual labels
        ),
        width = unit(2, "cm")
      ),
      "perc_text" = anno_text(txt, gp = gpar(fontsize = 8), which = "row"),
      show_annotation_name = FALSE
    )
  } else {
    right_anno <- NULL
  }
  # ==========================================
  # CUSTOM LOG-SCALE BARS IMPLEMENTATION END
  # ==========================================
  
  if (bin_anno) {
    cell_fun <- function(j, i, x, y, ...) grid.text(gp = gpar(fontsize = 8), 
                                                    sprintf("%.2f", z[i, j]), x, y)
  } else {
    cell_fun <- NULL
  }
  
  a <- ifelse(assay == "exprs", "expression", assay)
  f <- switch(fun, median = "med", fun)
  hm_title <- switch(scale, first = sprintf("%s %s\n%s", fun, 
                                            "scaled", a), last = sprintf("%s %s\n%s", "scaled", 
                                                                         fun, a), never = paste(fun, a, sep = "\n"))
  
  if (length(by) == 2) {
    col_title <- features
  } else if (length(features) == 1 && features %in% c("type", "state")) {
    col_title <- paste0(features, "_markers")
  } else {
    col_title <- ""
  }
  
  # Draw the Heatmap
  Heatmap(matrix = z, name = hm_title, col = colorRamp2(seq(min(z), 
                                                            max(z), l = n <- 100), colorRampPalette(hm_pal)(n)), 
          column_title = col_title, column_title_side = ifelse(length(by) == 
                                                                 2, "top", "bottom"), cell_fun = cell_fun, cluster_rows = row_clust, 
          cluster_columns = col_clust, show_row_dend = row_dend, 
          show_column_dend = col_dend, clustering_distance_rows = distance, 
          clustering_method_rows = linkage, clustering_distance_columns = distance, 
          clustering_method_columns = linkage, show_row_names = (is.null(left_anno) || 
                                                                   isTRUE(by == "sample_id")) && !perc, row_names_side = ifelse(by[1] == 
                                                                                                                                  "cluster_id" || isFALSE(row_anno) && !row_dend || 
                                                                                                                                  isFALSE(row_clust), "left", "right"), top_annotation = top_anno, 
          left_annotation = left_anno, right_annotation = right_anno, 
          rect_gp = gpar(col = "white"), heatmap_legend_param = lgd_aes)
}

export_clustering_median_matrix <- function(sce, output_path, meta_key) {
  expr_mat <- assay(sce, "exprs")
  cluster_ids <- cluster_ids(sce, meta_key)
  
  cluster_medians <- t(sapply(levels(cluster_ids), function(cl) {
    matrixStats::rowMedians(expr_mat[, cluster_ids == cl, drop = FALSE])
  }))
  colnames(cluster_medians) <- rownames(expr_mat)
  
  cluster_medians <- round(cluster_medians, 2)
  
  write.csv(cluster_medians, paste0(output_path, "/matrix_median_clustering_", meta_key, ".csv"))
}


clustering_cytof_dataset <- function(sce, pdf_infix_name, xdim=10, ydim=10, minK=10, maxK=30) {
  output_path <- paste0("/Volumes/HD_rafael/Cytof/OUTPUTS_RAFA_v3/", pdf_infix_name,"/clustering_results")
  
  meta_key_seq=seq(10, maxK+1, 2)
  
  cat("Running FlowSOM clustering...\n")
  type_markers <- rownames(sce)[rowData(sce)$marker_class == "type"]
  
  set.seed(1234)
  sce <- cluster(
    sce,
    features = type_markers,
    xdim = xdim,
    ydim = ydim,
    maxK = maxK,
    seed = 1234
  )
  
  dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
  for (meta_int in meta_key_seq) {
    meta_key <- paste0("meta", meta_int)
    print(paste0("meta_key: ", meta_key))

    p_map <- plotDR_small(
      sce,
      "UMAP",
      color_by = meta_key,
      point_size = annotation_umap_point_size
    ) +
      theme_classic() +
      ggtitle(paste0("UMAP - ", meta_key))

    if (identical(meta_key, "meta22")) {
      p_map <- p_map + scale_color_manual(values = meta22_cols, drop = FALSE)
    }
    
    pdf(paste0(output_path, "/map_", meta_key, "_", pdf_infix_name,".pdf"), width = 8, height = 6)
    print(p_map)
    dev.off()
    
    csv_file_dir <- paste0(output_path, "/csv_files")
    dir.create(csv_file_dir, showWarnings = FALSE, recursive = TRUE)
    export_clustering_median_matrix(sce, csv_file_dir, meta_key)
    
    pdf(paste0(output_path, "/heatmap_", meta_key, "_", pdf_infix_name,".pdf"), width = 8, height = 6)
    print(plotExprHeatmap_log(sce,
                              features = rownames(assay(sce, "exprs")), 
                              by = "cluster_id",     
                              k = meta_key,          
                              scale = "last",
                              row_clust = FALSE,
                              row_anno = FALSE,
                              bars=TRUE,
                              perc=TRUE))
    dev.off()
  }
  
  delta_plot <- metadata(sce)$delta_area
  
  zoomed_delta_plot <- delta_plot + 
    coord_cartesian(ylim = c(0, 0.07)) + 
    scale_y_continuous(breaks = seq(0, 0.07, by = 0.01)) +
    scale_x_continuous(breaks = seq(2, 50, by = 2)) +
    ggtitle(paste0("Consensus Clustering Stability (Zoomed Y-Axis) - ", pdf_infix_name)) +
    theme_classic() 
  
  pdf(paste0(output_path, "/delta_area_stability_", pdf_infix_name, ".pdf"), 
      width = 6, height = 5)
  print(zoomed_delta_plot)
  dev.off()
  
  
  
  return(sce)
}

save_plot_by_condition <- function(sce, meta_key, facet_col, pdf_infix_name) {
  
  if (is.null(meta_key) && !is.null(panel$antigen)) {
    output_path <- paste0("/Volumes/HD_rafael/Cytof/OUTPUTS_RAFA_v3/", pdf_infix_name,"/per_gene_expression")
    dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
    for (gene in rownames(assay(sce, "exprs"))) {
      pdf(paste0(output_path, "/UMAP_", gene, "_by_condition_", pdf_infix_name,".pdf"), width = 8, height = 6)
      print(plotDR_small(sce, "UMAP", color_by = gene, facet_by = facet_col) + theme_classic())
      dev.off()
    }
  } else {
    output_path <- paste0("/Volumes/HD_rafael/Cytof/OUTPUTS_RAFA_v3/", pdf_infix_name,"/clustering_results")

    p <- plotDR_small(
      sce,
      "UMAP",
      color_by = meta_key,
      facet_by = facet_col,
      point_size = annotation_umap_point_size
    ) +
      theme_classic() +
      ggtitle(paste0("UMAP - ", meta_key))

    if (identical(meta_key, "meta22")) {
      p <- p + scale_color_manual(values = meta22_cols, drop = FALSE)
    }

    pdf(paste0(output_path, "/map_", meta_key, "_by_condition_", pdf_infix_name,".pdf"), width = 8, height = 6)
    print(p)
    dev.off()
  }
}

save_cluster_umap_by_condition <- function(
    sce,
    pdf_infix_name,
    cluster_col = "meta22_cluster",
    condition_col = "condition",
    cluster_cols = meta22_cols,
    output_dir = "/Volumes/HD_rafael/Cytof/OUTPUTS_RAFA_v3"
) {
  colData(sce)[[condition_col]] <- factor(
    as.character(colData(sce)[[condition_col]]),
    levels = c("WT", "KO")
  )

  output_path <- file.path(
    output_dir,
    pdf_infix_name,
    "clusters_by_condition"
  )
  dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
  
  p <- plotDR_small(
    sce,
    "UMAP",
    color_by = cluster_col,
    facet_by = condition_col,
    point_size = annotation_umap_point_size
  ) +
    theme_classic() +
    scale_color_manual(values = cluster_cols, drop = FALSE) +
    ggtitle(paste0("UMAP - ", cluster_col, " split by ", condition_col))
  
  pdf(
    paste0(
      output_path,
      "/umap_",
      cluster_col,
      "_split_by_",
      condition_col,
      "_",
      pdf_infix_name,
      ".pdf"
    ),
    width = 14,
    height = 7
  )
  print(p)
  dev.off()
  
  return(p)
}

calculate_markers_distributions <- function(sce, pdf_infix_name) {
  pdf(paste0("/Volumes/HD_rafael/Cytof/OUTPUTS_RAFA_v3/", pdf_infix_name,"/markers_distribution_", pdf_infix_name,".pdf"), width = 8, height = 6)
  
  for (gene in rownames(assay(sce, "exprs"))) {
    cat("Plotting distribution of:", gene, "\n")
    df <- data.frame(
      expression = sce@assays@data$exprs[gene, ],
      condition = sce$condition
    )
    
    p1 <- ggplot(df, aes(x = expression)) + 
      geom_density() + 
      theme_classic() + 
      ggtitle(paste0(gene, " - ", pdf_infix_name, " (Overall)"))
    
    print(p1)
    
    p2 <- ggplot(df, aes(x = expression, fill = condition)) + 
      geom_density(alpha = 0.5) + 
      facet_wrap(~ condition) + 
      theme_classic() + 
      ggtitle(paste0(gene, " - ", pdf_infix_name, " (By Condition)")) +
      theme(legend.position = "none")
    
    print(p2)
  }
  
  dev.off()
}

#### Annotations Cory ####
global_annotations <- c(
  "1"  = "CD4 T cells",
  "2"  = "CD4 T cells",
  "3"  = "CD4 T cells",
  "4"  = "Monocytes",
  "5"  = "Transitional microglia",
  "6"  = "Monocytes",
  "7"  = "Monocytes",
  "8"  = "Homeostatic microglia",
  "9"  = "DCs",
  "10" = "CD8 T cells",
  "11" = "DAM",
  "12" = "DCs",
  "13" = "DAM",
  "14" = "CD8 T cells",
  "15" = "Monocytes",
  "16" = "Macrophages",
  "17" = "Macrophages",
  "18" = "Monocytes",
  "19" = "Macrophages",
  "20" = "Monocyte-derived infiltrating cells",
  "21" = "Neutrophils",
  "22" = "Neutrophils"
)

acute_annotations <- c(
  "1"  = "DCs",
  "2"  = "DAM",
  "3"  = "Macrophages",
  "4"  = "DAM",
  "5"  = "Transitional microglia",
  "6"  = "DAM",
  "7"  = "DAM",
  "8"  = "DAM",
  "9"  = "DAM",
  "10" = "Macrophages",
  "11" = "DCs",
  "12" = "DAM",
  "13" = "Monocytes",
  "14" = "CD4 T cells",
  "15" = "Macrophages",
  "16" = "Neutrophils",
  "17" = "CD8 T cells",
  "18" = "CD4 T cells",
  "19" = "Neutrophils",
  "20" = "B cells",
  "21" = "CD4 T cells",
  "22" = "Monocytes"
)

chronic_annotations <- c(
  "1"  = "DAM",
  "2"  = "Macrophages",
  "3"  = "DCs",
  "4"  = "CD4 T cells",
  "5"  = "Unassigned",
  "6"  = "DCs",
  "7"  = "CD4 T cells",
  "8"  = "DCs",
  "9"  = "DCs",
  "10" = "Unassigned",
  "11" = "Transitional microglia",
  "12" = "DCs",
  "13" = "Unassigned",
  "14" = "DAM",
  "15" = "CD8 T cells",
  "16" = "Monocytes",
  "17" = "CD8 T cells",
  "18" = "CD8 T cells",
  "19" = "Macrophages",
  "20" = "Monocytes",
  "21" = "Neutrophils",
  "22" = "Homeostatic microglia"
)


#### Broad cell-type colours ####
all_cory_cols <- c(
  "Homeostatic microglia" = "#4A86B8", 
  "Transitional microglia" = "#3E8F8A", 
  "DAM" = "#78C9C7",                  
  "DCs" = "#C8B1CC",                   
  "Macrophages" = "#FF3B30",           
  "Monocytes" = "#F4B35E",              
  "Monocyte-derived infiltrating cells" = "#123D2A",
  "Neutrophils" = "#E600E6",            
  "CD4 T cells" = "#8B4A16",            
  "CD8 T cells" = "#B06A2B",            
  "B cells" = "#7B1E3A",                
  "Unassigned" = "grey70"
)

#### meta22 cluster colours ####
# Colours 1-18, 19 and 21 reproduce Cory's supplied palette. Clusters 20 and
# 22 use higher-contrast replacements for the original light grey/pale brown.
meta22_cols <- c(
  "1"  = "#CA2C22",
  "2"  = "#EB8777",
  "3"  = "#3163AB",
  "4"  = "#86AEDA",
  "5"  = "#7E346F",
  "6"  = "#A97DA3",
  "7"  = "#EF8632",
  "8"  = "#F2B770",
  "9"  = "#D43E88",
  "10" = "#DA8EC0",
  "11" = "#549E3F",
  "12" = "#BBDE93",
  "13" = "#689FAF",
  "14" = "#9DD1C7",
  "15" = "#9F7831",
  "16" = "#DDAD3B",
  "17" = "#7470AE",
  "18" = "#BBAFD1",
  "19" = "#666666",
  "20" = "#00A6D6",
  "21" = "#A48483",
  "22" = "#5B1A1A"
)

make_cols <- function(annotation_map, cols = all_cory_cols) {
  missing <- setdiff(unique(unname(annotation_map)), names(cols))
  
  if (length(missing) > 0) {
    stop("Missing colours for: ", paste(missing, collapse = ", "))
  }
  
  cols[names(cols) %in% unique(unname(annotation_map))]
}

global_cols  <- make_cols(global_annotations)
acute_cols   <- make_cols(acute_annotations)
chronic_cols <- make_cols(chronic_annotations)


add_bcell_umap_gate_annotation <- function(
    sce,
    pdf_infix_name,
    source_col = "meta_cory_annotation",
    output_col = "meta_cory_annotation_b_cells_included",
    gate_col = "manual_bcell_gate",
    dr_name = "UMAP",
    gate_type = c("ellipse", "rectangle"),
    
    # Ellipse gate
    x_center = NULL,
    y_center = NULL,
    x_radius = NULL,
    y_radius = NULL,
    
    # Rectangle gate
    x_min = NULL,
    x_max = NULL,
    y_min = NULL,
    y_max = NULL,
    
    bcell_label = "B cells",
    cols_map = NULL,
    output_dir = "/Volumes/HD_rafael/Cytof/OUTPUTS_RAFA_v3",
    facet_col = "condition",
    point_size = annotation_umap_point_size,
    alpha_val = 0.6
) {
  
  gate_type <- match.arg(gate_type)
  
  if (!source_col %in% colnames(colData(sce))) {
    stop("source_col not found in colData(sce): ", source_col)
  }
  
  umap <- as.data.frame(reducedDim(sce, dr_name))
  colnames(umap)[1:2] <- c("UMAP_1", "UMAP_2")
  
  if (gate_type == "ellipse") {
    
    required_vals <- c(x_center, y_center, x_radius, y_radius)
    
    if (any(is.null(required_vals))) {
      stop(
        "For ellipse gate, provide x_center, y_center, x_radius, and y_radius."
      )
    }
    
    in_gate <- (
      ((umap$UMAP_1 - x_center) / x_radius)^2 +
        ((umap$UMAP_2 - y_center) / y_radius)^2
    ) <= 1
  }
  
  if (gate_type == "rectangle") {
    
    required_vals <- c(x_min, x_max, y_min, y_max)
    
    if (any(is.null(required_vals))) {
      stop(
        "For rectangle gate, provide x_min, x_max, y_min, and y_max."
      )
    }
    
    in_gate <- (
      umap$UMAP_1 >= x_min &
        umap$UMAP_1 <= x_max &
        umap$UMAP_2 >= y_min &
        umap$UMAP_2 <= y_max
    )
  }
  
  base_annotation <- as.character(colData(sce)[[source_col]])
  new_annotation <- base_annotation
  new_annotation[in_gate] <- bcell_label
  
  base_levels <- levels(colData(sce)[[source_col]])
  
  if (is.null(base_levels)) {
    base_levels <- unique(base_annotation)
  }
  
  new_levels <- unique(c(base_levels, bcell_label))
  
  colData(sce)[[output_col]] <- factor(
    new_annotation,
    levels = new_levels
  )
  
  colData(sce)[[gate_col]] <- in_gate
  
  metadata(sce)$manual_bcell_gate <- list(
    source_col = source_col,
    output_col = output_col,
    gate_col = gate_col,
    dr_name = dr_name,
    gate_type = gate_type,
    x_center = x_center,
    y_center = y_center,
    x_radius = x_radius,
    y_radius = y_radius,
    x_min = x_min,
    x_max = x_max,
    y_min = y_min,
    y_max = y_max,
    n_cells = sum(in_gate)
  )
  
  cat("B-cell gated cells:", sum(in_gate), "\n")
  
  output_path <- file.path(
    output_dir,
    pdf_infix_name,
    "cell_types_by_condition"
  )
  dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
  
  if (!is.null(cols_map)) {
    
    cols_map2 <- cols_map
    
    if (!bcell_label %in% names(cols_map2)) {
      cols_map2 <- c(cols_map2, setNames(all_cory_cols[[bcell_label]], bcell_label))
    }
    
    cols_map2 <- cols_map2[names(cols_map2) %in% levels(colData(sce)[[output_col]])]
    
    p_ann <- plotDR_small(
      sce,
      "UMAP",
      color_by = output_col,
      point_size = point_size
    ) +
      theme_classic() +
      scale_color_manual(values = cols_map2, drop = FALSE) +
      ggtitle(paste0(output_col, " - ", pdf_infix_name))
    
    p_ann_facet <- plotDR_small(
      sce,
      "UMAP",
      color_by = output_col,
      facet_by = facet_col,
      point_size = point_size
    ) +
      theme_classic() +
      scale_color_manual(values = cols_map2, drop = FALSE) +
      ggtitle(paste0(output_col, " by ", facet_col, " - ", pdf_infix_name))
    
    pdf(
      file.path(
        output_path,
        paste0("umap_", output_col, "_", pdf_infix_name, ".pdf")
      ),
      width = 12,
      height = 10
    )
    print(p_ann)
    dev.off()
    
    pdf(
      file.path(
        output_path,
        paste0("umap_", output_col, "_by_", facet_col, "_", pdf_infix_name, ".pdf")
      ),
      width = 14,
      height = 7
    )
    print(p_ann_facet)
    dev.off()
  }
  
  gate_df <- umap
  gate_df$in_gate <- in_gate
  
  p_gate <- ggplot(gate_df, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(
      aes(colour = in_gate),
      size = point_size,
      alpha = alpha_val
    ) +
    scale_colour_manual(
      values = c("FALSE" = "grey80", "TRUE" = all_cory_cols[[bcell_label]]),
      name = "B-cell gate"
    ) +
    coord_equal() +
    theme_classic() +
    ggtitle(paste0("Manual B-cell gate - ", pdf_infix_name))
  
  if (gate_type == "ellipse") {
    theta <- seq(0, 2 * pi, length.out = 300)
    
    ellipse_df <- data.frame(
      UMAP_1 = x_center + x_radius * cos(theta),
      UMAP_2 = y_center + y_radius * sin(theta)
    )
    
    p_gate <- p_gate +
      geom_path(
        data = ellipse_df,
        aes(x = UMAP_1, y = UMAP_2),
        inherit.aes = FALSE,
        colour = "black",
        linewidth = 0.8
      )
  }
  
  if (gate_type == "rectangle") {
    p_gate <- p_gate +
      geom_rect(
        xmin = x_min,
        xmax = x_max,
        ymin = y_min,
        ymax = y_max,
        inherit.aes = FALSE,
        fill = NA,
        colour = "black",
        linewidth = 0.8
      )
  }
  
  pdf(
    file.path(
      output_path,
      paste0("umap_manual_bcell_gate_check_", pdf_infix_name, ".pdf")
    ),
    width = 8,
    height = 6
  )
  print(p_gate)
  dev.off()
  
  return(sce)
}


save_proportion_plots_by_condition <- function(
    sce,
    pdf_infix_name,
    cluster_col = "meta22_cluster",
    celltype_col = NULL,
    condition_col = "condition",
    cluster_cols = meta22_cols,
    celltype_cols = all_cory_cols,
    output_dir = "/Volumes/HD_rafael/Cytof/OUTPUTS_RAFA_v3"
) {
  if (is.null(celltype_col)) {
    preferred_celltype_col <- "meta_cory_annotation_b_cells_included"
    celltype_col <- if (preferred_celltype_col %in% colnames(colData(sce))) {
      preferred_celltype_col
    } else {
      "meta_cory_annotation"
    }
  }

  required_cols <- c(cluster_col, celltype_col, condition_col)
  missing_cols <- required_cols[!required_cols %in% colnames(colData(sce))]

  if (length(missing_cols) > 0) {
    stop("Columns not found in colData(sce): ", paste(missing_cols, collapse = ", "))
  }

  plot_data <- data.frame(
    condition = factor(
      as.character(colData(sce)[[condition_col]]),
      levels = c("WT", "KO")
    ),
    cluster = factor(
      as.character(colData(sce)[[cluster_col]]),
      levels = names(cluster_cols)
    ),
    cell_type = factor(
      as.character(colData(sce)[[celltype_col]]),
      levels = names(celltype_cols)
    )
  ) |>
    dplyr::filter(!is.na(condition))

  plot_data$cell_type <- droplevels(plot_data$cell_type)

  cluster_proportions <- plot_data |>
    dplyr::filter(!is.na(cluster)) |>
    dplyr::count(condition, cluster, name = "n_cells") |>
    tidyr::complete(condition, cluster, fill = list(n_cells = 0)) |>
    dplyr::group_by(condition) |>
    dplyr::mutate(proportion = n_cells / sum(n_cells)) |>
    dplyr::ungroup()

  celltype_proportions <- plot_data |>
    dplyr::filter(!is.na(cell_type)) |>
    dplyr::count(condition, cell_type, name = "n_cells") |>
    tidyr::complete(condition, cell_type, fill = list(n_cells = 0)) |>
    dplyr::group_by(condition) |>
    dplyr::mutate(proportion = n_cells / sum(n_cells)) |>
    dplyr::ungroup()

  shared_theme <- theme_classic() +
    theme(
      axis.title.x = element_blank(),
      legend.title = element_blank()
    )

  p_clusters <- ggplot(
    cluster_proportions,
    aes(x = condition, y = proportion, fill = cluster)
  ) +
    geom_col(width = 0.8, colour = "white", linewidth = 0.15) +
    scale_fill_manual(values = cluster_cols, drop = FALSE) +
    scale_y_continuous(
      labels = scales::label_percent(accuracy = 1),
      expand = expansion(mult = c(0, 0))
    ) +
    coord_cartesian(ylim = c(0, 1), expand = FALSE) +
    labs(
      title = "Cluster composition by condition",
      y = "Proportion of cells"
    ) +
    shared_theme

  p_celltypes <- ggplot(
    celltype_proportions,
    aes(x = condition, y = proportion, fill = cell_type)
  ) +
    geom_col(width = 0.8, colour = "white", linewidth = 0.15) +
    scale_fill_manual(values = celltype_cols, drop = TRUE) +
    scale_y_continuous(
      labels = scales::label_percent(accuracy = 1),
      expand = expansion(mult = c(0, 0))
    ) +
    coord_cartesian(ylim = c(0, 1), expand = FALSE) +
    labs(
      title = "Cell-type composition by condition",
      y = "Proportion of cells"
    ) +
    shared_theme

  output_path <- file.path(output_dir, pdf_infix_name, "proportions_by_condition")
  dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

  ggsave(
    file.path(
      output_path,
      paste0("proportion_clusters_by_condition_", pdf_infix_name, ".pdf")
    ),
    p_clusters,
    width = 6,
    height = 7
  )

  ggsave(
    file.path(
      output_path,
      paste0("proportion_cell_types_by_condition_", pdf_infix_name, ".pdf")
    ),
    p_celltypes,
    width = 8,
    height = 7
  )

  invisible(list(
    cluster_proportions = cluster_proportions,
    celltype_proportions = celltype_proportions,
    plots = list(clusters = p_clusters, cell_types = p_celltypes)
  ))
}


annotate_based_on_maps <- function(
    sce,
    meta_key,
    pdf_infix_name,
    annotation_map,
    cols_map,
    facet_col = "condition"
) {
  output_path <- paste0("/Volumes/HD_rafael/Cytof/OUTPUTS_RAFA_v3/", pdf_infix_name, "/cell_types_by_condition")
  dir.create(output_path, recursive=T)
  
  meta_ids <- cluster_ids(sce, k = meta_key)
  
  sce$meta_cory_annotation <- annotation_map[as.character(meta_ids)]
  
  sce$meta_cory_annotation <- factor(
    sce$meta_cory_annotation,
    levels = names(cols_map)
  )
  
  p <- plotDR_small(
    sce,
    "UMAP",
    color_by = "meta_cory_annotation",
    point_size = annotation_umap_point_size
  ) +
    theme_classic() +
    scale_color_manual(values = cols_map, drop = FALSE)
  
  p_facet <- plotDR_small(
    sce,
    "UMAP",
    color_by = "meta_cory_annotation",
    facet_by = facet_col,
    point_size = annotation_umap_point_size
  ) +
    theme_classic() +
    scale_color_manual(values = cols_map, drop = FALSE)
  
  pdf(
    paste0(output_path, "/umap_annotations_cory_", pdf_infix_name, ".pdf"),
    width = 12, height = 10
  )
  print(p)
  dev.off()
  
  pdf(
    paste0(output_path, "/umap_annotations_cory_by_", facet_col, "_", pdf_infix_name, ".pdf"),
    width = 14, height = 7
  )
  print(p_facet)
  dev.off()
  
  return(sce)
}

plot_marker_heatmap_by_celltype <- function(
    sce,
    pdf_infix_name,
    annotation_col = "meta_cory_annotation",
    assay_name = "exprs",
    markers = panel$antigen
) {
  output_path <- paste0("/Volumes/HD_rafael/Cytof/OUTPUTS_RAFA_v3/", pdf_infix_name, "/marker_expression_by_celltype")
  dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
  
  if (!annotation_col %in% colnames(colData(sce))) {
    stop("Column not found in colData(sce): ", annotation_col)
  }
  
  markers <- markers[markers %in% rownames(assay(sce, assay_name))]
  
  ann <- as.character(colData(sce)[[annotation_col]])
  keep <- !is.na(ann)
  
  expr_mat <- assay(sce, assay_name)[markers, keep, drop = FALSE]
  ann <- ann[keep]
  
  celltypes <- sort(unique(ann))
  
  # median_mat <- sapply(celltypes, function(ct) {
  #   matrixStats::rowMedians(expr_mat[, ann == ct, drop = FALSE], na.rm = TRUE)
  # })
  # 
  # rownames(median_mat) <- markers
  
  # write.csv(
  #   median_mat,
  #   paste0(output_path, "/median_marker_expression_by_celltype_", pdf_infix_name, ".csv")
  # )
  # 
  # # Raw median expression heatmap
  # pdf(
  #   paste0(output_path, "/heatmap_raw_marker_expression_by_celltype_", pdf_infix_name, ".pdf"),
  #   width = 14, height = 8
  # )
  # 
  # print(
  #   Heatmap(
  #     median_mat,
  #     name = "Median\nexprs",
  #     col = colorRamp2(
  #       seq(min(median_mat), max(median_mat), length.out = 100),
  #       colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100)
  #     ),
  #     cluster_rows = TRUE,
  #     cluster_columns = TRUE,
  #     rect_gp = gpar(col = "white"),
  #     column_names_rot = 45
  #   )
  # )
  # 
  # dev.off()
  
  # Row-scaled heatmap, usually easier biologically
  scaled_mat <- t(scale(t(median_mat)))
  scaled_mat[is.na(scaled_mat)] <- 0
  
  pdf(
    paste0(output_path, "/heatmap_scaled_marker_expression_by_celltype_", pdf_infix_name, ".pdf"),
    width = 14, height = 8
  )
  
  print(
    Heatmap(
      scaled_mat,
      name = "Row-scaled\nexprs",
      col = colorRamp2(
        c(-2, 0, 2),
        c("#2166AC", "white", "#B2182B")
      ),
      cluster_rows = TRUE,
      cluster_columns = TRUE,
      rect_gp = gpar(col = "white"),
      column_names_rot = 45
    )
  )
  
  dev.off()
  
  return(median_mat)
}


plot_marker_heatmap_by_condition <- function(
    sce,
    pdf_infix_name,
    annotation_col = "meta_cory_annotation",
    condition_col = "condition",
    assay_name = "exprs",
    markers = panel$antigen,
    folder_save = "marker_expression_by_celltype"
) {
  
  output_path <- paste0("/Volumes/HD_rafael/Cytof/OUTPUTS_RAFA_v3/", pdf_infix_name, "/", folder_save)
  dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
  
  if (!annotation_col %in% colnames(colData(sce))) {
    stop("annotation_col not found in colData(sce): ", annotation_col)
  }
  
  if (!condition_col %in% colnames(colData(sce))) {
    stop("condition_col not found in colData(sce): ", condition_col)
  }
  
  annotation_tag <- gsub("[^A-Za-z0-9_]+", "_", annotation_col)
  
  # Keep only markers present in the SCE
  markers <- markers[markers %in% rownames(assay(sce, assay_name))]
  
  ann_raw <- colData(sce)[[annotation_col]]
  ann <- as.character(ann_raw)
  condition <- as.character(colData(sce)[[condition_col]])
  
  keep <- !is.na(ann) & !is.na(condition) & condition %in% c("WT", "KO")
  
  expr_mat <- assay(sce, assay_name)[markers, keep, drop = FALSE]
  ann <- ann[keep]
  condition <- condition[keep]
  
  group <- paste(ann, condition, sep = " | ")
  
  # Keeps columns ordered by annotation level, then WT/KO.
  # This preserves numeric cluster order if annotation_col is a factor, e.g. meta20_cluster.
  if (is.factor(ann_raw)) {
    celltypes <- levels(ann_raw)
    celltypes <- celltypes[celltypes %in% ann]
  } else {
    celltypes <- sort(unique(ann))
  }
  
  condition_order <- c("WT", "KO")
  
  groups <- as.vector(
    t(outer(celltypes, condition_order, paste, sep = " | "))
  )
  groups <- groups[groups %in% unique(group)]
  
  median_mat <- vapply(groups, function(g) {
    matrixStats::rowMedians(
      expr_mat[, group == g, drop = FALSE],
      na.rm = TRUE
    )
  }, numeric(length(markers)))
  
  rownames(median_mat) <- markers
  
  # Split labels over two lines to reduce clipping
  column_labels <- colnames(median_mat)
  
  # Dynamic sizing, useful for meta20 because it has 40 WT/KO columns
  pdf_width <- max(16, length(groups) * 0.55)
  pdf_height <- max(8, length(markers) * 0.25 + 2)
  
  # Row-scaled z-score per marker across annotation-condition groups
  scaled_mat <- t(scale(t(median_mat)))
  scaled_mat[is.na(scaled_mat)] <- 0
  
  # Save scaled matrix for Cory
  write.csv(
    round(scaled_mat, 3),
    paste0(
      output_path,
      "/scaled_marker_expression_by_",
      annotation_tag,
      "_condition_",
      pdf_infix_name,
      ".csv"
    )
  )
  
  # Scaled expression heatmap, no clustering / no dendrogram
  pdf(
    paste0(
      output_path,
      "/heatmap_scaled_marker_expression_by_",
      annotation_tag,
      "_condition_",
      pdf_infix_name,
      ".pdf"
    ),
    width = pdf_width,
    height = pdf_height
  )
  
  ht_scaled <- Heatmap(
    scaled_mat,
    name = "Row-scaled\nexprs",
    col = colorRamp2(
      c(-2, 0, 2),
      c("#2166AC", "white", "#B2182B")
    ),
    cluster_rows = T,
    cluster_columns = T,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    rect_gp = gpar(col = "white"),
    column_labels = column_labels,
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = 10),
    row_names_gp = gpar(fontsize = 10),
    column_names_max_height = unit(70, "mm"),
    column_title = paste0(
      "Scaled marker expression by ",
      annotation_col,
      " and condition: ",
      pdf_infix_name
    )
  )
  
  draw(
    ht_scaled,
    heatmap_legend_side = "right",
    padding = unit(c(2, 10, 8, 2), "mm")
  )
  
  dev.off()
  
  return(list(
    median_mat = median_mat,
    scaled_mat = scaled_mat
  ))
}


plot_marker_delta_heatmap_by_condition <- function(
    sce,
    pdf_infix_name,
    annotation_col = "meta_cory_annotation",
    condition_col = "condition",
    assay_name = "exprs",
    markers = panel$antigen,
    annotation_subset = NULL,
    folder_save = "marker_expression_by_celltype",
    clip_quantile = 0.98,
    output_suffix = NULL
) {
  
  output_path <- paste0("/Volumes/HD_rafael/Cytof/OUTPUTS_RAFA_v3/", pdf_infix_name, "/", folder_save)
  dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
  
  if (!annotation_col %in% colnames(colData(sce))) {
    stop("annotation_col not found in colData(sce): ", annotation_col)
  }
  
  if (!condition_col %in% colnames(colData(sce))) {
    stop("condition_col not found in colData(sce): ", condition_col)
  }
  
  annotation_tag <- gsub("[^A-Za-z0-9_]+", "_", annotation_col)
  output_suffix_tag <- if (is.null(output_suffix) || !nzchar(output_suffix)) {
    ""
  } else {
    paste0("_", gsub("[^A-Za-z0-9_]+", "_", output_suffix))
  }
  
  expr_mat <- assay(sce, assay_name)
  requested_markers <- unique(as.character(markers))
  missing_markers <- setdiff(requested_markers, rownames(expr_mat))

  if (length(missing_markers) > 0 && !missing(markers)) {
    stop(
      "Requested markers not found in assay '",
      assay_name,
      "': ",
      paste(missing_markers, collapse = ", ")
    )
  }

  markers <- markers[markers %in% rownames(expr_mat)]

  if (length(markers) == 0) {
    stop("No requested markers were found in assay '", assay_name, "'")
  }

  ann_raw <- colData(sce)[[annotation_col]]
  ann <- as.character(ann_raw)
  condition <- as.character(colData(sce)[[condition_col]])
  
  keep <- !is.na(ann) & !is.na(condition) & condition %in% c("WT", "KO")
  
  expr_mat <- expr_mat[markers, keep, drop = FALSE]
  ann <- ann[keep]
  condition <- condition[keep]
  
  if (is.factor(ann_raw)) {
    annotations <- levels(ann_raw)
    annotations <- annotations[annotations %in% ann]
  } else {
    annotations <- sort(unique(ann))
  }
  
  if (!is.null(annotation_subset)) {
    requested_annotations <- unique(as.character(annotation_subset))
    missing_annotations <- setdiff(requested_annotations, unique(ann))

    if (length(missing_annotations) > 0) {
      stop(
        "Requested annotations not found in '",
        annotation_col,
        "': ",
        paste(missing_annotations, collapse = ", ")
      )
    }

    annotations <- requested_annotations
  }

  if (length(annotations) == 0) {
    stop("No annotations are available for the delta heatmap")
  }

  delta_mat <- matrix(
    NA_real_,
    nrow = length(markers),
    ncol = length(annotations),
    dimnames = list(markers, annotations)
  )
  
  median_wt_mat <- delta_mat
  median_ko_mat <- delta_mat
  
  for (ct in annotations) {
    
    wt_idx <- ann == ct & condition == "WT"
    ko_idx <- ann == ct & condition == "KO"
    
    if (sum(wt_idx) == 0 || sum(ko_idx) == 0) next
    
    wt_medians <- matrixStats::rowMedians(
      expr_mat[, wt_idx, drop = FALSE],
      na.rm = TRUE
    )
    
    ko_medians <- matrixStats::rowMedians(
      expr_mat[, ko_idx, drop = FALSE],
      na.rm = TRUE
    )
    
    median_wt_mat[, ct] <- wt_medians
    median_ko_mat[, ct] <- ko_medians
    delta_mat[, ct] <- ko_medians - wt_medians
  }
  
  contrast_df <- data.frame(
    marker = rep(rownames(delta_mat), times = ncol(delta_mat)),
    annotation = rep(colnames(delta_mat), each = nrow(delta_mat)),
    median_WT = as.vector(median_wt_mat),
    median_KO = as.vector(median_ko_mat),
    delta_KO_minus_WT = as.vector(delta_mat),
    stringsAsFactors = FALSE
  )
  
  write.csv(
    contrast_df,
    paste0(
      output_path,
      "/delta_KO_minus_WT_marker_expression_by_",
      annotation_tag,
      "_",
      pdf_infix_name,
      output_suffix_tag,
      ".csv"
    ),
    row.names = FALSE
  )
  
  # Raw KO - WT delta heatmap
  raw_max_abs <- as.numeric(
    quantile(abs(delta_mat), probs = clip_quantile, na.rm = TRUE)
  )
  
  if (!is.finite(raw_max_abs) || raw_max_abs == 0) {
    raw_max_abs <- max(abs(delta_mat), na.rm = TRUE)
  }
  
  if (!is.finite(raw_max_abs) || raw_max_abs == 0) {
    raw_max_abs <- 1
  }
  
  delta_mat_clipped <- pmax(pmin(delta_mat, raw_max_abs), -raw_max_abs)
  
  pdf_width <- max(12, length(annotations) * 0.85 + 4)
  pdf_height <- max(8, length(markers) * 0.28 + 2.5)
  
  pdf(
    paste0(
      output_path,
      "/heatmap_delta_KO_minus_WT_marker_expression_by_",
      annotation_tag,
      "_",
      pdf_infix_name,
      output_suffix_tag,
      ".pdf"
    ),
    width = pdf_width,
    height = pdf_height
  )
  
  ht_delta <- Heatmap(
    delta_mat_clipped,
    name = "KO - WT\nmedian exprs",
    col = colorRamp2(
      c(-raw_max_abs, 0, raw_max_abs),
      c("#2166AC", "white", "#B2182B")
    ),
    na_col = "grey90",
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    rect_gp = gpar(col = "white"),
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = 9),
    row_names_gp = gpar(fontsize = 10),
    column_names_max_height = unit(90, "mm"),
    column_title = paste0(
      "Within-cluster marker expression difference: KO - WT, ",
      pdf_infix_name
    )
  )
  
  draw(
    ht_delta,
    heatmap_legend_side = "right",
    padding = unit(c(25, 25, 8, 12), "mm")
  )
  
  dev.off()
  
  # Row-normalised KO - WT delta heatmap
  # This scales each marker row across clusters/cell types.
  print("row scaling")
  delta_mat_row_scaled <- t(scale(t(delta_mat)))
  delta_mat_row_scaled[is.na(delta_mat_row_scaled)] <- 0
  
  write.csv(
    round(delta_mat_row_scaled, 3),
    paste0(
      output_path,
      "/row_scaled_delta_KO_minus_WT_marker_expression_by_",
      annotation_tag,
      "_",
      pdf_infix_name,
      output_suffix_tag,
      ".csv"
    )
  )
  
  row_scaled_max_abs <- as.numeric(
    quantile(abs(delta_mat_row_scaled), probs = clip_quantile, na.rm = TRUE)
  )
  
  if (!is.finite(row_scaled_max_abs) || row_scaled_max_abs == 0) {
    row_scaled_max_abs <- max(abs(delta_mat_row_scaled), na.rm = TRUE)
  }
  
  if (!is.finite(row_scaled_max_abs) || row_scaled_max_abs == 0) {
    row_scaled_max_abs <- 1
  }
  
  delta_mat_row_scaled_clipped <- pmax(
    pmin(delta_mat_row_scaled, row_scaled_max_abs),
    -row_scaled_max_abs
  )
  
  pdf(
    paste0(
      output_path,
      "/heatmap_row_scaled_delta_KO_minus_WT_marker_expression_by_",
      annotation_tag,
      "_",
      pdf_infix_name,
      output_suffix_tag,
      ".pdf"
    ),
    width = pdf_width,
    height = pdf_height
  )
  
  ht_delta_row_scaled <- Heatmap(
    delta_mat_row_scaled_clipped,
    name = "Row-scaled\nKO - WT",
    col = colorRamp2(
      c(-row_scaled_max_abs, 0, row_scaled_max_abs),
      c("#2166AC", "white", "#B2182B")
    ),
    na_col = "grey90",
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    rect_gp = gpar(col = "white"),
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = 9),
    row_names_gp = gpar(fontsize = 10),
    column_names_max_height = unit(90, "mm"),
    column_title = paste0(
      "Row-normalised within-cluster marker expression difference: KO - WT, ",
      pdf_infix_name
    )
  )
  
  draw(
    ht_delta_row_scaled,
    heatmap_legend_side = "right",
    padding = unit(c(25, 25, 8, 12), "mm")
  )
  
  dev.off()
  
  return(list(
    contrast_df = contrast_df,
    delta_mat = delta_mat,
    delta_mat_clipped = delta_mat_clipped,
    delta_mat_row_scaled = delta_mat_row_scaled,
    delta_mat_row_scaled_clipped = delta_mat_row_scaled_clipped
  ))
}


plot_marker_delta_bubble_by_condition <- function(
    sce,
    pdf_infix_name,
    annotation_col = "meta_cory_annotation",
    condition_col = "condition",
    assay_name = "exprs",
    markers = panel$antigen,
    annotation_subset = NULL,
    folder_save = "marker_expression_by_celltype",
    clip_quantile = 0.98,
    output_suffix = NULL,
    max_bubble_size = 14
) {
  output_path <- paste0(
    "/Volumes/HD_rafael/Cytof/OUTPUTS_RAFA_v3/",
    pdf_infix_name,
    "/",
    folder_save
  )
  dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

  if (!annotation_col %in% colnames(colData(sce))) {
    stop("annotation_col not found in colData(sce): ", annotation_col)
  }

  if (!condition_col %in% colnames(colData(sce))) {
    stop("condition_col not found in colData(sce): ", condition_col)
  }

  annotation_tag <- gsub("[^A-Za-z0-9_]+", "_", annotation_col)
  output_suffix_tag <- if (is.null(output_suffix) || !nzchar(output_suffix)) {
    ""
  } else {
    paste0("_", gsub("[^A-Za-z0-9_]+", "_", output_suffix))
  }

  expr_mat <- assay(sce, assay_name)
  requested_markers <- unique(as.character(markers))
  missing_markers <- setdiff(requested_markers, rownames(expr_mat))

  if (length(missing_markers) > 0 && !missing(markers)) {
    stop(
      "Requested markers not found in assay '",
      assay_name,
      "': ",
      paste(missing_markers, collapse = ", ")
    )
  }

  markers <- requested_markers[requested_markers %in% rownames(expr_mat)]

  if (length(markers) == 0) {
    stop("No requested markers were found in assay '", assay_name, "'")
  }

  ann_raw <- colData(sce)[[annotation_col]]
  ann <- as.character(ann_raw)
  condition <- as.character(colData(sce)[[condition_col]])
  keep <- !is.na(ann) & !is.na(condition) & condition %in% c("WT", "KO")

  expr_mat <- expr_mat[markers, keep, drop = FALSE]
  ann <- ann[keep]
  condition <- condition[keep]

  if (is.null(annotation_subset)) {
    if (is.factor(ann_raw)) {
      annotations <- levels(ann_raw)
      annotations <- annotations[annotations %in% ann]
    } else {
      annotations <- sort(unique(ann))
    }
  } else {
    annotations <- unique(as.character(annotation_subset))
    missing_annotations <- setdiff(annotations, unique(ann))

    if (length(missing_annotations) > 0) {
      stop(
        "Requested annotations not found in '",
        annotation_col,
        "': ",
        paste(missing_annotations, collapse = ", ")
      )
    }
  }

  bubble_df <- do.call(
    rbind,
    lapply(annotations, function(annotation) {
      wt_idx <- ann == annotation & condition == "WT"
      ko_idx <- ann == annotation & condition == "KO"

      median_wt <- if (any(wt_idx)) {
        matrixStats::rowMedians(expr_mat[, wt_idx, drop = FALSE], na.rm = TRUE)
      } else {
        rep(NA_real_, length(markers))
      }

      median_ko <- if (any(ko_idx)) {
        matrixStats::rowMedians(expr_mat[, ko_idx, drop = FALSE], na.rm = TRUE)
      } else {
        rep(NA_real_, length(markers))
      }

      data.frame(
        marker = markers,
        annotation = annotation,
        median_WT = median_wt,
        median_KO = median_ko,
        delta_KO_minus_WT = median_ko - median_wt,
        n_cells_WT = sum(wt_idx),
        n_cells_KO = sum(ko_idx),
        n_cells_total = sum(wt_idx) + sum(ko_idx),
        stringsAsFactors = FALSE
      )
    })
  )

  delta_mat <- matrix(
    bubble_df$delta_KO_minus_WT,
    nrow = length(markers),
    ncol = length(annotations),
    dimnames = list(markers, annotations)
  )

  delta_mat_row_scaled <- t(scale(t(delta_mat)))
  delta_mat_row_scaled[is.na(delta_mat_row_scaled)] <- 0
  bubble_df$delta_KO_minus_WT_row_scaled <- as.vector(
    delta_mat_row_scaled
  )

  get_symmetric_limit <- function(values) {
    max_abs <- as.numeric(
      quantile(abs(values), probs = clip_quantile, na.rm = TRUE)
    )

    if (!is.finite(max_abs) || max_abs == 0) {
      max_abs <- max(abs(values), na.rm = TRUE)
    }

    if (!is.finite(max_abs) || max_abs == 0) {
      max_abs <- 1
    }

    max_abs
  }

  max_abs_delta <- get_symmetric_limit(bubble_df$delta_KO_minus_WT)
  max_abs_row_scaled <- get_symmetric_limit(
    bubble_df$delta_KO_minus_WT_row_scaled
  )

  bubble_df$marker <- factor(bubble_df$marker, levels = rev(markers))
  bubble_df$annotation <- factor(
    bubble_df$annotation,
    levels = annotations
  )

  raw_output_file_stem <- paste0(
    "bubbleplot_delta_KO_minus_WT_marker_expression_by_",
    annotation_tag,
    "_",
    pdf_infix_name,
    output_suffix_tag
  )

  row_scaled_output_file_stem <- paste0(
    "bubbleplot_row_scaled_delta_KO_minus_WT_marker_expression_by_",
    annotation_tag,
    "_",
    pdf_infix_name,
    output_suffix_tag
  )

  write.csv(
    bubble_df,
    paste0(output_path, "/", raw_output_file_stem, ".csv"),
    row.names = FALSE
  )

  make_bubble_plot <- function(
      colour_col,
      colour_limit,
      colour_title,
      plot_title
  ) {
    plot_df <- bubble_df
    plot_df$colour_value <- plot_df[[colour_col]]

    ggplot(
      plot_df,
      aes(x = annotation, y = marker)
    ) +
      geom_point(
        aes(size = n_cells_total, colour = colour_value),
        shape = 16,
        na.rm = TRUE
      ) +
      scale_colour_gradient2(
        low = "#2166AC",
        mid = "white",
        high = "#B2182B",
        midpoint = 0,
        limits = c(-colour_limit, colour_limit),
        oob = scales::squish,
        name = colour_title
      ) +
      scale_size_area(
        max_size = max_bubble_size,
        name = "Cells\n(WT + KO)"
      ) +
      labs(
        x = annotation_col,
        y = NULL,
        title = paste0(plot_title, ": ", pdf_infix_name)
      ) +
      theme_classic() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)
      )
  }

  p_raw <- make_bubble_plot(
    colour_col = "delta_KO_minus_WT",
    colour_limit = max_abs_delta,
    colour_title = "KO - WT\nmedian exprs",
    plot_title = "Marker-expression difference and cluster size"
  )

  p_row_scaled <- make_bubble_plot(
    colour_col = "delta_KO_minus_WT_row_scaled",
    colour_limit = max_abs_row_scaled,
    colour_title = "Row-scaled\nKO - WT",
    plot_title = "Row-scaled marker-expression difference and cluster size"
  )

  save_bubble_plot <- function(plot, output_file_stem) {
    pdf(
      paste0(output_path, "/", output_file_stem, ".pdf"),
      width = max(10, length(annotations) * 0.75 + 4),
      height = max(6, length(markers) * 0.65 + 2.5)
    )
    print(plot)
    dev.off()
  }

  save_bubble_plot(p_raw, raw_output_file_stem)
  save_bubble_plot(p_row_scaled, row_scaled_output_file_stem)

  return(list(
    data = bubble_df,
    raw_plot = p_raw,
    row_scaled_plot = p_row_scaled
  ))
}


# stats
extract_bio_rep_id <- function(x) {
  x <- sub("_.*$", "", x)                           # c07_Lymphocyte -> c07
  x <- sub("^([A-Za-z]+)0+([0-9]+)$", "\\1\\2", x)  # c07 -> c7
  x
}

safe_wilcox <- function(x, y) {
  tryCatch(
    suppressWarnings(wilcox.test(x, y, exact = FALSE)),
    error = function(e) NULL
  )
}

run_diffcyt_da_glmm <- function(
    sce,
    dataset_name,
    category_col,
    category_col_name = category_col,
    condition_col = "condition",
    sample_col = "sample_id",
    min_cells = 3,
    min_samples = NULL
) {
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("The lme4 package is required")
  }

  if (!requireNamespace("multcomp", quietly = TRUE)) {
    stop("The multcomp package is required")
  }

  required_cols <- c(category_col, condition_col, sample_col)
  missing_cols <- required_cols[!required_cols %in% colnames(colData(sce))]

  if (length(missing_cols) > 0) {
    stop("Columns not found in colData(sce): ", paste(missing_cols, collapse = ", "))
  }

  category_raw <- colData(sce)[[category_col]]

  meta <- data.frame(
    category = as.character(category_raw),
    condition = as.character(colData(sce)[[condition_col]]),
    bio_rep_id = extract_bio_rep_id(as.character(colData(sce)[[sample_col]])),
    stringsAsFactors = FALSE
  )

  meta <- meta[
    meta$condition %in% c("WT", "KO") &
      !is.na(meta$category) &
      !is.na(meta$bio_rep_id),
  ]

  if (is.factor(category_raw)) {
    categories <- levels(category_raw)
    categories <- categories[categories %in% meta$category]
  } else {
    categories <- sort(unique(meta$category))
  }

  sample_info <- unique(meta[c("bio_rep_id", "condition")])

  if (any(table(sample_info$bio_rep_id) > 1)) {
    stop("Each biological replicate must belong to only one condition")
  }

  sample_info$condition <- factor(
    sample_info$condition,
    levels = c("WT", "KO")
  )
  sample_info$bio_rep_id <- factor(sample_info$bio_rep_id)
  sample_info <- sample_info[
    order(sample_info$condition, sample_info$bio_rep_id),
  ]

  n_WT <- sum(sample_info$condition == "WT")
  n_KO <- sum(sample_info$condition == "KO")

  if (n_WT < 2 || n_KO < 2) {
    stop("At least two biological replicates per condition are required")
  }

  count_matrix <- table(
    factor(meta$category, levels = categories),
    factor(meta$bio_rep_id, levels = sample_info$bio_rep_id)
  )
  count_matrix <- unclass(count_matrix)
  storage.mode(count_matrix) <- "numeric"

  if (is.null(min_samples)) {
    min_samples <- ncol(count_matrix) / 2
  }

  tested <- rowSums(count_matrix >= min_cells) >= min_samples

  if (!any(tested)) {
    stop("No categories passed the minimum cell/sample filter")
  }

  tested_counts <- count_matrix[tested, , drop = FALSE]
  tested_sample_totals <- colSums(tested_counts)
  contrast <- matrix(c(0, 1), nrow = 1)

  glmm_results <- lapply(seq_len(nrow(tested_counts)), function(i) {
    y <- tested_counts[i, ] / tested_sample_totals
    data_i <- data.frame(
      y = y,
      n_cells_smp = tested_sample_totals,
      condition = sample_info$condition,
      bio_rep_id = sample_info$bio_rep_id
    )

    tryCatch({
      fit <- lme4::glmer(
        y ~ condition + (1 | bio_rep_id),
        data = data_i,
        family = "binomial",
        weights = n_cells_smp
      )
      test <- multcomp::glht(fit, linfct = contrast)
      test_summary <- summary(test)$test

      data.frame(
        category = rownames(tested_counts)[i],
        log_odds_KO_vs_WT = unname(test_summary$coefficients),
        p_value = unname(test_summary$pvalues),
        singular_fit = lme4::isSingular(fit),
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      data.frame(
        category = rownames(tested_counts)[i],
        log_odds_KO_vs_WT = NA_real_,
        p_value = NA_real_,
        singular_fit = NA,
        stringsAsFactors = FALSE
      )
    })
  }) %>%
    bind_rows()

  glmm_results$p_adj_BH <- p.adjust(glmm_results$p_value, method = "BH")

  sample_totals <- colSums(count_matrix)
  proportions <- sweep(count_matrix, 2, sample_totals, "/")
  wt_columns <- sample_info$condition == "WT"
  ko_columns <- sample_info$condition == "KO"

  results <- data.frame(
    dataset = dataset_name,
    category_col = category_col_name,
    category = rownames(count_matrix),
    test = "diffcyt-DA-GLMM",
    n_WT = n_WT,
    n_KO = n_KO,
    n_cells_WT = rowSums(count_matrix[, wt_columns, drop = FALSE]),
    n_cells_KO = rowSums(count_matrix[, ko_columns, drop = FALSE]),
    mean_WT = rowMeans(proportions[, wt_columns, drop = FALSE]),
    mean_KO = rowMeans(proportions[, ko_columns, drop = FALSE]),
    median_WT = apply(proportions[, wt_columns, drop = FALSE], 1, median),
    median_KO = apply(proportions[, ko_columns, drop = FALSE], 1, median),
    tested = tested,
    stringsAsFactors = FALSE
  )

  results$diff_KO_minus_WT <- results$mean_KO - results$mean_WT
  results$log_odds_KO_vs_WT <- NA_real_
  results$odds_ratio_KO_vs_WT <- NA_real_
  results$p_value <- NA_real_
  results$p_adj_BH <- NA_real_
  results$singular_fit <- NA

  tested_rows <- match(glmm_results$category, results$category)
  results$log_odds_KO_vs_WT[tested_rows] <- glmm_results$log_odds_KO_vs_WT
  results$odds_ratio_KO_vs_WT[tested_rows] <- exp(glmm_results$log_odds_KO_vs_WT)
  results$p_value[tested_rows] <- glmm_results$p_value
  results$p_adj_BH[tested_rows] <- glmm_results$p_adj_BH
  results$singular_fit[tested_rows] <- glmm_results$singular_fit
  results$direction <- ifelse(
    results$mean_KO >= results$mean_WT,
    "enriched_in_ko",
    "enriched_in_wt"
  )

  results <- results[order(results$p_value, na.last = TRUE), ]
  rownames(results) <- NULL

  return(results)
}

export_per_replicate_marker_expression <- function(
    sce,
    pdf_infix_name,
    annotation_col = "meta_cory_annotation",
    condition_col = "condition",
    sample_col = "sample_id",
    assay_name = "exprs",
    markers = panel$antigen,
    output_dir = "/Volumes/HD_rafael/Cytof/OUTPUTS_RAFA_v3"
) {
  required_cols <- c(annotation_col, condition_col, sample_col)
  missing_cols <- required_cols[!required_cols %in% colnames(colData(sce))]
  
  if (length(missing_cols) > 0) {
    stop("Columns not found in colData(sce): ", paste(missing_cols, collapse = ", "))
  }
  
  expr_mat <- assay(sce, assay_name)
  markers <- markers[markers %in% rownames(expr_mat)]
  
  meta <- as.data.frame(colData(sce))
  meta$cell_index <- seq_len(ncol(sce))
  meta$bio_rep_id <- extract_bio_rep_id(meta[[sample_col]])
  meta$condition <- factor(
    as.character(meta[[condition_col]]),
    levels = c("WT", "KO")
  )
  meta$annotation <- as.character(meta[[annotation_col]])
  
  meta <- meta %>%
    filter(
      condition %in% c("WT", "KO"),
      !is.na(bio_rep_id),
      !is.na(annotation)
    )
  
  replicate_expression <- lapply(markers, function(mk) {
    meta %>%
      mutate(expression = expr_mat[mk, cell_index]) %>%
      group_by(bio_rep_id, condition, annotation) %>%
      summarise(
        median_expression = median(expression, na.rm = TRUE),
        n_cells = sum(!is.na(expression)),
        .groups = "drop"
      ) %>%
      mutate(marker = mk)
  }) %>%
    bind_rows() %>%
    mutate(annotation_type = annotation_col) %>%
    select(
      annotation_type,
      annotation,
      marker,
      bio_rep_id,
      condition,
      median_expression,
      n_cells
    ) %>%
    arrange(annotation, marker, condition, bio_rep_id)
  
  output_path <- file.path(
    output_dir,
    pdf_infix_name,
    "per_replicate_marker_expression"
  )
  dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
  
  annotation_tag <- gsub("[^A-Za-z0-9_]+", "_", annotation_col)
  
  write.csv(
    replicate_expression,
    file.path(
      output_path,
      paste0(
        "per_replicate_marker_expression_by_",
        annotation_tag,
        "_",
        pdf_infix_name,
        ".csv"
      )
    ),
    row.names = FALSE
  )
  
  return(replicate_expression)
}

run_statistical_tests <- function(
    sce,
    pdf_infix_name,
    cluster_col = "meta22_cluster",
    condition_col = "condition",
    sample_col = "sample_id",
    assay_name = "exprs",
    markers = panel$antigen,
    output_dir = "/Volumes/HD_rafael/Cytof/OUTPUTS_RAFA_v3"
) {
  required_cols <- c(cluster_col, condition_col, sample_col)
  missing_cols <- required_cols[!required_cols %in% colnames(colData(sce))]
  
  if (length(missing_cols) > 0) {
    stop("Columns not found in colData(sce): ", paste(missing_cols, collapse = ", "))
  }
  
  output_path <- file.path(output_dir, pdf_infix_name, "stats")
  dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
  output_prefix <- file.path(output_path, pdf_infix_name)
  
  cluster_raw <- colData(sce)[[cluster_col]]
  
  if (is.factor(cluster_raw)) {
    clusters <- levels(cluster_raw)
    clusters <- clusters[clusters %in% as.character(cluster_raw)]
  } else {
    clusters <- sort(unique(as.character(cluster_raw)))
  }
  
  meta <- as.data.frame(colData(sce))
  meta$cell_index <- seq_len(ncol(sce))
  meta$bio_rep_id <- extract_bio_rep_id(meta[[sample_col]])
  meta$condition <- factor(
    as.character(meta[[condition_col]]),
    levels = c("WT", "KO")
  )
  meta$cluster <- as.character(meta[[cluster_col]])
  
  meta <- meta %>%
    filter(
      condition %in% c("WT", "KO"),
      !is.na(bio_rep_id),
      !is.na(cluster)
    )
  
  # Cluster proportions compared across biological replicates
  bio_rep_totals <- meta %>%
    count(bio_rep_id, condition, name = "total_cells")
  
  cluster_proportions_by_bioreplicate <- meta %>%
    count(bio_rep_id, condition, cluster, name = "n_cells") %>%
    complete(
      nesting(bio_rep_id, condition),
      cluster = clusters,
      fill = list(n_cells = 0)
    ) %>%
    left_join(bio_rep_totals, by = c("bio_rep_id", "condition")) %>%
    mutate(proportion = n_cells / total_cells)
  
  cluster_proportion_stats <- cluster_proportions_by_bioreplicate %>%
    group_by(cluster) %>%
    group_modify(~ {
      wt <- .x$proportion[.x$condition == "WT"]
      ko <- .x$proportion[.x$condition == "KO"]
      test_result <- safe_wilcox(ko, wt)
      
      tibble(
        test = "Wilcoxon rank-sum test",
        n_WT = length(wt),
        n_KO = length(ko),
        mean_WT = mean(wt, na.rm = TRUE),
        mean_KO = mean(ko, na.rm = TRUE),
        median_WT = median(wt, na.rm = TRUE),
        median_KO = median(ko, na.rm = TRUE),
        diff_KO_minus_WT = mean(ko, na.rm = TRUE) - mean(wt, na.rm = TRUE),
        statistic_W = if (is.null(test_result)) NA_real_ else unname(test_result$statistic),
        p_value = if (is.null(test_result)) NA_real_ else test_result$p.value
      )
    }) %>%
    ungroup() %>%
    mutate(
      p_adj_BH = p.adjust(p_value, method = "BH"),
      significant_p_0.05 = p_value < 0.05,
      significant_FDR_0.05 = p_adj_BH < 0.05
    ) %>%
    arrange(p_adj_BH)
  
  p_cluster_proportions <- ggplot(
    cluster_proportions_by_bioreplicate,
    aes(x = condition, y = proportion, fill = condition)
  ) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.15, size = 1.8, alpha = 0.85) +
    facet_wrap(~ cluster, scales = "free_y") +
    theme_classic() +
    labs(
      title = "Biological-replicate-level cluster proportions",
      subtitle = "Each point is one biological replicate",
      x = NULL,
      y = "Proportion of cells"
    ) +
    theme(legend.position = "none")
  
  # Marker expression compared between individual cells in each cluster
  expr_mat <- assay(sce, assay_name)
  markers <- markers[markers %in% rownames(expr_mat)]
  percell_results <- list()
  
  for (cl in clusters) {
    meta_cl <- meta %>% filter(cluster == cl)
    
    for (mk in markers) {
      df <- meta_cl %>%
        mutate(expression = expr_mat[mk, cell_index])
      
      wt <- df$expression[df$condition == "WT"]
      ko <- df$expression[df$condition == "KO"]
      test_result <- safe_wilcox(ko, wt)
      
      percell_results[[paste(cl, mk, sep = "_")]] <- data.frame(
        cluster = cl,
        marker = mk,
        test = "Wilcoxon rank-sum test",
        n_cells_WT = sum(!is.na(wt)),
        n_cells_KO = sum(!is.na(ko)),
        mean_WT = mean(wt, na.rm = TRUE),
        mean_KO = mean(ko, na.rm = TRUE),
        median_WT = median(wt, na.rm = TRUE),
        median_KO = median(ko, na.rm = TRUE),
        diff_KO_minus_WT = mean(ko, na.rm = TRUE) - mean(wt, na.rm = TRUE),
        statistic_W = if (is.null(test_result)) NA_real_ else unname(test_result$statistic),
        p_value = if (is.null(test_result)) NA_real_ else test_result$p.value
      )
    }
  }
  
  percell_expression_stats <- bind_rows(percell_results) %>%
    mutate(
      p_adj_BH = p.adjust(p_value, method = "BH"),
      significant_p_0.05 = p_value < 0.05,
      significant_FDR_0.05 = p_adj_BH < 0.05
    ) %>%
    arrange(p_adj_BH)
  
  write.csv(
    cluster_proportions_by_bioreplicate,
    paste0(output_prefix, "_cluster_proportions_by_bioreplicate.csv"),
    row.names = FALSE
  )
  
  write.csv(
    cluster_proportion_stats,
    paste0(output_prefix, "_cluster_proportion_wilcox_stats.csv"),
    row.names = FALSE
  )
  
  write.csv(
    percell_expression_stats,
    paste0(output_prefix, "_cluster_percell_marker_expression_wilcox_stats.csv"),
    row.names = FALSE
  )
  
  ggsave(
    paste0(output_prefix, "_cluster_proportion_plot.pdf"),
    p_cluster_proportions,
    width = 12,
    height = 8,
    dpi = 300
  )
  
  return(list(
    cluster_proportions_by_bioreplicate = cluster_proportions_by_bioreplicate,
    cluster_proportion_stats = cluster_proportion_stats,
    percell_expression_stats = percell_expression_stats,
    plots = list(cluster_proportion_plot = p_cluster_proportions)
  ))
}

plot_marker_violin_by_selected_metaclusters <- function(
    sce,
    pdf_infix_name,
    marker_cluster_map,
    meta_key = "meta22",
    condition_col = "condition",
    assay_name = "exprs",
    output_dir = "/Volumes/HD_rafael/Cytof/OUTPUTS_RAFA_v3",
    max_cells_per_cluster_condition = 1000,
    width = 7,
    height = 5
) {
  
  out <- file.path(
    output_dir,
    pdf_infix_name,
    "violins_selected_marker_by_meta_cluster"
  )
  dir.create(out, showWarnings = FALSE, recursive = TRUE)
  
  expr_mat <- assay(sce, assay_name)
  available_markers <- rownames(expr_mat)
  
  clean_name <- function(x) {
    gsub("[^a-z0-9]", "", tolower(x))
  }
  
  resolve_marker <- function(marker) {
    hit <- available_markers[clean_name(available_markers) == clean_name(marker)]
    
    if (length(hit) == 0) {
      stop("Marker not found in assay: ", marker)
    }
    
    if (length(hit) > 1) {
      stop("Marker name is ambiguous: ", marker, " matched ", paste(hit, collapse = ", "))
    }
    
    hit
  }
  
  n_clusters <- as.integer(gsub("meta", "", meta_key))
  
  meta_cluster <- factor(
    as.character(cluster_ids(sce, k = meta_key)),
    levels = as.character(seq_len(n_clusters))
  )
  
  condition <- factor(
    as.character(colData(sce)[[condition_col]]),
    levels = c("WT", "KO")
  )
  
  for (requested_marker in names(marker_cluster_map)) {
    
    mk <- resolve_marker(requested_marker)
    selected_clusters <- as.character(marker_cluster_map[[requested_marker]])
    
    cat("Plotting:", mk, "clusters:", paste(selected_clusters, collapse = ", "), "\n")
    
    df <- data.frame(
      expression = expr_mat[mk, ],
      meta_cluster = meta_cluster,
      condition = condition
    ) |>
      dplyr::filter(
        !is.na(meta_cluster),
        !is.na(condition),
        meta_cluster %in% selected_clusters
      )
    
    df$meta_cluster <- factor(
      as.character(df$meta_cluster),
      levels = selected_clusters
    )
    
    if (!is.null(max_cells_per_cluster_condition)) {
      set.seed(1234)
      
      df <- df |>
        dplyr::group_by(meta_cluster, condition) |>
        dplyr::mutate(.rand = runif(dplyr::n())) |>
        dplyr::arrange(.rand, .by_group = TRUE) |>
        dplyr::slice_head(n = max_cells_per_cluster_condition) |>
        dplyr::ungroup() |>
        dplyr::select(-.rand)
    }
    
    p <- ggplot(df, aes(x = meta_cluster, y = expression, fill = condition)) +
      geom_violin(
        scale = "width",
        trim = TRUE,
        position = position_dodge(width = 0.85),
        linewidth = 0.2
      ) +
      geom_boxplot(
        width = 0.12,
        outlier.shape = NA,
        position = position_dodge(width = 0.85),
        alpha = 0.6
      ) +
      scale_fill_manual(
        values = c("WT" = "#7A7A4A", "KO" = "#fff176"),
        drop = FALSE
      ) +
      theme_classic() +
      labs(
        title = paste0(mk, " expression in selected ", meta_key, " clusters"),
        subtitle = paste0(pdf_infix_name, " | WT vs KO"),
        x = meta_key,
        y = paste0(mk, " expression"),
        fill = NULL
      ) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top"
      )
    
    safe_mk <- gsub("[^A-Za-z0-9_]+", "_", mk)
    
    ggsave(
      filename = file.path(
        out,
        paste0(
          "violin_",
          safe_mk,
          "_clusters_",
          paste(selected_clusters, collapse = "_"),
          "_",
          meta_key,
          "_",
          pdf_infix_name,
          ".pdf"
        )
      ),
      plot = p,
      width = width,
      height = height
    )
  }
  
  invisible(out)
}

plot_marker_boxplot_by_selected_metaclusters <- function(
    sce,
    pdf_infix_name,
    marker_cluster_map,
    meta_key = "meta22",
    condition_col = "condition",
    assay_name = "exprs",
    output_dir = "/Volumes/HD_rafael/Cytof/OUTPUTS_RAFA_v3",
    max_cells_per_cluster_condition = 1000,
    point_size = 0.3,
    point_alpha = 0.22,
    box_width = 0.5,
    box_alpha = 0.85,
    dodge_width = 0.75,
    jitter_width = 0.12,
    width = 7,
    height = 5
) {
  
  out <- file.path(
    output_dir,
    pdf_infix_name,
    "boxplots_selected_marker_by_meta_cluster"
  )
  dir.create(out, showWarnings = FALSE, recursive = TRUE)
  
  expr_mat <- assay(sce, assay_name)
  available_markers <- rownames(expr_mat)
  
  clean_name <- function(x) {
    gsub("[^a-z0-9]", "", tolower(x))
  }
  
  resolve_marker <- function(marker) {
    hit <- available_markers[clean_name(available_markers) == clean_name(marker)]
    
    if (length(hit) == 0) {
      stop("Marker not found in assay: ", marker)
    }
    
    if (length(hit) > 1) {
      stop(
        "Marker name is ambiguous: ",
        marker,
        " matched ",
        paste(hit, collapse = ", ")
      )
    }
    
    hit
  }
  
  n_clusters <- as.integer(gsub("[^0-9]", "", meta_key))
  
  meta_cluster <- factor(
    as.character(cluster_ids(sce, k = meta_key)),
    levels = as.character(seq_len(n_clusters))
  )
  
  condition <- factor(
    as.character(colData(sce)[[condition_col]]),
    levels = c("WT", "KO")
  )
  
  for (requested_marker in names(marker_cluster_map)) {
    
    mk <- resolve_marker(requested_marker)
    selected_clusters <- as.character(marker_cluster_map[[requested_marker]])
    
    cat(
      "Plotting:",
      mk,
      "clusters:",
      paste(selected_clusters, collapse = ", "),
      "\n"
    )
    
    df <- data.frame(
      expression = expr_mat[mk, ],
      meta_cluster = meta_cluster,
      condition = condition
    ) |>
      dplyr::filter(
        !is.na(meta_cluster),
        !is.na(condition),
        meta_cluster %in% selected_clusters
      )
    
    df$meta_cluster <- factor(
      as.character(df$meta_cluster),
      levels = selected_clusters
    )
    
    if (!is.null(max_cells_per_cluster_condition)) {
      set.seed(1234)
      
      df <- df |>
        dplyr::group_by(meta_cluster, condition) |>
        dplyr::mutate(.rand = runif(dplyr::n())) |>
        dplyr::arrange(.rand, .by_group = TRUE) |>
        dplyr::slice_head(n = max_cells_per_cluster_condition) |>
        dplyr::ungroup() |>
        dplyr::select(-.rand)
    }
    
    p <- ggplot(df, aes(x = meta_cluster, y = expression)) +
      
      # Grey cell-level points, dodged around each WT / KO box
      geom_point(
        aes(colour = condition),
        alpha = point_alpha,
        size = point_size,
        position = position_jitterdodge(
          dodge.width = dodge_width,
          jitter.width = jitter_width,
          jitter.height = 0,
          seed = 1234
        ),
        show.legend = FALSE
      ) +
      
      scale_colour_manual(
        values = c(
          "WT" = "grey55",
          "KO" = "grey55"
        ),
        guide = "none"
      ) +
      
      # WT / KO boxplots using black-to-yellow style
      geom_boxplot(
        aes(fill = condition),
        width = box_width,
        outlier.shape = NA,
        position = position_dodge(width = dodge_width),
        alpha = box_alpha,
        linewidth = 0.45,
        colour = "black"
      ) +
      
      scale_fill_manual(
        values = c(
          "WT" = "#7A7A4A",
          "KO" = "#fff176"
        ),
        drop = FALSE
      ) +
      
      theme_classic() +
      labs(
        title = paste0(mk, " expression in selected ", meta_key, " clusters"),
        subtitle = paste0(pdf_infix_name, " | WT vs KO"),
        x = meta_key,
        y = paste0(mk, " expression"),
        fill = NULL
      ) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top",
        plot.title = element_text(face = "bold"),
        panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.25),
        panel.grid.major.x = element_blank()
      )
    
    safe_mk <- gsub("[^A-Za-z0-9_]+", "_", mk)
    
    ggsave(
      filename = file.path(
        out,
        paste0(
          "boxplot_",
          safe_mk,
          "_clusters_",
          paste(selected_clusters, collapse = "_"),
          "_",
          meta_key,
          "_",
          pdf_infix_name,
          ".pdf"
        )
      ),
      plot = p,
      width = width,
      height = height
    )
  }
  
  invisible(out)
}
