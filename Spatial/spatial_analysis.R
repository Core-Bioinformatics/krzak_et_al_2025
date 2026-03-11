library(ggplot2)
library(dplyr)
library(Matrix)
library(matrixStats)
library(FNN)
library(tibble)
library(viridis)
library(ggpubr)

setwd("/Users/rafael/Desktop/RA/krzak_et_al_2025/Spatial")
source("spatial_helper_functions.R")

tissues_inputs <- readRDS("objects/tissues_inputs.rds")

names(tissues_inputs)

genes_friends <- c(
  "ARG1", "FIZZ1", "MGL1", "MGL2", "HIF1A", "IL6", "IL1B", "IL2", "CCR7", "CCR5",
  "NLRP1", "P2RY12", "FUCA1", "KLF4", "TLR7", "HRH1", "MMP12", "CCL7", "MMP1", "CD1E",
  "IL1R2", "CLEC5A", "ADORA1", "AI182371", "AI464131", "AQP11", "DLL1", "ESRP2",
  "FAM131B", "FUZ", "GDPD4", "GGACT", "H1FX", "IGSF11", "KCNA6", "KCTD11", "KLHL25",
  "LCA5L", "LGI3", "MMGT2", "MPP2", "PAQR8", "RGS9BP", "SERPINB8", "SLC16A6", "SORBS2",
  "TXNIP", "PTGS2", "NPAS4", "CYR61", "PHXR4", "KCNQ1OT1", "NFIL3", "TRIM29", "RGS16",
  "FOXD4", "4931406H21RIK", "CALHM2", "KMT2D", "NLRP1C-PS", "IL2", "GM438", "KCNK18",
  "HPCAL4", "OGN", "DLC1", "EGR2", "ZFP42", "ANKRD50", "TMEM132C", "APLNR", "TDRD1"
)

theme_large_text <- theme(
  legend.title = element_text(size = 16),
  legend.text = element_text(size = 13),
  axis.title = element_text(size = 16),
  axis.text = element_text(size = 13),
  title = element_text(size=16)
)

pathway_gene_lists <- list(sucnr1_friends = list(sucnr1_friends = genes_friends))

tissues_to_run <- names(tissues_inputs)

k <- 10
pct_thresh_active_neighbors <- 0.9
min_hotspots <- 5
aggregation <- "mean"

res_all <- run_activation_neighborhood_approach(
  tissues = tissues_to_run,
  pathway_genes_list = pathway_gene_lists$sucnr1_friends,
  k = k,
  pct_thresh_active_neighbors = pct_thresh_active_neighbors,
  min_hotspots = min_hotspots,
  tissue_inputs_map = tissues_inputs,
  output_dir = NULL,
  save_pdf = FALSE
)

# Plot mean expression per niche and then mean expression per tissue (facetted by niche)
for (aggregation in c("mean", "sum")) {
  pdf(paste0("figures/", aggregation,"_expression_summary_plot.pdf"), width = 16, height = 8)

  print(plot_compare_avg_niche_by_sample(
    tissues_selected = tissues_to_run,
    res_all = res_all,
    pathway_name = "sucnr1_friends",
    tissue_inputs_map = tissues_inputs,
    aggregation = aggregation,
    harmonize_genes = "per_tissue"
  ) + theme_large_text)

  print(plot_compare_avg_niche_by_niche(
    tissues_selected = tissues_to_run,
    res_all = res_all,
    pathway_name = "sucnr1_friends",
    tissue_inputs_map = tissues_inputs,
    aggregation = aggregation,
  ) + theme_large_text)
  dev.off()
}


# Statistical analysis
get_sample_level_agg_df <- function(
    tissues_selected,
    res_all,
    pathway_name,
    tissue_inputs_map,
    aggregation = c("mean", "sum", "median"),
    harmonize_genes = c("per_tissue", "common", "union"),
    drop_na_niche = TRUE
) {
  aggregation     <- match.arg(aggregation)
  harmonize_genes <- match.arg(harmonize_genes)

  # 1. Which genes passed in each tissue
  genes_passed_by_tissue <- lapply(tissues_selected, function(t) {
    gp <- character(0)
    if (!is.null(res_all[[t]]) &&
        !is.null(res_all[[t]][[pathway_name]]) &&
        !is.null(res_all[[t]][[pathway_name]]$summary)) {
      gp <- res_all[[t]][[pathway_name]]$summary$gene
    }
    gp
  })
  names(genes_passed_by_tissue) <- tissues_selected

  print(genes_passed_by_tissue)

  gene_set_common <- NULL
  gene_set_union  <- NULL
  if (harmonize_genes == "common") {
    gene_set_common <- Reduce(intersect, genes_passed_by_tissue)
  } else if (harmonize_genes == "union") {
    gene_set_union <- Reduce(union, genes_passed_by_tissue)
  }

  build_sample_niche <- function(t) {
    gp <- genes_passed_by_tissue[[t]]
    ti <- tissue_inputs_map[[t]]
    if (is.null(ti) || is.null(ti$expr_mat) || is.null(ti$metadata)) return(NULL)

    expr_mat <- ti$expr_mat
    md       <- ti$metadata

    genes_use <- switch(
      harmonize_genes,
      per_tissue = gp,
      common     = gene_set_common,
      union      = gene_set_union
    )
    genes_use <- genes_use[genes_use %in% rownames(expr_mat)]
    if (length(genes_use) == 0) return(NULL)

    expr_sub <- expr_mat[genes_use, , drop = FALSE]
    if (ncol(expr_sub) == 0) return(NULL)

    # Aggregate genes per spot
    agg_value_per_cell <- switch(
      aggregation,
      mean   = as.numeric(Matrix::colMeans(expr_sub)),
      sum    = as.numeric(Matrix::colSums(expr_sub)),
      median = as.numeric(matrixStats::colMedians(as.matrix(expr_sub)))
    )

    df_cells <- data.frame(
      cell     = colnames(expr_sub),
      agg_value = agg_value_per_cell,
      tissue   = t,
      stringsAsFactors = FALSE
    )

    # Attach niche, handling list-column cases
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

    if (drop_na_niche) {
      df_cells <- df_cells[!is.na(df_cells$niche), , drop = FALSE]
      if (nrow(df_cells) == 0) return(NULL)
    }

    # Collapse spots to one number per tissue × niche
    df_sample_niche <- df_cells %>%
      dplyr::group_by(niche) %>%
      dplyr::summarise(
        agg_value = mean(agg_value, na.rm = TRUE),
        n_cells   = dplyr::n(),
        .groups   = "drop"
      ) %>%
      dplyr::mutate(tissue = t) %>%
      dplyr::select(tissue, niche, agg_value, n_cells)

    df_sample_niche
  }

  parts <- lapply(tissues_selected, build_sample_niche)
  parts <- Filter(Negate(is.null), parts)
  if (length(parts) == 0) {
    return(data.frame())
  }

  df_samples <- dplyr::bind_rows(parts)

  df_samples <- df_samples %>% dplyr::filter(niche != "GM")

  df_samples
}


df_samples_mean <- get_sample_level_agg_df(
  tissues_selected  = tissues_to_run,
  res_all           = res_all,
  pathway_name      = "sucnr1_friends",
  tissue_inputs_map = tissues_inputs,
  aggregation       = "mean"
)
df_samples_mean

df_samples_sum <- get_sample_level_agg_df(
  tissues_selected  = tissues_to_run,
  res_all           = res_all,
  pathway_name      = "sucnr1_friends",
  tissue_inputs_map = tissues_inputs,
  aggregation       = "sum"
)

df_samples_sum


lesion_niches <- c("PLWM", "core", "rim")

lesion_wide <- df_samples_sum %>%
  filter(niche %in% lesion_niches) %>%
  tidyr::pivot_wider(
    id_cols = tissue,
    names_from = niche,
    values_from = agg_value
  )

# Keep tissues that have all three values
lesion_wide <- lesion_wide %>%
  filter(!is.na(PLWM) & !is.na(core) & !is.na(rim))

head(lesion_wide)
friedman_global <- friedman.test(as.matrix(lesion_wide[, c("PLWM", "core", "rim")]))

friedman_global

wilcox_plwm_core <- wilcox.test(lesion_wide$PLWM, lesion_wide$core, paired = TRUE)
wilcox_plwm_rim  <- wilcox.test(lesion_wide$PLWM, lesion_wide$rim,  paired = TRUE)
wilcox_rim_core  <- wilcox.test(lesion_wide$rim,  lesion_wide$core, paired = TRUE)

wilcox_plwm_core
wilcox_plwm_rim
wilcox_rim_core

p_raw <- c(
  plwm_vs_core = wilcox_plwm_core$p.value,
  plwm_vs_rim  = wilcox_plwm_rim$p.value,
  rim_vs_core  = wilcox_rim_core$p.value
)

# Apply Benjamini-Hochberg (FDR) correction
p_adj <- p.adjust(p_raw, method = "BH")

results <- data.frame(
  comparison = names(p_raw),
  p_raw      = unname(p_raw),
  p_adj      = unname(p_adj),
  stringsAsFactors = FALSE
)

results$signif <- ifelse(results$p_adj < 0.001, "***",
                         ifelse(results$p_adj < 0.01, "**",
                                ifelse(results$p_adj < 0.05, "*", "ns")))

print(results)

