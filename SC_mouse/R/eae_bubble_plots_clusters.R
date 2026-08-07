library(ggplot2)
library(Seurat)

options(stringsAsFactors = FALSE)

# ============================================================
# EAE WT vs KO bubble plots for requested cluster panels
#
# Dot size: percent expressing
# Dot color: row-scaled mean SCT expression across displayed groups
# ============================================================

local({
  ca <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", grep("^--file=", ca, value = TRUE))
  paths_file <- if (length(f) == 1) {
    file.path(dirname(normalizePath(f)), "00_paths.R")
  } else {
    file.path("00_paths.R")
  }
  if (!file.exists(paths_file)) {
    for (cand in c("R/00_paths.R", "00_paths.R")) {
      if (file.exists(cand)) {
        paths_file <- cand
        break
      }
    }
  }
  source(paths_file)
})

OUTPUT_DIR <- ensure_out_dir("eae_bubble_plots_clusters")
PROCESSED_DIR <- OUTPUT_DIR

source_methods("load_annotated_object.R")

REQUESTED_GENES <- c(
    "Ccl2",
    "Ccl3",
    "Ccl6",
    "Ccl12",
    "Cxcl2",
    "Cxcl9",
    "Cxcl10",
    "Il1b",
    "Il1a",
    "Il1ra",
    "Il6",
    "Arg1",
    "Nos"
)

GENE_ALIASES <- c(
    "Il1ra" = "Il1rn",
    "Nos" = "Nos2"
)

PLOT_SPECS <- list(
    list(
        panel_id = "microglia",
        title = "EAE microglia clusters: WT vs KO",
        clusters = c("1", "2", "5", "13", "16")
    ),
    list(
        panel_id = "monocyte_macrophages",
        title = "EAE monocyte/macrophage clusters: WT vs KO",
        clusters = c("6", "15", "18")
    )
)

resolve_requested_genes <- function(genes, available_genes, aliases = NULL) {
    actual_genes <- vapply(genes, function(gene_name) {
        if (gene_name %in% available_genes) {
            return(gene_name)
        }
        if (!is.null(aliases) && gene_name %in% names(aliases) && aliases[[gene_name]] %in% available_genes) {
            return(aliases[[gene_name]])
        }
        NA_character_
    }, character(1))

    missing_genes <- genes[is.na(actual_genes)]
    if (length(missing_genes) > 0) {
        stop("Missing requested genes after alias lookup: ", paste(missing_genes, collapse = ", "))
    }

    data.frame(
        requested_gene = genes,
        actual_gene = unname(actual_genes),
        stringsAsFactors = FALSE
    )
}

build_group_levels <- function(clusters) {
    as.vector(rbind(
        paste(clusters, "WT"),
        paste(clusters, "KO")
    ))
}

build_bubble_data <- function(so, gene_lookup, group_levels, z_clip = 1) {
    expr_mat <- LayerData(so, assay = "SCT", layer = "data")[gene_lookup$actual_gene, , drop = FALSE]
    rownames(expr_mat) <- gene_lookup$actual_gene

    group_ids <- factor(as.character(so$plot_group), levels = group_levels)
    mean_expr <- matrix(
        0,
        nrow = nrow(expr_mat),
        ncol = length(group_levels),
        dimnames = list(rownames(expr_mat), group_levels)
    )
    pct_expr <- mean_expr

    for (group_id in group_levels) {
        group_mask <- group_ids == group_id
        if (!any(group_mask)) next
        group_mat <- expr_mat[, group_mask, drop = FALSE]
        mean_expr[, group_id] <- rowMeans(group_mat)
        pct_expr[, group_id] <- rowMeans(group_mat > 0)
    }

    scaled_expr <- t(scale(t(mean_expr)))
    scaled_expr[is.na(scaled_expr)] <- 0
    scaled_expr[scaled_expr > z_clip] <- z_clip
    scaled_expr[scaled_expr < -z_clip] <- -z_clip

    bubble_df <- expand.grid(
        gene = rownames(mean_expr),
        plot_group = colnames(mean_expr),
        stringsAsFactors = FALSE
    )
    bubble_df$avg_expression <- as.vector(mean_expr)
    bubble_df$scaled_expression <- as.vector(scaled_expr)
    bubble_df$percentage_expressed <- as.vector(pct_expr) * 100
    bubble_df$gene <- factor(
        bubble_df$gene,
        levels = rev(gene_lookup$actual_gene)
    )
    bubble_df$plot_group <- factor(bubble_df$plot_group, levels = group_levels)
    bubble_df
}

plot_bubble <- function(bubble_df, title_text) {
    ggplot(
        bubble_df,
        aes(
            x = .data$plot_group,
            y = .data$gene,
            colour = .data$scaled_expression,
            size = .data$percentage_expressed
        )
    ) +
        geom_point() +
        scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
        scale_size_continuous(range = c(1, 10)) +
        labs(
            title = title_text,
            x = "",
            y = "",
            colour = "Row-scaled\nexpression",
            size = "% cells"
        ) +
        theme_classic() +
        theme(
            axis.text = element_text(size = 14),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title = element_text(size = 14),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 13),
            plot.title = element_text(size = 16, face = "bold")
        )
}

write_plot_outputs <- function(plot_obj, bubble_df, gene_lookup, panel_id) {
    pdf_path <- file.path(OUTPUT_DIR, paste0(panel_id, "_bubbleplot.pdf"))
    png_path <- file.path(OUTPUT_DIR, paste0(panel_id, "_bubbleplot.png"))
    csv_path <- file.path(OUTPUT_DIR, paste0(panel_id, "_bubbleplot_data.csv"))
    lookup_path <- file.path(OUTPUT_DIR, paste0(panel_id, "_gene_lookup.csv"))

    ggsave(pdf_path, plot_obj, width = 10, height = 7)
    ggsave(png_path, plot_obj, width = 10, height = 7, dpi = 300)
    write.csv(bubble_df, csv_path, row.names = FALSE)
    write.csv(data.frame(gene = gene_lookup$actual_gene), lookup_path, row.names = FALSE)

    file.copy(pdf_path, file.path(PROCESSED_DIR, basename(pdf_path)), overwrite = TRUE)
    file.copy(png_path, file.path(PROCESSED_DIR, basename(png_path)), overwrite = TRUE)
    file.copy(csv_path, file.path(PROCESSED_DIR, basename(csv_path)), overwrite = TRUE)
    file.copy(lookup_path, file.path(PROCESSED_DIR, basename(lookup_path)), overwrite = TRUE)
}

cat("Preparing annotated object...\n")
so <- load_annotated_object(
    add_marker_celltypes = FALSE
)

available_genes <- rownames(LayerData(so, assay = "SCT", layer = "data"))
gene_lookup <- resolve_requested_genes(
    genes = REQUESTED_GENES,
    available_genes = available_genes,
    aliases = GENE_ALIASES
)

for (spec in PLOT_SPECS) {
    so_plot <- subset(
        so,
        subset = eae == "EAE" & stable_20_clusters %in% spec$clusters & condition %in% c("ctrl", "ko")
    )

    group_levels <- build_group_levels(spec$clusters)
    so_plot$plot_group <- factor(
        paste(
            as.character(so_plot$stable_20_clusters),
            ifelse(so_plot$condition == "ctrl", "WT", "KO")
        ),
        levels = group_levels
    )

    bubble_df <- build_bubble_data(
        so = so_plot,
        gene_lookup = gene_lookup,
        group_levels = group_levels
    )
    plot_obj <- plot_bubble(bubble_df, spec$title)
    write_plot_outputs(plot_obj, bubble_df, gene_lookup, spec$panel_id)
}

readme_path <- file.path(OUTPUT_DIR, "README.md")
writeLines(
    c(
        "# EAE Bubble Plots",
        "",
        sprintf("**Date**: %s", Sys.Date()),
        "",
        "## What This Folder Contains",
        "",
        "This folder contains two EAE-only bubble plots:",
        "",
        "- microglia clusters: `1`, `2`, `5`, `13`, `16`",
        "- monocyte/macrophage clusters: `6`, `15`, `18`",
        "",
        "## Data Used",
        "",
        "- Base object: `seurat_object.rds`",
        "- Stable clusters: from `clustassess_object.rds` (Most Abundant / 1,950 / SLM / 20)",
        "- Expression layer: SCT `data`",
        "",
        "## How To Read The Plots",
        "",
        "- subset: `EAE` only",
        "- x-axis: stable cluster plus genotype (`WT`, `KO`)",
        "- dot size: percent of cells expressing the gene",
        "- dot color: mean expression scaled by gene across the displayed groups",
        "",
        "## Files",
        "",
        "- `microglia_bubbleplot.pdf` and `microglia_bubbleplot.png`",
        "- `monocyte_macrophages_bubbleplot.pdf` and `monocyte_macrophages_bubbleplot.png`",
        "- `*_bubbleplot_data.csv`: plotting data used to make each figure",
        "- `*_gene_lookup.csv`: gene symbols used in each figure"
    ),
    con = readme_path
)
file.copy(readme_path, file.path(PROCESSED_DIR, "README.md"), overwrite = TRUE)

cat("Done.\n")
