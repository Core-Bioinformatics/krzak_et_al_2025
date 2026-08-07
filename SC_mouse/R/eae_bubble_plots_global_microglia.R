library(ggplot2)
library(Seurat)

options(stringsAsFactors = FALSE)

# ============================================================
# Global homeostatic microglia / DAM bubble plots
#
# Updated request: keep NO EAE and EAE separate, and compare
# Ctrl vs KO within each disease-state group.
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

OUTPUT_BASENAME <- "eae_bubble_plots_global_microglia"
OUTPUT_DIR <- ensure_out_dir(OUTPUT_BASENAME)
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
    "Il6"
)

GENE_ALIASES <- c(
    "Il1ra" = "Il1rn"
)

CELLTYPE_SPECS <- list(
    list(
        panel_id = "microglia_homeostatic_global",
        label = "Global homeostatic microglia",
        celltype = "Microglia_Homeostatic"
    ),
    list(
        panel_id = "microglia_dam_global",
        label = "Global DAM",
        celltype = "Microglia_DAM"
    )
)

EAE_SPECS <- list(
    list(panel_suffix = "no_eae", value = "NO EAE", label = "NO EAE"),
    list(panel_suffix = "eae", value = "EAE", label = "EAE")
)

GENOTYPE_LEVELS <- c("Ctrl", "KO")
GENOTYPE_LABELS <- c("ctrl" = "Ctrl", "ko" = "KO")

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

build_bubble_data <- function(so, gene_lookup, group_levels, z_clip = 1) {
    expr_mat <- LayerData(so, assay = "SCT", layer = "data")[gene_lookup$actual_gene, , drop = FALSE]
    rownames(expr_mat) <- gene_lookup$actual_gene

    group_ids <- factor(as.character(so$plot_group), levels = group_levels)
    group_n <- table(group_ids)
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
    bubble_df$requested_gene <- gene_lookup$requested_gene[match(bubble_df$gene, gene_lookup$actual_gene)]
    bubble_df$avg_expression <- as.vector(mean_expr)
    bubble_df$scaled_expression <- as.vector(scaled_expr)
    bubble_df$percentage_expressed <- as.vector(pct_expr) * 100
    bubble_df$n_cells <- as.integer(group_n[as.character(bubble_df$plot_group)])
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
            x = "Genotype",
            y = "",
            colour = "Row-scaled\nexpression",
            size = "% cells"
        ) +
        theme_classic() +
        theme(
            axis.text = element_text(size = 14),
            axis.text.x = element_text(size = 14),
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

    ggsave(pdf_path, plot_obj, width = 7, height = 7)
    ggsave(png_path, plot_obj, width = 7, height = 7, dpi = 300)
    write.csv(bubble_df, csv_path, row.names = FALSE)
    write.csv(gene_lookup, lookup_path, row.names = FALSE)

    file.copy(pdf_path, file.path(PROCESSED_DIR, basename(pdf_path)), overwrite = TRUE)
    file.copy(png_path, file.path(PROCESSED_DIR, basename(png_path)), overwrite = TRUE)
    file.copy(csv_path, file.path(PROCESSED_DIR, basename(csv_path)), overwrite = TRUE)
    file.copy(lookup_path, file.path(PROCESSED_DIR, basename(lookup_path)), overwrite = TRUE)

    data.frame(
        panel_id = panel_id,
        pdf = basename(pdf_path),
        png = basename(png_path),
        csv = basename(csv_path),
        stringsAsFactors = FALSE
    )
}

html_escape <- function(x) {
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    x <- gsub(">", "&gt;", x, fixed = TRUE)
    x <- gsub("\"", "&quot;", x, fixed = TRUE)
    x
}

write_index_html <- function(plot_files) {
    panels <- unlist(lapply(seq_len(nrow(plot_files)), function(i) {
        title <- html_escape(plot_files$title[[i]])
        png <- html_escape(plot_files$png[[i]])
        pdf <- html_escape(plot_files$pdf[[i]])
        csv <- html_escape(plot_files$csv[[i]])
        c(
            "<section>",
            sprintf("<h2>%s</h2>", title),
            sprintf("<p><a href=\"%s\">PDF</a> | <a href=\"%s\">plot data CSV</a></p>", pdf, csv),
            sprintf("<img src=\"%s\" alt=\"%s bubble plot\">", png, title),
            "</section>"
        )
    }))

    html <- c(
        "<!doctype html>",
        "<html lang=\"en\">",
        "<head>",
        "<meta charset=\"utf-8\">",
        "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">",
        "<title>EAE global microglia Ctrl vs KO bubble plots</title>",
        "<style>",
        "body{font-family:Arial,sans-serif;max-width:1100px;margin:32px auto;padding:0 24px;color:#222;}",
        "h1{font-size:28px;margin-bottom:4px;} h2{font-size:20px;margin-top:36px;}",
        "p{line-height:1.45;} img{max-width:100%;height:auto;border:1px solid #ddd;}",
        "a{color:#0645ad;}",
        "</style>",
        "</head>",
        "<body>",
        "<h1>EAE global microglia Ctrl vs KO bubble plots</h1>",
        "<p>Condition-specific bubble plots for marker-defined global homeostatic microglia and DAM cells. Each plot compares Ctrl and KO within one disease state.</p>",
        "<p><a href=\"README.md\">README</a> | <a href=\"all_bubbleplot_data.csv\">all plot data CSV</a> | <a href=\"cell_counts.csv\">cell counts CSV</a></p>",
        panels,
        "</body>",
        "</html>"
    )

    index_path <- file.path(OUTPUT_DIR, "index.html")
    writeLines(html, con = index_path)
    file.copy(index_path, file.path(PROCESSED_DIR, "index.html"), overwrite = TRUE)
}

cat("Preparing annotated object...\n")
so <- load_annotated_object()

available_genes <- rownames(LayerData(so, assay = "SCT", layer = "data"))
gene_lookup <- resolve_requested_genes(
    genes = REQUESTED_GENES,
    available_genes = available_genes,
    aliases = GENE_ALIASES
)

all_plot_files <- list()
all_bubble_data <- list()
all_cell_counts <- list()

for (celltype_spec in CELLTYPE_SPECS) {
    for (eae_spec in EAE_SPECS) {
        panel_id <- paste(celltype_spec$panel_id, eae_spec$panel_suffix, "ctrl_ko", sep = "_")
        title_text <- sprintf("%s: %s Ctrl vs KO", celltype_spec$label, eae_spec$label)

        so_plot <- subset(
            so,
            subset = eae == eae_spec$value &
                condition %in% names(GENOTYPE_LABELS) &
                celltype_markers == celltype_spec$celltype
        )

        so_plot$plot_group <- factor(
            unname(GENOTYPE_LABELS[as.character(so_plot$condition)]),
            levels = GENOTYPE_LEVELS
        )

        bubble_df <- build_bubble_data(
            so = so_plot,
            gene_lookup = gene_lookup,
            group_levels = GENOTYPE_LEVELS
        )
        bubble_df$panel_id <- panel_id
        bubble_df$celltype <- celltype_spec$celltype
        bubble_df$celltype_label <- celltype_spec$label
        bubble_df$eae <- eae_spec$label
        bubble_df$genotype <- as.character(bubble_df$plot_group)

        plot_obj <- plot_bubble(bubble_df, title_text)
        plot_files <- write_plot_outputs(plot_obj, bubble_df, gene_lookup, panel_id)
        plot_files$title <- title_text
        all_plot_files[[panel_id]] <- plot_files
        all_bubble_data[[panel_id]] <- bubble_df

        count_table <- table(so_plot$plot_group)
        all_cell_counts[[panel_id]] <- data.frame(
            panel_id = panel_id,
            celltype = celltype_spec$celltype,
            celltype_label = celltype_spec$label,
            eae = eae_spec$label,
            genotype = GENOTYPE_LEVELS,
            n_cells = as.integer(count_table[GENOTYPE_LEVELS]),
            stringsAsFactors = FALSE
        )
    }
}

plot_files_df <- do.call(rbind, all_plot_files)
all_bubble_df <- do.call(rbind, all_bubble_data)
cell_counts_df <- do.call(rbind, all_cell_counts)

write.csv(all_bubble_df, file.path(OUTPUT_DIR, "all_bubbleplot_data.csv"), row.names = FALSE)
write.csv(cell_counts_df, file.path(OUTPUT_DIR, "cell_counts.csv"), row.names = FALSE)
file.copy(file.path(OUTPUT_DIR, "all_bubbleplot_data.csv"), file.path(PROCESSED_DIR, "all_bubbleplot_data.csv"), overwrite = TRUE)
file.copy(file.path(OUTPUT_DIR, "cell_counts.csv"), file.path(PROCESSED_DIR, "cell_counts.csv"), overwrite = TRUE)

readme_path <- file.path(OUTPUT_DIR, "README.md")
writeLines(
    c(
        "# EAE Global Microglia Ctrl-vs-KO Bubble Plots",
        "",
        sprintf("**Date**: %s", Sys.Date()),
        "",
        "## What This Folder Contains",
        "",
        "This folder contains condition-specific Ctrl-vs-KO bubble plots for two global marker-celltype groups:",
        "",
        "- homeostatic microglia: `celltype_markers == \"Microglia_Homeostatic\"`",
        "- DAM: `celltype_markers == \"Microglia_DAM\"`",
        "",
        "For each marker-celltype group, the script writes one plot for `NO EAE` and one plot for `EAE`.",
        "",
        "## Data Used",
        "",
        "- Base object: `seurat_object.rds`",
        "- Stable clusters: from `clustassess_object.rds` (Most Abundant / 1,950 / SLM / 20)",
        "- Expression layer: SCT `data`",
        "- Cell-type labels: `celltype_markers` from `marker_genes.R`",
        "",
        "## How To Read The Plots",
        "",
        "- subset: cells assigned to the target marker-based cell type and disease state",
        "- x-axis: genotype (`Ctrl`, `KO`) within that disease state",
        "- dot size: percent of cells expressing the gene",
        "- dot color: mean expression scaled by gene across the displayed Ctrl and KO groups",
        "- removed genes from the earlier broader request: `Arg1`, `Nos2`",
        "",
        "## Files",
        "",
        "- `index.html`: browser index with all PNG plots and links",
        "- `*_bubbleplot.pdf` and `*_bubbleplot.png`: one plot per marker-celltype and disease state",
        "- `*_bubbleplot_data.csv`: plotting data used to make each figure",
        "- `all_bubbleplot_data.csv`: combined plotting data",
        "- `cell_counts.csv`: cells per marker-celltype, disease state, and genotype",
        "- `*_gene_lookup.csv`: requested and actual gene symbols used in each figure",
        "",
        "## Cell Counts",
        "",
        paste(capture.output(print(cell_counts_df, row.names = FALSE)), collapse = "\n"),
        "",
        "## Script",
        "",
        "Generated by the corresponding script in `R/`",
        "",
        "```bash",
        "Rscript R/<script>.R",
        "```"
    ),
    con = readme_path
)
file.copy(readme_path, file.path(PROCESSED_DIR, "README.md"), overwrite = TRUE)
write_index_html(plot_files_df)


cat("\nDone.\n")

