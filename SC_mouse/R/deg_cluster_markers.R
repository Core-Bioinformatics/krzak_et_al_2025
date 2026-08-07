library(Seurat)
library(dplyr)

options(future.globals.maxSize = 4 * 1024^3)  # 4 GiB

# ============================================================
# Cluster marker genes: cluster vs complement
#
# For each EAE filter level (ALL, EAE, NO EAE):
#   For each cluster in CLUSTER_COL (default: stable_20_clusters):
#     FindMarkers(ident.1 = cluster) — one cluster vs all other cells
#     Keep only positive log2FC (genes upregulated in the cluster)
#     Report top N by log2FC
#
# Output: one CSV per cluster/filter + combined top-N CSV + Markdown report.
# ============================================================

# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

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

OUTPUT_DIR <- ensure_out_dir("deg_cluster_markers")

# FindMarkers parameters
TEST_USE <- "wilcox"
LOGFC_THRESHOLD <- 0.25
MIN_PCT <- 0.1

# Top N marker genes per cluster (by positive log2FC)
TOP_N <- 5

# Minimum cells in cluster to run (below this, flag but still run)
MIN_CELLS_WARNING <- 20

# Clustering column to use
CLUSTER_COL <- "stable_20_clusters"

# Clusters to analyse: NULL = auto-detect from data after loading
CLUSTERS_OF_INTEREST <- NULL

# EAE filter levels: "ALL" means no filter, then EAE-only, then NO EAE-only
EAE_FILTERS <- c("ALL", "EAE", "NO EAE")

# Sanitise filter name for filenames
sanitise_filter <- function(filter_name) gsub(" ", "", filter_name)

# ------------------------------------------------------------
# Load data (annotated Seurat object + marker voting)
# ------------------------------------------------------------

source_methods("load_annotated_object.R")
cat("Loading annotated Seurat object...\n")
so <- load_annotated_object()
cat(sprintf("  %d cells x %d features\n", ncol(so), nrow(so)))
gc()

# ------------------------------------------------------------
# Helper: run cluster-vs-complement for one cluster
# ------------------------------------------------------------

run_cluster_markers <- function(so_sub, cluster_id, filter_label) {
    Idents(so_sub) <- CLUSTER_COL

    n_cluster <- sum(Idents(so_sub) == cluster_id, na.rm = TRUE)
    n_complement <- sum(Idents(so_sub) != cluster_id, na.rm = TRUE)

    label <- sprintf("%s_cluster_%s", sanitise_filter(filter_label), cluster_id)
    cat(sprintf("  %s: cluster=%d, complement=%d", label, n_cluster, n_complement))

    if (n_cluster == 0) {
        cat(" -> SKIPPED (0 cells in cluster)\n")
        return(NULL)
    }

    low_n <- n_cluster < MIN_CELLS_WARNING
    if (low_n) cat(sprintf(" [WARNING: low n=%d]", n_cluster))
    cat("\n")

    tryCatch({
        degs <- FindMarkers(
            so_sub,
            ident.1 = cluster_id,
            test.use = TEST_USE,
            logfc.threshold = LOGFC_THRESHOLD,
            min.pct = MIN_PCT,
            only.pos = TRUE,
            verbose = FALSE
        )

        if (nrow(degs) == 0) {
            cat("    -> 0 markers found\n")
            return(list(
                full = NULL,
                summary = data.frame(
                    filter = filter_label, cluster = cluster_id,
                    n_cluster = n_cluster, n_complement = n_complement,
                    n_markers_total = 0, n_markers_sig = 0,
                    low_n = low_n, stringsAsFactors = FALSE
                )
            ))
        }

        degs$gene <- rownames(degs)
        degs$filter <- filter_label
        degs$cluster <- cluster_id
        degs$n_cluster <- n_cluster
        degs$n_complement <- n_complement
        degs$low_n_warning <- low_n

        # Save full CSV
        fname <- sprintf("markers_%s.csv", label)
        fpath <- file.path(OUTPUT_DIR, fname)
        write.csv(degs, fpath, row.names = FALSE)
        cat(sprintf("    -> %d markers (pos log2FC), saved to %s\n",
            nrow(degs), basename(fpath)))

        sig <- degs[degs$p_val_adj < 0.05, ]

        list(
            full = degs,
            summary = data.frame(
                filter = filter_label, cluster = cluster_id,
                n_cluster = n_cluster, n_complement = n_complement,
                n_markers_total = nrow(degs), n_markers_sig = nrow(sig),
                low_n = low_n, stringsAsFactors = FALSE
            )
        )
    }, error = function(e) {
        cat(sprintf("    -> ERROR: %s\n", conditionMessage(e)))
        list(
            full = NULL,
            summary = data.frame(
                filter = filter_label, cluster = cluster_id,
                n_cluster = n_cluster, n_complement = n_complement,
                n_markers_total = NA, n_markers_sig = NA,
                low_n = low_n, stringsAsFactors = FALSE
            )
        )
    })
}

# ============================================================
# Main loop: iterate over EAE filter levels
# ============================================================

summary_rows <- list()
top_n_rows <- list()

for (eae_filter in EAE_FILTERS) {
    filter_tag <- sanitise_filter(eae_filter)

    cat(sprintf("\n##########################################################\n"))
    cat(sprintf("# Filter: %s\n", eae_filter))
    cat(sprintf("##########################################################\n"))

    if (eae_filter == "ALL") {
        so_sub <- so
    } else {
        so_sub <- subset(so, eae == eae_filter)
    }
    cat(sprintf("  Subset: %d cells\n", ncol(so_sub)))

    for (cl in CLUSTERS_OF_INTEREST) {
        res <- run_cluster_markers(so_sub, cl, eae_filter)
        if (!is.null(res)) {
            summary_rows[[length(summary_rows) + 1]] <- res$summary

            if (!is.null(res$full) && nrow(res$full) > 0) {
                top <- res$full %>%
                    filter(p_val_adj < 0.05) %>%
                    arrange(desc(avg_log2FC)) %>%
                    head(TOP_N)

                if (nrow(top) > 0) {
                    top_n_rows[[length(top_n_rows) + 1]] <- top
                }
            }
        }
    }

    if (eae_filter != "ALL") {
        rm(so_sub)
        gc()
    }
}

rm(so)
gc()

# ============================================================
# Summary table
# ============================================================

summary_df <- bind_rows(summary_rows)
summary_path <- file.path(OUTPUT_DIR, "markers_summary.csv")
write.csv(summary_df, summary_path, row.names = FALSE)

cat("\n========================================\n")
cat("SUMMARY\n")
cat("========================================\n\n")
print(summary_df, right = FALSE)

# ============================================================
# Combined top-N table
# ============================================================

if (length(top_n_rows) > 0) {
    top_n_df <- bind_rows(top_n_rows)
    top_n_path <- file.path(OUTPUT_DIR, "markers_top_genes.csv")
    write.csv(top_n_df, top_n_path, row.names = FALSE)
    cat(sprintf("\nTop %d markers per cluster saved to %s\n", TOP_N, top_n_path))
} else {
    top_n_df <- data.frame()
    cat("\nNo significant markers found.\n")
}

# ============================================================
# Markdown report
# ============================================================

cat("\nGenerating Markdown report...\n")

report_path <- file.path(OUTPUT_DIR, "markers_report.md")
report <- file(report_path, open = "w")

writeLines(c(
    "# Cluster Marker Genes Report (cluster vs complement)",
    "",
    sprintf("**Date**: %s", Sys.Date()),
    "",
    "**Dataset**: `seurat_object.rds`",
    "",
    "**Metadata source**: ClustAssess stable partition + marker-gene voting",
    "",
    "## Methods",
    "",
    sprintf("- **Clustering**: `%s`", CLUSTER_COL),
    sprintf("- **Test**: Wilcoxon rank-sum (`test.use = \"%s\"`)", TEST_USE),
    sprintf("- **Log2FC threshold**: %s", LOGFC_THRESHOLD),
    sprintf("- **Min percent expressed**: %s", MIN_PCT),
    "- **Direction**: only positive log2FC (genes upregulated in the cluster vs all other cells)",
    sprintf("- **Top N**: %d genes per cluster, ranked by log2FC, filtered to adj. p-value < 0.05", TOP_N),
    "",
    "### EAE filter levels",
    "",
    "- **ALL**: no filtering, all cells included",
    "- **EAE**: only cells from EAE mice",
    "- **NO EAE**: only cells from non-EAE mice",
    ""
), report)

for (eae_filter in EAE_FILTERS) {
    filter_tag <- sanitise_filter(eae_filter)
    filter_rows <- summary_df[summary_df$filter == eae_filter, ]
    filter_top <- if (nrow(top_n_df) > 0) {
        top_n_df[top_n_df$filter == eae_filter, ]
    } else {
        data.frame()
    }

    writeLines(c(
        "---",
        "",
        sprintf("# Filter: %s", eae_filter),
        "",
        "## Summary",
        "",
        "| Cluster | n cluster | n complement | Markers total | Sig (adj.p<0.05) | Low n |",
        "|---------|-----------|--------------|---------------|------------------|-------|"
    ), report)

    for (i in seq_len(nrow(filter_rows))) {
        r <- filter_rows[i, ]
        writeLines(sprintf("| %s | %d | %d | %s | %s | %s |",
            r$cluster, r$n_cluster, r$n_complement,
            ifelse(is.na(r$n_markers_total), "ERROR", as.character(r$n_markers_total)),
            ifelse(is.na(r$n_markers_sig), "-", as.character(r$n_markers_sig)),
            ifelse(r$low_n, "YES", "no")
        ), report)
    }
    writeLines("", report)

    # Top genes table per filter
    writeLines(c(
        sprintf("## Top %d marker genes per cluster [%s]", TOP_N, eae_filter),
        ""
    ), report)

    for (cl in CLUSTERS_OF_INTEREST) {
        cl_top <- if (nrow(filter_top) > 0) {
            filter_top[filter_top$cluster == cl, ]
        } else {
            data.frame()
        }

        if (nrow(cl_top) > 0) {
            writeLines(c(
                sprintf("### Cluster %s", cl),
                "",
                "| Gene | avg_log2FC | pct.1 (cluster) | pct.2 (complement) | p_val_adj |",
                "|------|------------|-----------------|--------------------|-----------| "
            ), report)

            for (j in seq_len(nrow(cl_top))) {
                d <- cl_top[j, ]
                writeLines(sprintf("| %s | %.3f | %.3f | %.3f | %.2e |",
                    d$gene, d$avg_log2FC, d$pct.1, d$pct.2, d$p_val_adj
                ), report)
            }
            writeLines("", report)
        }
    }
}

writeLines(c(
    "---",
    "",
    "## Output files",
    "",
    "| File | Description |",
    "|------|-------------|",
    "| `markers_summary.csv` | Summary table (one row per cluster per filter) |",
    sprintf("| `markers_top_genes.csv` | Combined top %d marker genes per cluster (sig only) |", TOP_N),
    "| `markers_<filter>_cluster_<N>.csv` | Full marker gene list for cluster N under filter |",
    "| `markers_report.md` | This report |",
    ""
), report)

close(report)
cat(sprintf("Report saved to %s\n", report_path))



cat("\nDone.\n")
