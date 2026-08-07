library(Seurat)
library(dplyr)

# Allow large objects to be passed to future workers (needed for large subsets)
options(future.globals.maxSize = 4 * 1024^3)  # 4 GiB

# ============================================================
# DEG analysis: CTRL vs KO (per EAE status)
#
# For each EAE level (EAE, NO EAE):
#   1. Cluster-level (stable_20_clusters): all clusters
#   2. Cell-type-level (celltype_markers): all viable cell types
#
# Uses Wilcoxon rank-sum test via Seurat::FindMarkers().
# Output: one CSV per comparison + a summary CSV + Markdown report.
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

OUTPUT_DIR <- ensure_out_dir("deg_ctrl_vs_ko")

# FindMarkers parameters
TEST_USE <- "wilcox"
LOGFC_THRESHOLD <- 0.25
MIN_PCT <- 0.1

# Minimum cells per group to run DEG (below this, flag but still run)
MIN_CELLS_WARNING <- 20

# Clusters: run all (sorted numerically)
CLUSTERS_OF_INTEREST <- as.character(1:20)

# Cell types to always skip
CELLTYPES_SKIP <- c("Unassigned")

# EAE levels to process
EAE_LEVELS <- c("EAE", "NO EAE")

# Sanitise EAE level for filenames: "EAE" -> "EAE", "NO EAE" -> "NOEAE"
sanitise_eae <- function(eae_level) gsub(" ", "", eae_level)

# ------------------------------------------------------------
# Load data (annotated Seurat object + marker voting)
# ------------------------------------------------------------

source_methods("load_annotated_object.R")
cat("Loading annotated Seurat object...\n")
so <- load_annotated_object()
cat(sprintf("  %d cells x %d features\n", ncol(so), nrow(so)))
cat(sprintf("  stable_20_clusters: %d levels; celltype_markers: %d levels\n",
    nlevels(factor(so$stable_20_clusters)), nlevels(factor(so$celltype_markers))))
gc()

# ------------------------------------------------------------
# Helper: run FindMarkers for a subset, return results + metadata
# ------------------------------------------------------------

run_deg <- function(so_sub, group_col, group_name, comparison_label) {
    Idents(so_sub) <- "condition"

    n_ctrl <- sum(so_sub$condition == "ctrl")
    n_ko   <- sum(so_sub$condition == "ko")

    cat(sprintf("  %s: ctrl=%d, ko=%d", comparison_label, n_ctrl, n_ko))

    if (n_ctrl == 0 || n_ko == 0) {
        cat(" -> SKIPPED (0 in one group)\n")
        return(NULL)
    }

    low_n <- min(n_ctrl, n_ko) < MIN_CELLS_WARNING
    if (low_n) cat(sprintf(" [WARNING: low n, min=%d]", min(n_ctrl, n_ko)))
    cat("\n")

    tryCatch({
        degs <- FindMarkers(
            so_sub,
            ident.1 = "ctrl",
            ident.2 = "ko",
            test.use = TEST_USE,
            logfc.threshold = LOGFC_THRESHOLD,
            min.pct = MIN_PCT,
            verbose = FALSE
        )

        if (nrow(degs) == 0) {
            cat("    -> 0 DEGs found\n")
            return(data.frame(
                comparison = comparison_label,
                group = group_name,
                n_ctrl = n_ctrl,
                n_ko = n_ko,
                n_deg_total = 0,
                n_deg_up_ctrl = 0,
                n_deg_up_ko = 0,
                n_deg_sig = 0,
                n_deg_sig_up_ctrl = 0,
                n_deg_sig_up_ko = 0,
                low_n = low_n,
                stringsAsFactors = FALSE
            ))
        }

        # Add gene column and metadata
        degs$gene <- rownames(degs)
        degs$comparison <- comparison_label
        degs$group <- group_name
        degs$n_ctrl <- n_ctrl
        degs$n_ko <- n_ko
        degs$low_n_warning <- low_n

        # Save individual CSV
        fname <- sprintf("DEG_%s.csv", gsub("[^A-Za-z0-9_]", "_", comparison_label))
        fpath <- file.path(OUTPUT_DIR, fname)
        write.csv(degs, fpath, row.names = FALSE)
        cat(sprintf("    -> %d DEGs, saved to %s\n", nrow(degs), basename(fpath)))

        # Summary row
        sig <- degs[degs$p_val_adj < 0.05, ]
        data.frame(
            comparison = comparison_label,
            group = group_name,
            n_ctrl = n_ctrl,
            n_ko = n_ko,
            n_deg_total = nrow(degs),
            n_deg_up_ctrl = sum(degs$avg_log2FC > 0),
            n_deg_up_ko = sum(degs$avg_log2FC < 0),
            n_deg_sig = nrow(sig),
            n_deg_sig_up_ctrl = sum(sig$avg_log2FC > 0),
            n_deg_sig_up_ko = sum(sig$avg_log2FC < 0),
            low_n = low_n,
            stringsAsFactors = FALSE
        )
    }, error = function(e) {
        cat(sprintf("    -> ERROR: %s\n", conditionMessage(e)))
        data.frame(
            comparison = comparison_label,
            group = group_name,
            n_ctrl = n_ctrl,
            n_ko = n_ko,
            n_deg_total = NA,
            n_deg_up_ctrl = NA,
            n_deg_up_ko = NA,
            n_deg_sig = NA,
            n_deg_sig_up_ctrl = NA,
            n_deg_sig_up_ko = NA,
            low_n = low_n,
            stringsAsFactors = FALSE
        )
    })
}

# ============================================================
# Main loop: iterate over EAE levels
# ============================================================

summary_rows <- list()

for (eae_level in EAE_LEVELS) {
    eae_tag <- sanitise_eae(eae_level)

    cat(sprintf("\n##########################################################\n"))
    cat(sprintf("# EAE status: %s\n", eae_level))
    cat(sprintf("##########################################################\n"))

    so_sub_eae <- subset(so, eae == eae_level)
    cat(sprintf("  Subset: %d cells\n", ncol(so_sub_eae)))

    # ----------------------------------------------------------
    # Cluster-level DEGs
    # ----------------------------------------------------------
    cat(sprintf("\n  ---- CLUSTER-LEVEL DEGs (stable_20_clusters) [%s] ----\n\n", eae_level))

    for (cl in CLUSTERS_OF_INTEREST) {
        n_cells <- sum(so_sub_eae$stable_20_clusters == cl, na.rm = TRUE)
        if (n_cells == 0) {
            cat(sprintf("  %s_cluster_%s_ctrl_vs_ko: 0 cells total -> SKIPPED\n", eae_tag, cl))
            next
        }
        so_cl <- subset(so_sub_eae, stable_20_clusters == cl)
        label <- sprintf("%s_cluster_%s_ctrl_vs_ko", eae_tag, cl)
        res <- run_deg(so_cl, "stable_20_clusters", cl, label)
        if (!is.null(res)) {
            res$eae_status <- eae_level
            res$comparison_type <- "cluster"
            summary_rows[[length(summary_rows) + 1]] <- res
        }
        rm(so_cl)
    }

    # ----------------------------------------------------------
    # Cell-type-level DEGs
    # ----------------------------------------------------------
    cat(sprintf("\n  ---- CELLTYPE-LEVEL DEGs (celltype_markers) [%s] ----\n\n", eae_level))

    all_celltypes <- levels(so_sub_eae$celltype_markers)
    celltypes_to_run <- setdiff(all_celltypes, CELLTYPES_SKIP)

    for (ct in celltypes_to_run) {
        n_cells <- sum(so_sub_eae$celltype_markers == ct, na.rm = TRUE)
        if (n_cells == 0) {
            cat(sprintf("  %s_celltype_%s_ctrl_vs_ko: 0 cells total -> SKIPPED\n", eae_tag, ct))
            next
        }
        so_ct <- subset(so_sub_eae, celltype_markers == ct)
        label <- sprintf("%s_celltype_%s_ctrl_vs_ko", eae_tag, ct)
        res <- run_deg(so_ct, "celltype_markers", ct, label)
        if (!is.null(res)) {
            res$eae_status <- eae_level
            res$comparison_type <- "celltype"
            summary_rows[[length(summary_rows) + 1]] <- res
        }
        rm(so_ct)
    }

    rm(so_sub_eae)
    gc()
}

# Free the full object
rm(so)
gc()

# ============================================================
# Summary table
# ============================================================

summary_df <- bind_rows(summary_rows)
summary_path <- file.path(OUTPUT_DIR, "DEG_summary.csv")
write.csv(summary_df, summary_path, row.names = FALSE)

cat("\n========================================\n")
cat("SUMMARY\n")
cat("========================================\n\n")
print(summary_df, right = FALSE)
cat(sprintf("\nSummary saved to %s\n", summary_path))
cat(sprintf("Individual DEG CSVs in %s\n", OUTPUT_DIR))

# ============================================================
# Generate Markdown report
# ============================================================

cat("\nGenerating Markdown report...\n")

report_path <- file.path(OUTPUT_DIR, "DEG_report.md")
report <- file(report_path, open = "w")

writeLines(c(
    "# DEG Analysis Report: CTRL vs KO (by EAE status)",
    "",
    sprintf("**Date**: %s", Sys.Date()),
    "",
    "**Dataset**: `seurat_object.rds`",
    "",
    "**Metadata source**: ClustAssess stable partition + marker-gene voting",
    "",
    "## Methods",
    "",
    sprintf("- **Test**: Wilcoxon rank-sum (`test.use = \"%s\"`)", TEST_USE),
    sprintf("- **Log2FC threshold**: %s (genes below this are not tested)", LOGFC_THRESHOLD),
    sprintf("- **Min percent expressed**: %s (gene must be detected in >= %s%% of cells in either group)", MIN_PCT, MIN_PCT * 100),
    "- **Comparison direction**: `ident.1 = ctrl`, `ident.2 = ko`",
    "  - Positive `avg_log2FC` = upregulated in **ctrl** (wildtype)",
    "  - Negative `avg_log2FC` = upregulated in **ko** (Sucnr1 knockout)",
    "- **Significance**: adjusted p-value < 0.05 (Bonferroni correction, Seurat default)",
    "",
    "- **Skipped cell types**: `Unassigned` (not a meaningful cell type). Any cluster/celltype with 0 cells in either ctrl or ko is also skipped.",
    ""
), report)

# Loop over EAE levels for report sections
for (eae_level in EAE_LEVELS) {
    eae_tag <- sanitise_eae(eae_level)
    eae_rows <- summary_df[summary_df$eae_status == eae_level, ]

    writeLines(c(
        "---",
        "",
        sprintf("# %s: CTRL vs KO", eae_level),
        ""
    ), report)

    # --- Cluster-level ---
    cluster_rows <- eae_rows[eae_rows$comparison_type == "cluster", ]
    writeLines(c(
        sprintf("## Cluster-level DEGs (`stable_20_clusters`) [%s]", eae_level),
        "",
        "Clustering derived from ClustAssess (Most_Abundant / 1950 features / SLM / 20 clusters).",
        ""
    ), report)

    if (nrow(cluster_rows) > 0) {
        writeLines(c(
            "| Cluster | n ctrl | n ko | DEGs total | Up in ctrl | Up in ko | Sig (adj.p<0.05) | Sig up ctrl | Sig up ko | Low n |",
            "|---------|--------|------|------------|------------|----------|------------------|-------------|-----------|-------|"
        ), report)
        for (i in seq_len(nrow(cluster_rows))) {
            r <- cluster_rows[i, ]
            writeLines(sprintf("| %s | %d | %d | %s | %s | %s | %s | %s | %s | %s |",
                r$group, r$n_ctrl, r$n_ko,
                ifelse(is.na(r$n_deg_total), "ERROR", as.character(r$n_deg_total)),
                ifelse(is.na(r$n_deg_up_ctrl), "-", as.character(r$n_deg_up_ctrl)),
                ifelse(is.na(r$n_deg_up_ko), "-", as.character(r$n_deg_up_ko)),
                ifelse(is.na(r$n_deg_sig), "-", as.character(r$n_deg_sig)),
                ifelse(is.na(r$n_deg_sig_up_ctrl), "-", as.character(r$n_deg_sig_up_ctrl)),
                ifelse(is.na(r$n_deg_sig_up_ko), "-", as.character(r$n_deg_sig_up_ko)),
                ifelse(r$low_n, "YES", "no")
            ), report)
        }
        writeLines("", report)

        # Top genes per cluster
        for (i in seq_len(nrow(cluster_rows))) {
            r <- cluster_rows[i, ]
            if (!is.na(r$n_deg_sig) && r$n_deg_sig > 0) {
                fname <- sprintf("DEG_%s_cluster_%s_ctrl_vs_ko.csv", eae_tag, r$group)
                fpath <- file.path(OUTPUT_DIR, fname)
                if (file.exists(fpath)) {
                    deg_data <- read.csv(fpath)
                    deg_sig <- deg_data[deg_data$p_val_adj < 0.05, ]
                    deg_sig <- deg_sig[order(deg_sig$p_val_adj), ]

                    writeLines(c(
                        sprintf("### Cluster %s: top significant DEGs", r$group),
                        ""
                    ), report)

                    n_show <- min(20, nrow(deg_sig))
                    writeLines(c(
                        "| Gene | avg_log2FC | pct.1 (ctrl) | pct.2 (ko) | p_val_adj | Direction |",
                        "|------|------------|--------------|------------|-----------|-----------|"
                    ), report)
                    for (j in seq_len(n_show)) {
                        d <- deg_sig[j, ]
                        direction <- ifelse(d$avg_log2FC > 0, "up in ctrl", "up in ko")
                        writeLines(sprintf("| %s | %.3f | %.3f | %.3f | %.2e | %s |",
                            d$gene, d$avg_log2FC, d$pct.1, d$pct.2, d$p_val_adj, direction
                        ), report)
                    }
                    if (nrow(deg_sig) > n_show) {
                        writeLines(sprintf("\n*... and %d more significant DEGs (see CSV)*\n",
                            nrow(deg_sig) - n_show), report)
                    }
                    writeLines("", report)
                }
            }
        }
    }

    # --- Cell-type-level ---
    ct_rows <- eae_rows[eae_rows$comparison_type == "celltype", ]
    writeLines(c(
        sprintf("## Cell-type-level DEGs (`celltype_markers`) [%s]", eae_level),
        "",
        "Cell types assigned by marker gene voting (`marker_genes.R`).",
        ""
    ), report)

    if (nrow(ct_rows) > 0) {
        writeLines(c(
            "| Cell type | n ctrl | n ko | DEGs total | Up in ctrl | Up in ko | Sig (adj.p<0.05) | Sig up ctrl | Sig up ko | Low n |",
            "|-----------|--------|------|------------|------------|----------|------------------|-------------|-----------|-------|"
        ), report)
        for (i in seq_len(nrow(ct_rows))) {
            r <- ct_rows[i, ]
            writeLines(sprintf("| %s | %d | %d | %s | %s | %s | %s | %s | %s | %s |",
                r$group, r$n_ctrl, r$n_ko,
                ifelse(is.na(r$n_deg_total), "ERROR", as.character(r$n_deg_total)),
                ifelse(is.na(r$n_deg_up_ctrl), "-", as.character(r$n_deg_up_ctrl)),
                ifelse(is.na(r$n_deg_up_ko), "-", as.character(r$n_deg_up_ko)),
                ifelse(is.na(r$n_deg_sig), "-", as.character(r$n_deg_sig)),
                ifelse(is.na(r$n_deg_sig_up_ctrl), "-", as.character(r$n_deg_sig_up_ctrl)),
                ifelse(is.na(r$n_deg_sig_up_ko), "-", as.character(r$n_deg_sig_up_ko)),
                ifelse(r$low_n, "YES", "no")
            ), report)
        }
        writeLines("", report)

        # Top genes per cell type
        for (i in seq_len(nrow(ct_rows))) {
            r <- ct_rows[i, ]
            if (!is.na(r$n_deg_sig) && r$n_deg_sig > 0) {
                fname <- sprintf("DEG_%s_celltype_%s_ctrl_vs_ko.csv", eae_tag, r$group)
                fpath <- file.path(OUTPUT_DIR, fname)
                if (file.exists(fpath)) {
                    deg_data <- read.csv(fpath)
                    deg_sig <- deg_data[deg_data$p_val_adj < 0.05, ]
                    deg_sig <- deg_sig[order(deg_sig$p_val_adj), ]

                    writeLines(c(
                        sprintf("### %s: top significant DEGs", r$group),
                        ""
                    ), report)

                    n_show <- min(20, nrow(deg_sig))
                    writeLines(c(
                        "| Gene | avg_log2FC | pct.1 (ctrl) | pct.2 (ko) | p_val_adj | Direction |",
                        "|------|------------|--------------|------------|-----------|-----------|"
                    ), report)
                    for (j in seq_len(n_show)) {
                        d <- deg_sig[j, ]
                        direction <- ifelse(d$avg_log2FC > 0, "up in ctrl", "up in ko")
                        writeLines(sprintf("| %s | %.3f | %.3f | %.3f | %.2e | %s |",
                            d$gene, d$avg_log2FC, d$pct.1, d$pct.2, d$p_val_adj, direction
                        ), report)
                    }
                    if (nrow(deg_sig) > n_show) {
                        writeLines(sprintf("\n*... and %d more significant DEGs (see CSV)*\n",
                            nrow(deg_sig) - n_show), report)
                    }
                    writeLines("", report)
                }
            }
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
    "| `DEG_summary.csv` | Summary table (one row per comparison) |",
    "| `DEG_<EAE/NOEAE>_cluster_<N>_ctrl_vs_ko.csv` | Per-gene DEG results for cluster N |",
    "| `DEG_<EAE/NOEAE>_celltype_<name>_ctrl_vs_ko.csv` | Per-gene DEG results for cell type |",
    "| `DEG_report.md` | This report |",
    ""
), report)

close(report)
cat(sprintf("Report saved to %s\n", report_path))



cat("\nDone.\n")
