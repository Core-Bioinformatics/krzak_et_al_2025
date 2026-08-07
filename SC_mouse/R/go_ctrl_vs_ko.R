library(gprofiler2)
library(dplyr)

# ============================================================
# GO / pathway enrichment on the CTRL vs KO DEGs
#
# Consumes the CSVs written by deg_ctrl_vs_ko.R:
#   DEG_<EAE|NOEAE>_<cluster|celltype>_<label>_ctrl_vs_ko.csv
#
# For each comparison, three queries are submitted:
#   up_ctrl  — significant genes with avg_log2FC > 0 (higher in ctrl/WT)
#   up_ko    — significant genes with avg_log2FC < 0 (higher in Sucnr1 KO)
#   all      — union of the two (direction-agnostic)
#
# Enrichment uses g:Profiler (gprofiler2::gost), g:SCS multiple-testing
# correction, against a custom background of the genes present in the
# analysed object (protein-coding, haemoglobin-depleted).
#
# Output: one CSV per query + a summary CSV + Markdown report.
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

DEG_DIR <- file.path(OUT_DIR, "deg_ctrl_vs_ko")
OUTPUT_DIR <- ensure_out_dir("go_ctrl_vs_ko")

# Background (gene universe) for the enrichment test.
#   "custom"    — genes retained in the analysed object (recommended)
#   "annotated" — g:Profiler default (all annotated genes for the organism)
BACKGROUND_MODE <- "custom"
BACKGROUND_PATH <- file.path(DATA_DIR, "gene_background.csv")

ORGANISM <- "mmusculus"

# Annotation sources to query. NULL = every source g:Profiler offers.
GO_SOURCES <- c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC")

# DEG selection thresholds (applied on top of the FindMarkers filtering)
FDR_CUTOFF <- 0.05        # p_val_adj (Bonferroni, as returned by Seurat)
LFC_CUTOFF <- 0.25        # |avg_log2FC|

# Skip queries with fewer than this many genes (enrichment is meaningless)
MIN_GENES <- 5

# Directions to test per comparison
DIRECTIONS <- c("up_ctrl", "up_ko", "all")

# Network retries for the g:Profiler API
MAX_RETRIES <- 3
RETRY_SLEEP <- 10         # seconds

# Terms listed per query in the Markdown report
REPORT_TOP_N <- 10


# ------------------------------------------------------------
# Background gene universe
# ------------------------------------------------------------

background <- NULL
domain_scope <- "annotated"

if (BACKGROUND_MODE == "custom") {
    if (!file.exists(BACKGROUND_PATH)) {
        stop(sprintf("Background gene list not found: %s", BACKGROUND_PATH))
    }
    bg <- read.csv(BACKGROUND_PATH, stringsAsFactors = FALSE)
    background <- unique(bg[[1]])
    background <- background[!is.na(background) & nzchar(background)]
    domain_scope <- "custom"
    cat(sprintf("Background: %d genes from %s\n", length(background),
        basename(BACKGROUND_PATH)))
} else {
    cat("Background: g:Profiler annotated domain (no custom background)\n")
}

# ------------------------------------------------------------
# Input DEG tables
# ------------------------------------------------------------

deg_files <- list.files(DEG_DIR,
    pattern = "^DEG_(EAE|NOEAE)_(cluster|celltype)_.*_ctrl_vs_ko\\.csv$",
    full.names = TRUE)

if (length(deg_files) == 0) {
    stop(sprintf("No DEG CSVs found in %s — run deg_ctrl_vs_ko.R first.", DEG_DIR))
}

cat(sprintf("Found %d DEG tables in %s\n", length(deg_files), DEG_DIR))

# Parse "DEG_EAE_cluster_1_ctrl_vs_ko.csv" -> stratum / type / label
parse_deg_filename <- function(path) {
    stem <- sub("\\.csv$", "", basename(path))
    stem <- sub("^DEG_", "", stem)
    stem <- sub("_ctrl_vs_ko$", "", stem)
    parts <- strsplit(stem, "_", fixed = TRUE)[[1]]
    list(
        eae_status = ifelse(parts[1] == "NOEAE", "NO EAE", "EAE"),
        comparison_type = parts[2],
        group = paste(parts[-(1:2)], collapse = "_")
    )
}

# gost() returns list-columns (parents, evidence_codes, intersection);
# collapse them so the table can be written as flat CSV.
flatten_list_cols <- function(df) {
    for (nm in names(df)) {
        if (is.list(df[[nm]])) {
            df[[nm]] <- vapply(df[[nm]],
                function(x) paste(unlist(x), collapse = "|"),
                character(1))
        }
    }
    df
}

# gost() returns NULL both when nothing is enriched and (via error) when the
# API call fails, so wrap it to tell the two apart:
#   list(ok = TRUE,  gost = <gost() return, may be NULL>)  — query ran
#   list(ok = FALSE, gost = NULL)                          — query never succeeded
run_gost <- function(genes, tag) {
    for (attempt in seq_len(MAX_RETRIES)) {
        res <- tryCatch(
            list(ok = TRUE, gost = gost(
                query = genes,
                organism = ORGANISM,
                ordered_query = FALSE,
                significant = TRUE,
                user_threshold = 0.05,
                correction_method = "g_SCS",
                domain_scope = domain_scope,
                custom_bg = background,
                sources = GO_SOURCES,
                evcodes = TRUE
            )),
            error = function(e) {
                cat(sprintf("    attempt %d failed: %s\n", attempt, conditionMessage(e)))
                NULL
            }
        )
        if (!is.null(res)) return(res)
        if (attempt < MAX_RETRIES) Sys.sleep(RETRY_SLEEP)
    }
    cat(sprintf("    WARNING: g:Profiler query failed after %d attempts (%s)\n",
        MAX_RETRIES, tag))
    list(ok = FALSE, gost = NULL)
}

# ------------------------------------------------------------
# Main loop
# ------------------------------------------------------------

summary_rows <- list()
term_store <- list()   # keyed by output stem, used by the report

for (deg_file in sort(deg_files)) {
    meta <- parse_deg_filename(deg_file)
    deg <- read.csv(deg_file, stringsAsFactors = FALSE)

    if (nrow(deg) == 0 || !("gene" %in% names(deg))) {
        cat(sprintf("\n%s: empty or malformed — skipped\n", basename(deg_file)))
        next
    }

    sig <- deg[!is.na(deg$p_val_adj) & deg$p_val_adj < FDR_CUTOFF &
               abs(deg$avg_log2FC) >= LFC_CUTOFF, ]

    label <- sprintf("%s / %s %s", meta$eae_status, meta$comparison_type, meta$group)
    cat(sprintf("\n%s: %d significant DEGs (of %d tested)\n",
        label, nrow(sig), nrow(deg)))

    for (direction in DIRECTIONS) {
        genes <- switch(direction,
            up_ctrl = sig$gene[sig$avg_log2FC > 0],
            up_ko   = sig$gene[sig$avg_log2FC < 0],
            all     = sig$gene
        )
        genes <- unique(genes[!is.na(genes) & nzchar(genes)])

        eae_tag <- gsub(" ", "", meta$eae_status)
        out_stem <- sprintf("GO_%s_%s_%s_%s",
            eae_tag, meta$comparison_type, meta$group, direction)

        row <- data.frame(
            eae_status = meta$eae_status,
            comparison_type = meta$comparison_type,
            group = meta$group,
            direction = direction,
            n_genes_query = length(genes),
            n_terms = NA_integer_,
            n_terms_GOBP = NA_integer_,
            top_term_id = NA_character_,
            top_term_name = NA_character_,
            top_term_p = NA_real_,
            status = NA_character_,
            output_file = NA_character_,
            stringsAsFactors = FALSE
        )

        if (length(genes) < MIN_GENES) {
            cat(sprintf("  %-8s %3d genes — skipped (< %d)\n",
                direction, length(genes), MIN_GENES))
            row$status <- "skipped_too_few_genes"
            summary_rows[[length(summary_rows) + 1]] <- row
            next
        }

        cat(sprintf("  %-8s %3d genes — querying g:Profiler...\n",
            direction, length(genes)))
        res <- run_gost(genes, out_stem)

        if (!res$ok) {
            row$status <- "query_failed"
            summary_rows[[length(summary_rows) + 1]] <- row
            next
        }

        if (is.null(res$gost) || is.null(res$gost$result) || nrow(res$gost$result) == 0) {
            cat("    no significant terms\n")
            row$n_terms <- 0L
            row$status <- "no_significant_terms"
            summary_rows[[length(summary_rows) + 1]] <- row
            next
        }

        terms <- res$gost$result
        terms <- terms[order(terms$p_value), ]
        terms$eae_status <- meta$eae_status
        terms$comparison_type <- meta$comparison_type
        terms$group <- meta$group
        terms$direction <- direction

        out_path <- file.path(OUTPUT_DIR, paste0(out_stem, ".csv"))
        write.csv(flatten_list_cols(terms), out_path, row.names = FALSE)

        row$n_terms <- nrow(terms)
        row$n_terms_GOBP <- sum(terms$source == "GO:BP")
        row$top_term_id <- terms$term_id[1]
        row$top_term_name <- terms$term_name[1]
        row$top_term_p <- terms$p_value[1]
        row$status <- "ok"
        row$output_file <- basename(out_path)
        summary_rows[[length(summary_rows) + 1]] <- row

        term_store[[out_stem]] <- terms[, c("source", "term_id", "term_name",
            "p_value", "term_size", "intersection_size")]

        cat(sprintf("    %d terms (%d GO:BP) -> %s\n",
            nrow(terms), row$n_terms_GOBP, basename(out_path)))
    }
}

# ------------------------------------------------------------
# Summary table
# ------------------------------------------------------------

summary_df <- bind_rows(summary_rows)
summary_path <- file.path(OUTPUT_DIR, "GO_summary.csv")
write.csv(summary_df, summary_path, row.names = FALSE)
cat(sprintf("\nSummary saved to %s (%d rows)\n", summary_path, nrow(summary_df)))

# ------------------------------------------------------------
# Markdown report
# ------------------------------------------------------------

report_path <- file.path(OUTPUT_DIR, "GO_report.md")
report <- file(report_path, "w")

writeLines(c(
    "# GO / pathway enrichment — CTRL vs KO DEGs",
    "",
    sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    "## Method",
    "",
    sprintf("- Input: `%s` (output of `deg_ctrl_vs_ko.R`)", DEG_DIR),
    sprintf("- DEG selection: adjusted *P* < %.2g and |log2FC| >= %.2f",
        FDR_CUTOFF, LFC_CUTOFF),
    "- Directions: `up_ctrl` (higher in ctrl/WT), `up_ko` (higher in Sucnr1 KO), `all` (both)",
    sprintf("- Enrichment: g:Profiler (`gprofiler2::gost`), organism `%s`, g:SCS correction, threshold 0.05",
        ORGANISM),
    sprintf("- Sources: %s", paste(GO_SOURCES, collapse = ", ")),
    sprintf("- Background: %s",
        if (domain_scope == "custom")
            sprintf("custom, %d genes (`%s`)", length(background), basename(BACKGROUND_PATH))
        else "g:Profiler annotated domain"),
    sprintf("- Queries with fewer than %d genes were skipped", MIN_GENES),
    "",
    "## Summary",
    "",
    sprintf("- Queries attempted: %d", nrow(summary_df)),
    sprintf("- Queries with enriched terms: %d", sum(summary_df$status == "ok")),
    sprintf("- Queries skipped (too few genes): %d",
        sum(summary_df$status == "skipped_too_few_genes")),
    sprintf("- Queries with no significant terms: %d",
        sum(summary_df$status == "no_significant_terms")),
    sprintf("- Queries failed: %d", sum(summary_df$status == "query_failed")),
    ""
), report)

for (eae_status in unique(summary_df$eae_status)) {
    writeLines(c(sprintf("## %s", eae_status), ""), report)

    for (ctype in c("cluster", "celltype")) {
        sub <- summary_df[summary_df$eae_status == eae_status &
                          summary_df$comparison_type == ctype, ]
        if (nrow(sub) == 0) next

        writeLines(c(sprintf("### %s-level", ctype), ""), report)

        for (grp in unique(sub$group)) {
            grp_rows <- sub[sub$group == grp, ]
            writeLines(c(sprintf("#### %s", grp), ""), report)

            for (i in seq_len(nrow(grp_rows))) {
                r <- grp_rows[i, ]
                writeLines(sprintf("**%s** — %d genes, %s", r$direction,
                    r$n_genes_query,
                    if (r$status == "ok") sprintf("%d enriched terms", r$n_terms)
                    else r$status), report)
                writeLines("", report)

                stem <- sprintf("GO_%s_%s_%s_%s", gsub(" ", "", r$eae_status),
                    r$comparison_type, r$group, r$direction)
                if (!is.null(term_store[[stem]])) {
                    tt <- head(term_store[[stem]], REPORT_TOP_N)
                    writeLines(c(
                        "| Source | Term ID | Term | Adj. P | Term size | Overlap |",
                        "|--------|---------|------|--------|-----------|---------|"
                    ), report)
                    for (j in seq_len(nrow(tt))) {
                        writeLines(sprintf("| %s | %s | %s | %.2e | %d | %d |",
                            tt$source[j], tt$term_id[j], tt$term_name[j],
                            tt$p_value[j], tt$term_size[j],
                            tt$intersection_size[j]), report)
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
    "| `GO_summary.csv` | Summary table (one row per comparison x direction) |",
    "| `GO_<EAE/NOEAE>_cluster_<N>_<direction>.csv` | Enriched terms for cluster N |",
    "| `GO_<EAE/NOEAE>_celltype_<name>_<direction>.csv` | Enriched terms for cell type |",
    "| `GO_report.md` | This report |",
    ""
), report)

close(report)
cat(sprintf("Report saved to %s\n", report_path))



cat("\nDone.\n")
