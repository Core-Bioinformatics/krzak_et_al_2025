library(dplyr)
library(Seurat)

options(stringsAsFactors = FALSE)

# ============================================================
# Dedicated pooled-cell Fisher exact analysis
#
# This script reproduces the project Fisher abundance-testing
# style using row-normalized 2x2 tables.
#
# Each condition row is rescaled to sum to 100 before the
# Fisher test is run. This keeps the test focused on relative
# composition differences rather than raw cell-yield scale,
# while avoiding the implicit rounding that occurred when
# non-integer percentages were passed directly to fisher.test().
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

OUTPUT_DIR <- ensure_out_dir("cell_abundance_fisher")
PROCESSED_DIR <- OUTPUT_DIR

source_methods("load_annotated_object.R")

sort_levels_mixed <- function(x) {
    x <- unique(as.character(x))
    x <- x[!is.na(x)]
    nums <- suppressWarnings(as.numeric(x))
    if (!any(is.na(nums))) {
        return(as.character(sort(nums)))
    }
    sort(x)
}

prepare_meta_subset <- function(meta, category_col, eae_value, exclude_values = NULL) {
    keep <- !is.na(meta[[category_col]]) & meta$eae == eae_value & meta$condition %in% c("ctrl", "ko")
    if (!is.null(exclude_values)) {
        keep <- keep & !(meta[[category_col]] %in% exclude_values)
    }
    meta[keep, , drop = FALSE]
}

scale_row_to_total_100 <- function(in_count, out_count) {
    total <- in_count + out_count
    if (total == 0) {
        return(c(in_scaled = NA_integer_, out_scaled = NA_integer_))
    }
    in_scaled <- as.integer(round(100 * in_count / total))
    in_scaled <- max(0L, min(100L, in_scaled))
    out_scaled <- 100L - in_scaled
    c(in_scaled = in_scaled, out_scaled = out_scaled)
}

run_scaled_fisher <- function(meta, category_col, eae_value, exclude_values = NULL) {
    df <- prepare_meta_subset(meta, category_col, eae_value, exclude_values)
    categories <- sort_levels_mixed(df[[category_col]])

    res <- lapply(categories, function(category_id) {
        ctrl_in <- sum(df$condition == "ctrl" & df[[category_col]] == category_id)
        ctrl_out <- sum(df$condition == "ctrl" & df[[category_col]] != category_id)
        ko_in <- sum(df$condition == "ko" & df[[category_col]] == category_id)
        ko_out <- sum(df$condition == "ko" & df[[category_col]] != category_id)

        ctrl_scaled <- scale_row_to_total_100(ctrl_in, ctrl_out)
        ko_scaled <- scale_row_to_total_100(ko_in, ko_out)

        tab <- matrix(
            c(ctrl_scaled[["in_scaled"]], ctrl_scaled[["out_scaled"]],
              ko_scaled[["in_scaled"]], ko_scaled[["out_scaled"]]),
            nrow = 2,
            byrow = TRUE,
            dimnames = list(condition = c("ctrl", "ko"), membership = c("in_category", "out_category"))
        )

        ft <- fisher.test(tab, alternative = "two.sided", conf.int = TRUE, conf.level = 0.95)

        ctrl_total <- ctrl_in + ctrl_out
        ko_total <- ko_in + ko_out
        ctrl_prop <- if (ctrl_total > 0) ctrl_in / ctrl_total else NA_real_
        ko_prop <- if (ko_total > 0) ko_in / ko_total else NA_real_

        direction <- if (isTRUE(all.equal(ctrl_prop, ko_prop))) {
            "no_change"
        } else if (ko_prop > ctrl_prop) {
            "enriched_in_ko"
        } else {
            "enriched_in_ctrl"
        }

        data.frame(
            category = category_id,
            ctrl_cells = ctrl_in,
            ko_cells = ko_in,
            ctrl_other_cells = ctrl_out,
            ko_other_cells = ko_out,
            ctrl_total_cells = ctrl_total,
            ko_total_cells = ko_total,
            ctrl_prop = ctrl_prop,
            ko_prop = ko_prop,
            ctrl_in_scaled100 = ctrl_scaled[["in_scaled"]],
            ctrl_out_scaled100 = ctrl_scaled[["out_scaled"]],
            ko_in_scaled100 = ko_scaled[["in_scaled"]],
            ko_out_scaled100 = ko_scaled[["out_scaled"]],
            odds_ratio = unname(ft$estimate),
            conf_low = unname(ft$conf.int[1]),
            conf_high = unname(ft$conf.int[2]),
            p_value = ft$p.value,
            direction = direction,
            stringsAsFactors = FALSE
        )
    })

    res <- bind_rows(res)
    res$padj <- p.adjust(res$p_value, method = "BH")
    res <- res %>%
        arrange(.data$padj, .data$p_value, .data$category)
    res
}

write_output <- function(df, filename) {
    output_path <- file.path(OUTPUT_DIR, filename)
    processed_path <- file.path(PROCESSED_DIR, filename)
    write.csv(df, output_path, row.names = FALSE)
    file.copy(output_path, processed_path, overwrite = TRUE)
    cat(sprintf("Saved %s\n", output_path))
}

cat("Preparing annotated object...\n")
so <- load_annotated_object()
meta <- so[[]]

cluster_eae <- run_scaled_fisher(meta, "stable_20_clusters", "EAE")
cluster_noeae <- run_scaled_fisher(meta, "stable_20_clusters", "NO EAE")
celltype_eae <- run_scaled_fisher(meta, "celltype_markers", "EAE", exclude_values = "Unassigned")
celltype_noeae <- run_scaled_fisher(meta, "celltype_markers", "NO EAE", exclude_values = "Unassigned")

write_output(cluster_eae, "cluster_fisher-stable20_vs_condition-EAE.csv")
write_output(cluster_noeae, "cluster_fisher-stable20_vs_condition-NO-EAE.csv")
write_output(celltype_eae, "celltype_fisher-celltype_vs_condition-EAE.csv")
write_output(celltype_noeae, "celltype_fisher-celltype_vs_condition-NO-EAE.csv")

combined <- bind_rows(
    mutate(cluster_eae, eae = "EAE", category_col = "stable_20_clusters"),
    mutate(cluster_noeae, eae = "NO EAE", category_col = "stable_20_clusters"),
    mutate(celltype_eae, eae = "EAE", category_col = "celltype_markers"),
    mutate(celltype_noeae, eae = "NO EAE", category_col = "celltype_markers")
) %>%
    relocate("eae", "category_col")
write_output(combined, "cell_abundance_fisher_all_results.csv")

readme_path <- file.path(OUTPUT_DIR, "README.md")
writeLines(
    c(
        "# Cell Abundance Fisher Analysis",
        "",
        sprintf("**Date**: %s", Sys.Date()),
        "",
        "## What This Folder Contains",
        "",
        "This folder contains pooled-cell Fisher exact test results for cell abundance comparisons.",
        "",
        "- compare `ctrl` (WT) vs `ko` within `EAE` and within `NO EAE`",
        "- do this separately for `stable_20_clusters` and `celltype_markers`",
        "- report counts, odds ratios, raw p-values, and BH-adjusted p-values",
        "",
        "## Data Used",
        "",
        "- Base object: `seurat_object.rds`",
        "- Stable clusters: from `clustassess_object.rds` (Most Abundant / 1,950 / SLM / 20)",
        "- Marker-based cell types: assigned with `marker_genes.R`",
        "- Cell-type tables exclude `Unassigned`, matching the current cell-type plotting workflow",
        "",
        "## Statistical Method",
        "",
        "These files use pooled-cell Fisher exact tests on row-normalized 2x2 tables:",
        "",
        "- in category vs all other cells",
        "- `ctrl` vs `ko`",
        "- each condition row is rescaled to sum to 100 before testing",
        "- the scaled values are stored explicitly in the output as `*_scaled100`",
        "- two-sided Fisher exact test",
        "- BH correction within each result table",
        "",
        "## Which Files To Open First",
        "",
        "- `cluster_fisher-stable20_vs_condition-EAE.csv`: cluster-level Fisher results for EAE",
        "- `cluster_fisher-stable20_vs_condition-NO-EAE.csv`: cluster-level Fisher results for NO EAE",
        "- `celltype_fisher-celltype_vs_condition-EAE.csv`: cell-type Fisher results for EAE",
        "- `celltype_fisher-celltype_vs_condition-NO-EAE.csv`: cell-type Fisher results for NO EAE",
        "- `cell_abundance_fisher_all_results.csv`: combined table across all four analyses",
        "",
        "## Interpretation Note",
        "",
        "These are pooled-cell Fisher tests. They are useful for composition-style comparisons, but they treat cells as independent observations."
    ),
    con = readme_path
)
file.copy(readme_path, file.path(PROCESSED_DIR, "README.md"), overwrite = TRUE)

cat("Done.\n")
