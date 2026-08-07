library(dplyr)
library(Seurat)

options(stringsAsFactors = FALSE)

# ============================================================
# Cell abundance statistics for WT(ctrl) vs KO
#
# Primary analysis:
#   Sample-level composition test with biological replication
#   using a quasibinomial GLM on per-sample in-category counts.
#
# Secondary analysis:
#   Legacy pooled-cell Fisher exact test, matching the March
#   notebook intent but using integer counts rather than rounded
#   percentages.
#
# Output:
#   CSV tables under output/
#   Markdown methods summary
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

OUTPUT_DIR <- ensure_out_dir("cell_abundance_stats")
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

filter_category_meta <- function(meta, category_col, eae_value, exclude_values = NULL) {
    keep <- !is.na(meta[[category_col]]) & meta$eae == eae_value
    if (!is.null(exclude_values)) {
        keep <- keep & !(meta[[category_col]] %in% exclude_values)
    }
    meta[keep, , drop = FALSE]
}

run_sample_level_composition <- function(
    meta,
    category_col,
    eae_value,
    exclude_values = NULL,
    sample_col = "sample_name",
    group_col = "condition"
) {
    df <- filter_category_meta(meta, category_col, eae_value, exclude_values)
    df[[group_col]] <- factor(df[[group_col]], levels = c("ctrl", "ko"))

    totals <- df %>%
        count(.data[[sample_col]], .data[[group_col]], name = "sample_total_cells")

    categories <- sort_levels_mixed(df[[category_col]])
    results <- lapply(categories, function(category_id) {
        in_cat <- df %>%
            filter(.data[[category_col]] == category_id) %>%
            count(.data[[sample_col]], name = "in_category_cells")

        sample_df <- totals %>%
            left_join(in_cat, by = sample_col) %>%
            mutate(
                in_category_cells = ifelse(is.na(.data$in_category_cells), 0L, .data$in_category_cells),
                out_category_cells = .data$sample_total_cells - .data$in_category_cells,
                proportion = .data$in_category_cells / .data$sample_total_cells,
                condition = factor(.data[[group_col]], levels = c("ctrl", "ko"))
            ) %>%
            arrange(.data$condition, .data[[sample_col]])

        fit_warning <- NA_character_
        fit <- withCallingHandlers(
            glm(
                cbind(in_category_cells, out_category_cells) ~ condition,
                family = quasibinomial(),
                data = sample_df
            ),
            warning = function(w) {
                fit_warning <<- conditionMessage(w)
                invokeRestart("muffleWarning")
            }
        )

        coef_mat <- summary(fit)$coefficients
        coef_name <- "conditionko"
        if (!coef_name %in% rownames(coef_mat)) {
            estimate <- NA_real_
            std_error <- NA_real_
            statistic <- NA_real_
            p_value <- NA_real_
        } else {
            estimate <- unname(coef_mat[coef_name, "Estimate"])
            std_error <- unname(coef_mat[coef_name, "Std. Error"])
            statistic <- unname(coef_mat[coef_name, "t value"])
            p_value <- unname(coef_mat[coef_name, "Pr(>|t|)"])
        }

        ctrl_mask <- sample_df$condition == "ctrl"
        ko_mask <- sample_df$condition == "ko"
        mean_prop_ctrl <- mean(sample_df$proportion[ctrl_mask])
        mean_prop_ko <- mean(sample_df$proportion[ko_mask])
        direction <- if (isTRUE(all.equal(mean_prop_ctrl, mean_prop_ko))) {
            "no_change"
        } else if (mean_prop_ko > mean_prop_ctrl) {
            "enriched_in_ko"
        } else {
            "enriched_in_ctrl"
        }

        data.frame(
            analysis = "sample_level_quasibinomial",
            eae = eae_value,
            category_col = category_col,
            category = category_id,
            ctrl_cells = sum(sample_df$in_category_cells[ctrl_mask]),
            ko_cells = sum(sample_df$in_category_cells[ko_mask]),
            ctrl_total_cells = sum(sample_df$sample_total_cells[ctrl_mask]),
            ko_total_cells = sum(sample_df$sample_total_cells[ko_mask]),
            n_ctrl_samples = sum(ctrl_mask),
            n_ko_samples = sum(ko_mask),
            mean_prop_ctrl = mean_prop_ctrl,
            mean_prop_ko = mean_prop_ko,
            log_odds_ko_vs_ctrl = estimate,
            std_error = std_error,
            test_statistic = statistic,
            p_value = p_value,
            direction = direction,
            warning = fit_warning,
            stringsAsFactors = FALSE
        )
    })

    res <- bind_rows(results)
    res$padj <- p.adjust(res$p_value, method = "BH")
    res <- res %>%
        arrange(.data$padj, .data$p_value, .data$category)
    res
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

run_pooled_fisher <- function(
    meta,
    category_col,
    eae_value,
    exclude_values = NULL,
    group_col = "condition"
) {
    df <- filter_category_meta(meta, category_col, eae_value, exclude_values)
    df[[group_col]] <- factor(df[[group_col]], levels = c("ctrl", "ko"))
    categories <- sort_levels_mixed(df[[category_col]])

    results <- lapply(categories, function(category_id) {
        a <- sum(df[[group_col]] == "ctrl" & df[[category_col]] == category_id)
        b <- sum(df[[group_col]] == "ctrl" & df[[category_col]] != category_id)
        c <- sum(df[[group_col]] == "ko" & df[[category_col]] == category_id)
        d <- sum(df[[group_col]] == "ko" & df[[category_col]] != category_id)

        ctrl_scaled <- scale_row_to_total_100(a, b)
        ko_scaled <- scale_row_to_total_100(c, d)

        ft <- fisher.test(
            matrix(
                c(ctrl_scaled[["in_scaled"]], ctrl_scaled[["out_scaled"]],
                  ko_scaled[["in_scaled"]], ko_scaled[["out_scaled"]]),
                nrow = 2,
                byrow = TRUE,
                dimnames = list(condition = c("ctrl", "ko"), membership = c("in_category", "out_category"))
            )
        )

        ctrl_total <- a + b
        ko_total <- c + d
        ctrl_prop <- if (ctrl_total > 0) a / ctrl_total else NA_real_
        ko_prop <- if (ko_total > 0) c / ko_total else NA_real_
        direction <- if (isTRUE(all.equal(ctrl_prop, ko_prop))) {
            "no_change"
        } else if (ko_prop > ctrl_prop) {
            "enriched_in_ko"
        } else {
            "enriched_in_ctrl"
        }

        data.frame(
            analysis = "pooled_cell_fisher",
            eae = eae_value,
            category_col = category_col,
            category = category_id,
            ctrl_cells = a,
            ko_cells = c,
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

    res <- bind_rows(results)
    res$padj <- p.adjust(res$p_value, method = "BH")
    res <- res %>%
        arrange(.data$padj, .data$p_value, .data$category)
    res
}

write_output_table <- function(df, filename) {
    output_path <- file.path(OUTPUT_DIR, filename)
    processed_path <- file.path(PROCESSED_DIR, filename)
    write.csv(df, output_path, row.names = FALSE)
    file.copy(output_path, processed_path, overwrite = TRUE)
    cat(sprintf("Saved %s\n", output_path))
}

cat("Preparing annotated object...\n")
so <- load_annotated_object()
meta <- so[[]]

sample_summary <- meta %>%
    count(sample_name, eae, condition, rep, name = "n_cells") %>%
    arrange(.data$eae, .data$condition, .data$rep)
write_output_table(sample_summary, "sample_cell_counts.csv")

analysis_specs <- list(
    list(category_col = "stable_20_clusters", label = "stable20", exclude_values = NULL),
    list(category_col = "celltype_markers", label = "celltypes", exclude_values = "Unassigned")
)

eae_levels <- c("EAE", "NO EAE")
all_results <- list()

for (spec in analysis_specs) {
    for (eae_value in eae_levels) {
        sample_res <- run_sample_level_composition(
            meta = meta,
            category_col = spec$category_col,
            eae_value = eae_value,
            exclude_values = spec$exclude_values
        )
        fisher_res <- run_pooled_fisher(
            meta = meta,
            category_col = spec$category_col,
            eae_value = eae_value,
            exclude_values = spec$exclude_values
        )

        sample_file <- sprintf(
            "%s_sample_level_%s.csv",
            tolower(gsub(" ", "", eae_value)),
            spec$label
        )
        fisher_file <- sprintf(
            "%s_pooled_fisher_%s.csv",
            tolower(gsub(" ", "", eae_value)),
            spec$label
        )

        write_output_table(sample_res, sample_file)
        write_output_table(fisher_res, fisher_file)

        all_results[[length(all_results) + 1]] <- sample_res
        all_results[[length(all_results) + 1]] <- fisher_res
    }
}

combined_results <- bind_rows(all_results)
write_output_table(combined_results, "cell_abundance_stats_all_results.csv")

methods_path <- file.path(OUTPUT_DIR, "README.md")
writeLines(
    c(
        "# Cell Abundance Statistics",
        "",
        sprintf("**Date**: %s", Sys.Date()),
        "",
        "## What This Folder Contains",
        "",
        "This folder contains sample-aware abundance analyses for the same cluster and cell-type comparisons.",
        "",
        "## Why This Analysis Is Useful",
        "",
        "The project has biological replication (`rep` 1 to 3 in each group within each EAE stratum). Because of that, it is useful to also test abundance changes at the sample level rather than treating every cell as independent.",
        "",
        "## Data Used",
        "",
        "- Base object: `seurat_object.rds`",
        "- Stable clusters: from `clustassess_object.rds` (Most Abundant / 1,950 / SLM / 20)",
        "- Marker-based cell types: assigned with `marker_genes.R`",
        "- Cell-type tables exclude `Unassigned`, matching the current cell-type plotting workflow",
        "",
        "## Main Method",
        "",
        "- unit of replication: `sample_name`",
        "- comparison: `ctrl` vs `ko`, separately within `EAE` and `NO EAE`",
        "- model: quasibinomial GLM on in-category vs out-category counts per sample",
        "- multiple testing: BH correction within each result table",
        "",
        "## How To Read This Folder",
        "",
        "- `sample_cell_counts.csv`: how many cells were available per sample",
        "- `eae_sample_level_stable20.csv`: primary sample-level cluster analysis for EAE",
        "- `noeae_sample_level_stable20.csv`: primary sample-level cluster analysis for NO EAE",
        "- `eae_sample_level_celltypes.csv`: primary sample-level cell-type analysis for EAE",
        "- `noeae_sample_level_celltypes.csv`: primary sample-level cell-type analysis for NO EAE",
        "- `eae_pooled_fisher_*.csv` and `noeae_pooled_fisher_*.csv`: pooled-cell Fisher outputs included for direct comparison",
        "- `cell_abundance_stats_all_results.csv`: combined export",
        "",
        "## Interpretation Note",
        "",
        "In practice these sample-level results are more conservative than the pooled-cell Fisher tables, because they respect biological replication rather than counting each cell as an independent replicate.",
        "",
        "The pooled-cell Fisher outputs in this folder use the same row-normalized-to-100 convention as the dedicated Fisher script."
    ),
    con = methods_path
)
file.copy(methods_path, file.path(PROCESSED_DIR, "README.md"), overwrite = TRUE)

cat("Done.\n")
