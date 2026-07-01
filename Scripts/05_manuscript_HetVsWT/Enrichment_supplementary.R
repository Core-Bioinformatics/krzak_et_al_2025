suppressPackageStartupMessages({
  library(dplyr)
  library(gprofiler2)
})

app_dir <- "C:/Users/Administrator/Documents/R/bulk_mRNA/NatImm_bulkAnalyseR_Macrophage_AEAE_HetVsWT"
de_dir  <- "C:/Users/Administrator/Documents/R/bulk_mRNA/NatImm_Macrophage_AEAE_HetVsWT_supplementary_outputs/02_DE"
out_dir <- "C:/Users/Administrator/Documents/R/bulk_mRNA/NatImm_Macrophage_AEAE_HetVsWT_supplementary_outputs/03_Enrichment"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

deg.csv    <- file.path(de_dir, "DE_genes_abslog2FC1_BH05_Het_vs_WT.csv")
all.de.csv <- file.path(de_dir, "DE_table_ALL_genes_Het_vs_WT.csv")

sources.app <- c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "TF", "MIRNA")
sources.go  <- c("GO:BP", "GO:MF", "GO:CC")

stopifnot(file.exists(deg.csv), file.exists(all.de.csv))

degs   <- read.csv(deg.csv,    stringsAsFactors = FALSE)
all.de <- read.csv(all.de.csv, stringsAsFactors = FALSE)

cat("DEGs (query):", nrow(degs), "\n")
cat("All tested genes (custom_bg):", nrow(all.de), "\n")

if (!"direction" %in% colnames(degs)) {
  cat("NOTE: 'direction' column not found -- computing from log2FC.\n")
  stopifnot("log2FC" %in% colnames(degs))
  degs <- degs %>%
    dplyr::mutate(
      direction = dplyr::case_when(
        log2FC > 0 ~ "higher_in_WT",
        log2FC < 0 ~ "higher_in_Het",
        TRUE ~ "no_change"
      )
    )
}

cat("\nDEG direction counts:\n")
print(table(degs$direction))

degs.wt  <- dplyr::filter(degs, direction == "higher_in_WT")
degs.het <- dplyr::filter(degs, direction == "higher_in_Het")

custom.bg <- all.de$gene_id

cat("\nQuery sizes: all=", nrow(degs),
    "| up-in-WT=", nrow(degs.wt),
    "| up-in-Het=", nrow(degs.het), "\n")
cat("Background:", length(custom.bg), "genes\n\n")

run_gost <- function(query.ids, label, sources) {

  cat("Enrichment:", label, "(n =", length(query.ids), ", sources:", paste(sources, collapse = "+"), ") ")

  if (length(query.ids) < 2) {
    cat("Too few genes (<2) - skipping.\n\n")
    return(NULL)
  }

  res <- tryCatch(
    gprofiler2::gost(
      query             = query.ids,
      organism          = "mmusculus",
      correction_method = "fdr",      
      custom_bg         = custom.bg,
      sources           = sources,
      evcodes           = TRUE          
    ),
    error = function(e) {
      cat("ERROR:", conditionMessage(e), "\n\n")
      return(NULL)
    }
  )

  if (is.null(res) || is.null(res$result) || nrow(res$result) == 0) {
    cat("No enriched terms found.\n\n")
    return(data.frame())
  }

  result <- res$result %>%
    dplyr::mutate(
      parents = sapply(.data$parents, toString),
      intersection_names = sapply(.data$intersection, function(x) {
        ensids <- strsplit(x, split = ",")[[1]]
        syms   <- all.de$gene_name[match(ensids, all.de$gene_id)]
        paste(syms, collapse = ",")
      }),
      query_label = label
    )

  cat("Enriched terms:", nrow(result), "\n\n")
  return(result)
}

save_result <- function(res, filename) {
  path <- file.path(out_dir, filename)
  if (is.null(res) || nrow(res) == 0) {
    writeLines(
      c(paste("No enriched terms found."),
        paste("File:", filename)),
      sub("\\.csv$", "_NO_RESULTS.txt", path)
    )
    cat("No results - wrote:", sub("\\.csv$", "_NO_RESULTS.txt", filename), "\n")
  } else {
    write.csv(res, path, row.names = FALSE)
    cat("Saved:", filename, "(", nrow(res), "terms)\n")
  }
}

res.all.sources <- run_gost(degs$gene_id, "all_DEGs_all_sources", sources.app)
save_result(res.all.sources, "Enrichment_master_all_DEGs_Het_vs_WT.csv")

res.go.all <- run_gost(degs$gene_id,    "all_DEGs_GO_only",      sources.go)
res.go.wt  <- run_gost(degs.wt$gene_id, "higher_in_WT_GO_only",  sources.go)
res.go.het <- run_gost(degs.het$gene_id,"higher_in_Het_GO_only", sources.go)

save_result(res.go.all, "GO_enrichment_all_DEGs_Het_vs_WT.csv")
save_result(res.go.wt,  "GO_enrichment_up_in_WT_Het_vs_WT.csv")
save_result(res.go.het, "GO_enrichment_up_in_Het_Het_vs_WT.csv")

count_terms <- function(x) if (is.null(x) || nrow(x) == 0) 0L else nrow(x)

writeLines(
  c(
    "Enrichment analysis - summary",
    "",
    "Comparison: Macrophage A-EAE Het vs WT",
    "Samples: Het = S56, S54, S58 - WT = S60, S50",
    "",
    "DEG threshold: |log2FC|>=1, pvalAdj<=0.05 (BH)",
    "log2FC direction: positive = higher in WT - negative = higher in Het",
    "",
    paste0("DEGs total: ", nrow(degs),
           " - up-in-WT: ", nrow(degs.wt),
           " - up-in-Het: ", nrow(degs.het)),
    paste0("Background (custom_bg): ", length(custom.bg), " genes (all tested in 5-sample DE)"),
    "",
    "gprofiler2::gost() parameters (verified against bulkAnalyseR source):",
    "  organism = mmusculus",
    "  correction_method = fdr",
    "  evcodes = TRUE",
    "",
    "Runs:",
    paste0("  All-source (", paste(sources.app, collapse=","), "): ", count_terms(res.all.sources), " terms"),
    paste0("  GO-only all DEGs: ", count_terms(res.go.all), " terms"),
    paste0("  GO-only up-in-WT: ", count_terms(res.go.wt), " terms"),
    paste0("  GO-only up-in-Het: ", count_terms(res.go.het), " terms"),
    "",
    "Note on GO-only separate runs:",
    "  GO enrichments (for Luca) were re-run with sources=GO:BP/MF/CC only,",
    "  NOT filtered from the all-source run. This gives GO-specific FDR",
    "  correction (more appropriate for a focused GO analysis).",
    "",
    "Note on background:",
    "  custom_bg = 9484 genes (our 5-sample renormalised analysis).",
    "  Luca's previous result used effective_domain_size=8936 (old 39-sample",
    "  app, un-renormalised). Differences in results are expected.",
    "",
    "S50 caveat: included; flagged low_mapping_or_low_assignment."
  ),
  file.path(out_dir, "Enrichment_summary_Het_vs_WT.txt")
)

print(list.files(out_dir, full.names = TRUE))
cat("\nDone.\n")
