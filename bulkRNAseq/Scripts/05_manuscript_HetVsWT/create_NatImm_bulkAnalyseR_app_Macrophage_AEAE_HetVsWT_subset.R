Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1",
  R_DEFAULT_NUM_THREADS = "1"
)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(bulkAnalyseR)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(shinythemes)
})

base_dir  <- "C:/Users/Administrator/Documents/R/bulk_mRNA"
input_dir <- base_dir
app_dir   <- file.path(base_dir, "NatImm_bulkAnalyseR_Macrophage_AEAE_HetVsWT")

counts_file   <- file.path(input_dir, "expression_counts_raw.tsv")
metadata_file <- file.path(input_dir, "metadata_clean_for_bulkAnalyseR.tsv")

dir.create(app_dir, recursive = TRUE, showWarnings = FALSE)

het_samples <- c("S56", "S54", "S58")
wt_samples  <- c("S60", "S50")
selected_samples <- c(het_samples, wt_samples)

group_lookup <- setNames(
  c(rep("Het", length(het_samples)), rep("WT", length(wt_samples))),
  selected_samples
)

print(group_lookup)
cat("\nGroup sizes: Het n =", length(het_samples),
    "| WT n =", length(wt_samples), "\n")
cat("NOTE: very small group sizes -- treat any DE result here as\n")
cat("exploratory, not statistically robust.\n")

flagged_low_quality <- c("S50")
flagged_in_selection <- intersect(selected_samples, flagged_low_quality)
if (length(flagged_in_selection) > 0) {
  cat("Selected sample previously flagged 'low_mapping_or_low_assignment':\n")
  print(flagged_in_selection)
}

counts_full <- read.delim(counts_file, check.names = FALSE, stringsAsFactors = FALSE)
meta_full   <- read.delim(metadata_file, check.names = FALSE, stringsAsFactors = FALSE)

cat("Full counts dim:", dim(counts_full), "\n")
cat("Full metadata dim:", dim(meta_full), "\n")
print(colnames(meta_full))

known_bad_cols <- c("genotype", "matched_bio_replicate", "bio_replicate_id", "assigned_gene_counts")
present_bad <- intersect(known_bad_cols, colnames(meta_full))


stopifnot(colnames(counts_full)[1] == "Geneid")
stopifnot(!anyDuplicated(counts_full$Geneid))
stopifnot(!any(is.na(counts_full$Geneid)))
stopifnot(colnames(meta_full)[1] == "sample")
stopifnot(all(selected_samples %in% colnames(counts_full)))
stopifnot(all(selected_samples %in% meta_full$sample))

expr <- as.matrix(counts_full[, selected_samples, drop = FALSE])
rownames(expr) <- counts_full$Geneid
storage.mode(expr) <- "numeric"

meta_sub <- meta_full[match(selected_samples, meta_full$sample), , drop = FALSE]
rownames(meta_sub) <- NULL
stopifnot(identical(meta_sub$sample, colnames(expr)))

cat("Expression subset dim:", dim(expr), "\n")
cat("Metadata subset dim:", dim(meta_sub), "\n")

if ("cre_status" %in% colnames(meta_sub)) {
  check_df <- data.frame(
    sample              = meta_sub$sample,
    luca_group          = group_lookup[meta_sub$sample],
    metadata_cre_status = meta_sub$cre_status
  )
  print(check_df)
  mismatches <- as.character(check_df$luca_group) != as.character(check_df$metadata_cre_status)
  }

metadata_app <- data.frame(
  sample           = meta_sub$sample,
  cell_type        = meta_sub$cell_type,
  time_point       = meta_sub$time_point,
  tissue_timepoint = meta_sub$condition,                       # original cell_type_time_point label, kept for reference
  driver           = meta_sub$driver,
  cre_status       = meta_sub$cre_status,
  genotype_group   = paste(meta_sub$driver, meta_sub$cre_status, sep = "_"),
  mapping_qc_flag  = meta_sub$mapping_qc_flag,
  condition        = unname(group_lookup[meta_sub$sample]),    # Het / WT -- last column, used by bulkAnalyseR for DE
  stringsAsFactors = FALSE
)

stopifnot(identical(metadata_app$sample, colnames(expr)))
stopifnot(tail(colnames(metadata_app), 1) == "condition")
stopifnot(!anyNA(metadata_app$condition))

other_cols <- setdiff(colnames(metadata_app), "condition")
for (col in other_cols) {
  metadata_app[[col]][is.na(metadata_app[[col]])] <- "Unknown"
  metadata_app[[col]][metadata_app[[col]] == ""] <- "Unknown"
  metadata_app[[col]][metadata_app[[col]] == "NA"] <- "Unknown"
}

print(metadata_app)
cat("\nCondition (DE) counts:\n")
print(table(metadata_app$condition))
cat("\nMapping QC flag counts:\n")
print(table(metadata_app$mapping_qc_flag))

before <- nrow(expr)
expr <- expr[rowSums(expr) > 0, , drop = FALSE]
after <- nrow(expr)
cat("Genes before:", before, "\n")
cat("Genes after removing zero-sum genes:", after, "\n")

expr.proc <- preprocessExpressionMatrix(
  expr,
  output.plot = TRUE,
  normalisation.method = "quantile"
)

cat("Processed expression matrix dim:", dim(expr.proc), "\n")
cat("NA values:", sum(is.na(expr.proc)), "\n")
cat("Duplicated rownames:", anyDuplicated(rownames(expr.proc)), "\n")
cat("Range:", range(expr.proc), "\n")

stopifnot(identical(colnames(expr.proc), metadata_app$sample))

anno <- AnnotationDbi::select(
  org.Mm.eg.db::org.Mm.eg.db,
  keys = rownames(expr.proc),
  keytype = "ENSEMBL",
  columns = "SYMBOL"
) %>%
  dplyr::distinct(ENSEMBL, .keep_all = TRUE) %>%
  dplyr::mutate(NAME = ifelse(is.na(SYMBOL), ENSEMBL, SYMBOL))

cat("Annotation dim:", dim(anno), "\n")

expression.matrix <- list(expr.proc)
metadata <- list(metadata_app)

save(expression.matrix, file = file.path(app_dir, "expression_matrix.rda"))
save(metadata, file = file.path(app_dir, "metadata.rda"))

write.table(
  data.frame(Geneid = rownames(expr), expr, check.names = FALSE),
  file = file.path(app_dir, "subset_raw_counts_before_preprocessing.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  metadata_app,
  file = file.path(app_dir, "metadata_used_in_app.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

app_text <- '
library(shiny)
library(dplyr)
library(ggplot2)
library(bulkAnalyseR)
library(shinythemes)
library(AnnotationDbi)
library(org.Mm.eg.db)

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

rda.files <- list.files(pattern = "\\\\.rda$")
for(fl in rda.files) load(fl)

anno <- list()
anno[[1]] <- AnnotationDbi::select(
  getExportedValue("org.Mm.eg.db", "org.Mm.eg.db"),
  keys = rownames(expression.matrix[[1]]),
  keytype = "ENSEMBL",
  columns = "SYMBOL"
) %>%
  dplyr::distinct(ENSEMBL, .keep_all = TRUE) %>%
  dplyr::mutate(NAME = ifelse(is.na(SYMBOL), ENSEMBL, SYMBOL))

panels <- c(
  "Landing", "SampleSelect", "QC", "GRN", "DE", "DEplot",
  "DEsummary", "Enrichment", "GRNenrichment", "Cross", "Patterns"
)

ui <- function(request){
  navbarPage(
    "NatImm Macrophage A-EAE Het vs WT (subset)",
    theme = shinythemes::shinytheme("flatly"),
    header = tags$head(tags$style("body {overflow-y: scroll;}")),
    footer = bookmarkButton(),
    tabPanel(
      title = "RNA",
      modalityPanelUI(
        id = "RNA",
        metadata = metadata[[1]],
        organism = "mmusculus",
        panels.default = panels
      )
    )
  )
}

server <- function(input, output, session){
  modalityPanelServer(
    id = "RNA",
    expression.matrix = expression.matrix[[1]],
    metadata = metadata[[1]],
    anno = anno[[1]],
    organism = "mmusculus",
    panels.default = panels
  )
}

shinyApp(ui, server, enableBookmarking = "url")
'

writeLines(app_text, file.path(app_dir, "app.R"))

print(list.files(app_dir, full.names = TRUE))
cat("\nDone.\n")
