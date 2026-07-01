Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(bulkAnalyseR)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
})

base_dir <- "/iss-corescratch/se524/NatImm"
input_dir <- file.path(base_dir, "bulkAnalyseR_app", "NatImm_bulkAnalyseR_input")
app_dir <- file.path(base_dir, "bulkAnalyseR_app", "NatImm_bulkAnalyseR")

counts_file <- file.path(input_dir, "expression_counts_raw.tsv")
metadata_file <- file.path(input_dir, "metadata_clean_for_bulkAnalyseR.tsv")

dir.create(app_dir, recursive = TRUE, showWarnings = FALSE)

counts <- read.delim(counts_file, check.names = FALSE)
meta <- read.delim(metadata_file, check.names = FALSE)

cat("Counts dim raw:", dim(counts), "\n")
cat("Metadata dim:", dim(meta), "\n")

stopifnot(colnames(counts)[1] == "Geneid")
stopifnot(!anyDuplicated(counts$Geneid))
stopifnot(!any(is.na(counts$Geneid)))

expr <- as.matrix(counts[, -1, drop = FALSE])
rownames(expr) <- counts$Geneid
storage.mode(expr) <- "numeric"

stopifnot(colnames(meta)[1] == "sample")
stopifnot(identical(colnames(expr), meta$sample))

stopifnot(tail(colnames(meta), 1) == "condition")

cat("Expression matrix dim:", dim(expr), "\n")
cat("Metadata samples:", nrow(meta), "\n")
cat("Condition counts:\n")
print(table(meta$condition))

before <- nrow(expr)
expr <- expr[rowSums(expr) > 0, , drop = FALSE]
after <- nrow(expr)
cat("Genes before:", before, "\n")
cat("Genes after removing zero-sum genes:", after, "\n")

expr.proc <- tryCatch(
  {
    preprocessExpressionMatrix(expr, output.plot = FALSE)
  },
  error = function(e) {
    cat("preprocessExpressionMatrix failed. Falling back to CPM-like normalisation.\n")
    cat("Error was:\n")
    print(e)

    lib <- colSums(expr)
    expr.cpm <- sweep(expr, 2, lib, FUN = "/") * 1e6
    expr.log <- log2(expr.cpm + 1)
    expr.log
  }
)

cat("Processed expression matrix dim:", dim(expr.proc), "\n")

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
metadata <- list(meta)

save(expression.matrix, file = file.path(app_dir, "expression_matrix.rda"))
save(metadata, file = file.path(app_dir, "metadata.rda"))

app_text <- '
library(shiny)
library(dplyr)
library(ggplot2)
library(bulkAnalyseR)
library(shinythemes)

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

r.files <- list.files(path = getwd(), pattern = "\\\\.R$")
r.files <- setdiff(r.files, "app.R")
for(fl in r.files) source(fl)

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
    "NatImm bulk RNA-seq",
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
