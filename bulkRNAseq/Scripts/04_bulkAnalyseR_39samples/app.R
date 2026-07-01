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

r.files <- list.files(path = getwd(), pattern = "\\.R$")
r.files <- setdiff(r.files, "app.R")
for(fl in r.files) source(fl)

rda.files <- list.files(pattern = "\\.rda$")
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
