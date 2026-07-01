suppressPackageStartupMessages({
  library(bulkAnalyseR)
  library(dplyr)
  library(ggplot2)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
})

app_dir <- "C:/Users/Administrator/Documents/R/bulk_mRNA/NatImm_bulkAnalyseR_Macrophage_AEAE_HetVsWT"
out_dir <- "C:/Users/Administrator/Documents/R/bulk_mRNA/NatImm_Macrophage_AEAE_HetVsWT_supplementary_outputs/02_DE"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

var1 <- "Het"   
var2 <- "WT"    

lfc.threshold  <- 1     # |log2FC| >= this
pval.threshold <- 0.05  # BH-adjusted p-value <= this

luca.previous.genes <- c("gpnmb", "dst", "mxd1", "anpep", "cytip", "lmna")

load(file.path(app_dir, "expression_matrix.rda")) 
load(file.path(app_dir, "metadata.rda"))          

expr <- expression.matrix[[1]]
meta <- metadata[[1]]

stopifnot(identical(colnames(expr), meta$sample))
stopifnot(tail(colnames(meta), 1) == "condition")
cat("Expression matrix dim:", dim(expr), "\n")
print(table(meta$condition))

anno <- AnnotationDbi::select(
  org.Mm.eg.db::org.Mm.eg.db,
  keys = rownames(expr),
  keytype = "ENSEMBL",
  columns = "SYMBOL"
) %>%
  dplyr::distinct(ENSEMBL, .keep_all = TRUE) %>%
  dplyr::mutate(NAME = ifelse(is.na(SYMBOL), ENSEMBL, SYMBOL))

condition.vec <- meta$condition
pos.var1 <- which(condition.vec == var1)[1]
pos.var2 <- which(condition.vec == var2)[1]

cat("condition vector:", paste(condition.vec, collapse = ", "), "\n")
cat("first position of var1 ('", var1, "'): ", pos.var1, "\n", sep = "")
cat("first position of var2 ('", var2, "'): ", pos.var2, "\n", sep = "")

if (pos.var1 < pos.var2) {
  cat("-> contrast <- c(-1, 1)  =>  log2FC = ", var2, " - ", var1, "\n", sep = "")
  direction.note <- paste0("POSITIVE log2FC = higher in ", var2, " | NEGATIVE log2FC = higher in ", var1)
} else {
  cat("-> contrast <- c(1, -1)  =>  log2FC = ", var1, " - ", var2, "\n", sep = "")
  direction.note <- paste0("POSITIVE log2FC = higher in ", var1, " | NEGATIVE log2FC = higher in ", var2)
}
cat(direction.note, "\n")

de.table <- bulkAnalyseR::DEanalysis_edger(
  expression.matrix = expr,
  condition          = condition.vec,
  var1               = var1,
  var2               = var2,
  anno               = anno
)

cat("Genes tested:", nrow(de.table), "\n")
cat("NA pvalAdj:", sum(is.na(de.table$pvalAdj)), "\n")

if (any(de.table$gene_name == "Mobp", na.rm = TRUE)) {
  mobp.row <- dplyr::filter(de.table, gene_name == "Mobp")
  cat("\nEmpirical check -- Mobp (myelin gene, expected high specifically in S50/WT",
      "if WT-direction is correctly positive):\n")
  print(mobp.row)
}

write.csv(de.table, file.path(out_dir, "DE_table_ALL_genes_Het_vs_WT.csv"), row.names = FALSE)

de.degs <- de.table %>%
  dplyr::filter(abs(log2FC) >= lfc.threshold, pvalAdj <= pval.threshold) %>%
  dplyr::mutate(
    direction = ifelse(log2FC > 0, paste0("higher_in_", var2), paste0("higher_in_", var1))
  ) %>%
  dplyr::arrange(pvalAdj)

cat("DEGs at |log2FC| >=", lfc.threshold, "& pvalAdj <=", pval.threshold)

write.csv(
  de.degs,
  file.path(out_dir, sprintf("DE_genes_abslog2FC%s_BH%s_Het_vs_WT.csv",
                              gsub("\\.", "", as.character(lfc.threshold)),
                              gsub("0\\.", "", as.character(pval.threshold)))),
  row.names = FALSE
)

luca.status <- lapply(luca.previous.genes, function(g) {
  idx <- which(tolower(de.table$gene_name) == tolower(g))
  if (length(idx) == 0) {
    return(data.frame(
      gene_queried = g, present_in_DE_table = FALSE, gene_name_matched = NA,
      log2FC = NA, direction = NA, pval = NA, pvalAdj = NA,
      passes_threshold = FALSE
    ))
  }
  row <- de.table[idx[1], ]
  passes <- abs(row$log2FC) >= lfc.threshold & row$pvalAdj <= pval.threshold
  data.frame(
    gene_queried = g, present_in_DE_table = TRUE, gene_name_matched = row$gene_name,
    log2FC = row$log2FC,
    direction = ifelse(row$log2FC > 0, paste0("higher_in_", var2), paste0("higher_in_", var1)),
    pval = row$pval, pvalAdj = row$pvalAdj,
    passes_threshold = passes
  )
}) %>% dplyr::bind_rows() %>%
  dplyr::mutate(
    provenance_note = paste0(
      "Luca's prior result (Gpnmb/Dst/Mxd1/Anpep/Cytip/Lmna) was generated on a GO-enrichment ",
      "run with effective_domain_size=8936, matching the OLD 39-sample app's processed gene ",
      "universe, not this renormalised 5-sample (9484-gene) subset. Differences from this ",
      "table are expected and not necessarily an error."
    )
  )



write.csv(luca.status, file.path(out_dir, "Luca_previous_genes_status_Het_vs_WT.csv"), row.names = FALSE)

genes.to.highlight <- luca.status %>%
  dplyr::filter(passes_threshold) %>%
  dplyr::pull(gene_name_matched)

cat("Genes that will actually be highlighted on the volcano",
    "(passed |log2FC|>=", lfc.threshold, "& pvalAdj<=", pval.threshold, ") ===\n")
if (length(genes.to.highlight) == 0) {
  cat("NONE of Luca's previously-mentioned genes pass the threshold in this analysi")
  cat("The volcano below will be produced WITHOUT custom highlights as a result")
} else {
  print(genes.to.highlight)
}

volcano <- bulkAnalyseR::volcano_plot(
  genes.de.results  = de.table,
  pval.threshold    = pval.threshold,
  lfc.threshold     = lfc.threshold,
  raster            = TRUE,                          
  log10pval.cap     = TRUE,                           
 (capPVal switch = off)
  add.labels.auto   = FALSE,                          
  add.labels.custom = length(genes.to.highlight) > 0,
  genes.to.label    = genes.to.highlight
) +
  ggplot2::labs(
    title = "Macrophage A-EAE: Het vs WT",
    subtitle = direction.note
  )

ggplot2::ggsave(file.path(out_dir, "Volcano_Het_vs_WT.pdf"), plot = volcano, width = 7, height = 5.5, dpi = 300)
ggplot2::ggsave(file.path(out_dir, "Volcano_Het_vs_WT.png"), plot = volcano, width = 7, height = 5.5, dpi = 300)

print(list.files(out_dir, full.names = TRUE))

cat("\nDone.\n")
