suppressPackageStartupMessages({
  library(bulkAnalyseR)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(dplyr)
  library(grid)
})

app_dir <- "C:/Users/Administrator/Documents/R/bulk_mRNA/NatImm_bulkAnalyseR_Macrophage_AEAE_HetVsWT"
out_dir <- "C:/Users/Administrator/Documents/R/bulk_mRNA/NatImm_Macrophage_AEAE_HetVsWT_supplementary_outputs/01_QC"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

n.abundant <- 500  
show.values <- FALSE 
load(file.path(app_dir, "expression_matrix.rda")) 
load(file.path(app_dir, "metadata.rda"))           
expr <- expression.matrix[[1]]
meta <- metadata[[1]]

cat("Expression matrix dim:", dim(expr), "\n")
cat("Metadata dim:", dim(meta), "\n")
stopifnot(identical(colnames(expr), meta$sample))

include.exclude <- apply(meta, 2, function(x) {
  l <- length(unique(x))
  (l > 1) & (l < length(x))
})

if (sum(include.exclude == TRUE) != 0) {
  items <- colnames(meta)[include.exclude]
  items <- items[c(length(items), seq_len(length(items) - 1))] # rotate: last column (condition) shown first
} else {
  items <- colnames(meta)[2:ncol(meta)]
}

print(items)

meta_f <- as.data.frame(
  lapply(meta, function(x) if (!is.factor(x)) factor(x, levels = unique(x)) else x),
  stringsAsFactors = FALSE
)
colnames(meta_f) <- colnames(meta)

meta_ordered <- meta_f %>% dplyr::arrange(dplyr::across(dplyr::all_of(items)))

print(as.character(meta_ordered[[1]]))

expr_ordered <- expr[, as.character(meta_ordered[[1]]), drop = FALSE]
stopifnot(identical(colnames(expr_ordered), as.character(meta_ordered[[1]])))


jsi_plot <- bulkAnalyseR::jaccard_heatmap(
  expression.matrix      = expr_ordered,
  metadata                = meta_ordered,
  top.annotation.ids      = match(items, colnames(meta_ordered)),
  n.abundant              = n.abundant,
  show.values             = show.values,
  show.row.column.names   = (nrow(meta_ordered) <= 20)
)

pdf(file.path(out_dir, "JSI_heatmap_Het_vs_WT.pdf"), width = 7, height = 5)
ComplexHeatmap::draw(jsi_plot)
dev.off()

png(file.path(out_dir, "JSI_heatmap_Het_vs_WT.png"), width = 7, height = 5, units = "in", res = 300)
ComplexHeatmap::draw(jsi_plot)
dev.off()

print(list.files(out_dir, pattern = "JSI", full.names = TRUE))

n.abundant.used <- min(n.abundant, nrow(expr_ordered))
samples <- colnames(expr_ordered)
jsi.matrix <- matrix(0, nrow = length(samples), ncol = length(samples),
                      dimnames = list(samples, samples))
for (i in seq_along(samples)) {
  for (j in seq_len(i)) {
    gi <- order(expr_ordered[, i], decreasing = TRUE)[1:n.abundant.used]
    gj <- order(expr_ordered[, j], decreasing = TRUE)[1:n.abundant.used]
    jsi.matrix[i, j] <- jsi.matrix[j, i] <- bulkAnalyseR::jaccard_index(gi, gj)
  }
}

write.csv(
  jsi.matrix,
  file.path(out_dir, "JSI_values_Het_vs_WT.csv"),
  row.names = TRUE
)

print(round(jsi.matrix, 3))

cat("\nDone.\n")
