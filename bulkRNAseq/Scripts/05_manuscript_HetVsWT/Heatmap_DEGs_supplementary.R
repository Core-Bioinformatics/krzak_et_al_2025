suppressPackageStartupMessages({
  library(dplyr)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(grid)
})

app_dir <- "C:/Users/Administrator/Documents/R/bulk_mRNA/NatImm_bulkAnalyseR_Macrophage_AEAE_HetVsWT"
de_dir  <- "C:/Users/Administrator/Documents/R/bulk_mRNA/NatImm_Macrophage_AEAE_HetVsWT_supplementary_outputs/02_DE"
out_dir <- "C:/Users/Administrator/Documents/R/bulk_mRNA/NatImm_Macrophage_AEAE_HetVsWT_supplementary_outputs/04_Heatmaps"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

deg.csv <- file.path(de_dir, "DE_genes_abslog2FC1_BH05_Het_vs_WT.csv")

MAX.GENES     <- NULL  # NULL = all DEGs (95 here); set eg 50 to cap, ranked by |log2FC|
CLUSTER.ROWS  <- FALSE # FALSE matches the app's convention; TRUE clusters within each direction block

stopifnot(file.exists(deg.csv))

load(file.path(app_dir, "expression_matrix.rda")) 
load(file.path(app_dir, "metadata.rda"))           

expr <- expression.matrix[[1]]
meta <- metadata[[1]]
stopifnot(identical(colnames(expr), meta$sample))

degs <- read.csv(deg.csv, stringsAsFactors = FALSE)
cat("DEGs (|log2FC|>=1, pvalAdj<=0.05):", nrow(degs), "\n")
stopifnot(all(degs$gene_id %in% rownames(expr)))
stopifnot(all(c("direction") %in% colnames(degs))) # written by the volcano/DE script

print(table(degs$direction))

degs <- degs %>% dplyr::arrange(dplyr::desc(abs(log2FC)))

if (!is.null(MAX.GENES) && nrow(degs) > MAX.GENES) {
  cat("Capping to top", MAX.GENES, "DEGs by |log2FC|.\n")
  degs <- utils::head(degs, MAX.GENES)
}

dup.names <- degs$gene_name[duplicated(degs$gene_name)]
if (length(dup.names) > 0) {
  cat("\n*** NOTE: duplicate gene_name, disambiguating with gene_id suffix:\n")
  print(unique(dup.names))
  degs$gene_name <- ifelse(degs$gene_name %in% dup.names,
                            paste0(degs$gene_name, " (", degs$gene_id, ")"),
                            degs$gene_name)
}

expr.subset <- expr[degs$gene_id, , drop = FALSE]
rownames(expr.subset) <- degs$gene_name

zmat <- t(scale(t(expr.subset)))
zmat[zmat > 3] <- 3
zmat[zmat < -3] <- -3

cat("Genes (rows):", nrow(zmat), " Samples (columns):", ncol(zmat), "\n")

include.exclude <- apply(meta, 2, function(x) {
  l <- length(unique(x))
  (l > 1) & (l < length(x))
})
if (sum(include.exclude == TRUE) != 0) {
  items <- colnames(meta)[include.exclude]
  items <- items[c(length(items), seq_len(length(items) - 1))]
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

zmat <- zmat[, as.character(meta_ordered[[1]]), drop = FALSE]
stopifnot(identical(colnames(zmat), as.character(meta_ordered[[1]])))

top.annotation.ids <- match(items, colnames(meta_ordered))

qual.col.pals <- dplyr::filter(RColorBrewer::brewer.pal.info, .data$category == "qual")
col.vector <- unique(unlist(mapply(RColorBrewer::brewer.pal,
                                    qual.col.pals$maxcolors,
                                    rownames(qual.col.pals))))
top.annotation.colour.list <- list()
colind <- 1
for (annos in seq_len(length(top.annotation.ids))) {
  values <- as.character(unique(meta_ordered[, top.annotation.ids[annos]]))
  vec <- vector(mode = "character")
  for (i in seq_len(length(values))) {
    vec <- c(vec, col.vector[colind])
    names(vec)[i] <- values[i]
    colind <- colind + 1
  }
  top.annotation.colour.list[[colnames(meta_ordered)[top.annotation.ids[annos]]]] <- vec
}

top.annotation.df <- as.data.frame(meta_ordered[, top.annotation.ids])
colnames(top.annotation.df) <- colnames(meta_ordered)[top.annotation.ids]
ha <- ComplexHeatmap::HeatmapAnnotation(
  df = top.annotation.df,
  col = top.annotation.colour.list,
  show_annotation_name = FALSE
)

print(top.annotation.colour.list)

direction.label <- ifelse(degs$log2FC > 0, "Higher in WT", "Higher in Het")
row_split <- factor(direction.label, levels = c("Higher in WT", "Higher in Het"))
print(table(row_split))

breaks  <- seq(-3, 3, 6 / 9)
colours <- rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu"))

deg.heatmap <- ComplexHeatmap::Heatmap(
  matrix              = zmat,
  name                = "Scale",
  col                 = circlize::colorRamp2(breaks = breaks, colors = colours),
  top_annotation      = ha,
  cluster_rows        = CLUSTER.ROWS,
  cluster_columns     = FALSE,
  row_split           = row_split,
  row_title_rot       = 0,
  row_names_side      = "left",
  show_row_names      = TRUE,
  show_column_names   = (ncol(zmat) <= 20),
  row_names_gp        = grid::gpar(fontsize = ifelse(nrow(zmat) > 40, 6, 8)),
  column_names_gp     = grid::gpar(fontsize = 10)
)

n.genes <- nrow(zmat)
fig.height <- max(5, min(30, 0.16 * n.genes + 2.5)) # a bit extra for the two row-title labels
fig.width  <- 7.5

pdf(file.path(out_dir, "Heatmap_DEGs_Het_vs_WT.pdf"), width = fig.width, height = fig.height)
ComplexHeatmap::draw(deg.heatmap)
dev.off()

png(file.path(out_dir, "Heatmap_DEGs_Het_vs_WT.png"),
    width = fig.width, height = fig.height, units = "in", res = 300)
ComplexHeatmap::draw(deg.heatmap)
dev.off()


write.csv(
  data.frame(gene = rownames(zmat), direction = direction.label, round(zmat, 4), check.names = FALSE),
  file.path(out_dir, "Heatmap_DEGs_Het_vs_WT_zscore_values.csv"),
  row.names = FALSE
)

cat("\nDone.\n")
