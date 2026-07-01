suppressPackageStartupMessages({
  library(bulkAnalyseR)
  library(dplyr)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(grid)
})

app_dir <- "C:/Users/Administrator/Documents/R/bulk_mRNA/NatImm_bulkAnalyseR_Macrophage_AEAE_HetVsWT"
out_dir <- "C:/Users/Administrator/Documents/R/bulk_mRNA/NatImm_Macrophage_AEAE_HetVsWT_supplementary_outputs/04_Heatmaps"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

raw.counts.file <- file.path(app_dir, "subset_raw_counts_before_preprocessing.tsv")

heatmap.type <- "Z-score"    are also saved separately below
CLUSTER.ROWS <- FALSE       
luca.genes.raw <- c(
  "Il1b", "Il1a", "il6", "il10", "il4", "Cxcr4", "Cxcr2", "Cxcr3", "Cxcr6",
  "Cxcl2", "Cxcl9", "Cxcl10", "Cxcl16", "Nos2", "Arg1", "Ccr8", "Ccr2",
  "Ccr5", "Ccr1", "Ccl1", "Ccl5", "Ccl9", "Ccl6", "Ccl4", "Ccl2",
  "Il1rn", "Mrc1", "Enpp2"
)
stopifnot(length(luca.genes.raw) == 28)

gene.category <- c(
  rep("Interleukins", 5),                    # Il1b, Il1a, il6, il10, il4
  rep("CXCR receptors", 4),                  # Cxcr4, Cxcr2, Cxcr3, Cxcr6
  rep("CXCL chemokines", 4),                 # Cxcl2, Cxcl9, Cxcl10, Cxcl16
  rep("Macrophage activation (M1/M2)", 2),   # Nos2, Arg1
  rep("CCR receptors", 4),                   # Ccr8, Ccr2, Ccr5, Ccr1
  rep("CCL chemokines", 6),                  # Ccl1, Ccl5, Ccl9, Ccl6, Ccl4, Ccl2
  rep("Regulatory / macrophage markers", 3)  # Il1rn, Mrc1, Enpp2
)
stopifnot(length(gene.category) == 28)
names(gene.category) <- luca.genes.raw

luca.genes <- paste0(toupper(substr(luca.genes.raw, 1, 1)),
                      tolower(substr(luca.genes.raw, 2, nchar(luca.genes.raw))))
names(luca.genes) <- luca.genes.raw # keep the original spelling for reporting

print(luca.genes)

stopifnot(file.exists(raw.counts.file))
raw.df <- read.delim(raw.counts.file, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(colnames(raw.df)[1] == "Geneid")
raw.counts <- as.matrix(raw.df[, -1, drop = FALSE])
rownames(raw.counts) <- raw.df$Geneid
storage.mode(raw.counts) <- "numeric"

load(file.path(app_dir, "expression_matrix.rda")) 
load(file.path(app_dir, "metadata.rda"))          
processed.matrix <- expression.matrix[[1]]
meta <- metadata[[1]]

stopifnot(identical(colnames(raw.counts), meta$sample))
stopifnot(identical(colnames(processed.matrix), meta$sample))

cat("Raw (pre-noisyR, post zero-sum filter):", nrow(raw.counts), "genes\n")
cat("Main processed (noisyR + qnorm):", nrow(processed.matrix), "genes\n")

qnorm.no.noisyr <- bulkAnalyseR::preprocessExpressionMatrix(
  raw.counts,
  denoise               = FALSE,
  output.plot           = FALSE,
  normalisation.method  = "quantile"
)
cat("qnorm-without-noisyR matrix dim:", dim(qnorm.no.noisyr), "\n")

symbol.map <- AnnotationDbi::select(
  org.Mm.eg.db::org.Mm.eg.db,
  keys     = unique(luca.genes),
  keytype  = "SYMBOL",
  columns  = "ENSEMBL"
)
print(symbol.map)

multi.map <- symbol.map %>% dplyr::count(SYMBOL) %>% dplyr::filter(n > 1)
if (nrow(multi.map) > 0) {
  
  print(multi.map$SYMBOL)
}

status.list <- lapply(names(luca.genes), function(orig.name) {
  sym <- luca.genes[[orig.name]]
  candidates <- symbol.map$ENSEMBL[symbol.map$SYMBOL == sym & !is.na(symbol.map$ENSEMBL)]

  if (length(candidates) == 0) {
    return(data.frame(
      gene_queried = orig.name, gene_symbol_used = sym, ensembl_id = NA,
      present_in_raw_counts = FALSE, present_in_qnorm_no_noisyR = FALSE,
      present_in_main_processed_noisyR_qnorm = FALSE,
      notes = "No ENSEMBL ID found for this symbol in org.Mm.eg.db"
    ))
  }

  found.in.raw <- candidates[candidates %in% rownames(raw.counts)]
  if (length(found.in.raw) == 0) {
    return(data.frame(
      gene_queried = orig.name, gene_symbol_used = sym,
      ensembl_id = paste(candidates, collapse = ";"),
      present_in_raw_counts = FALSE, present_in_qnorm_no_noisyR = FALSE,
      present_in_main_processed_noisyR_qnorm = FALSE,
      notes = "ENSEMBL ID(s) exist but none present in this 5-sample subset (zero-sum across all 5 samples, removed before noisyR even ran)"
    ))
  }

  used.id <- found.in.raw[1]
  multi.note <- if (length(found.in.raw) > 1) {
    paste0("Multiple ENSEMBL IDs present (", paste(found.in.raw, collapse = ";"), "); using ", used.id)
  } else {
    ""
  }

  data.frame(
    gene_queried = orig.name, gene_symbol_used = sym, ensembl_id = used.id,
    present_in_raw_counts = TRUE,
    present_in_qnorm_no_noisyR = used.id %in% rownames(qnorm.no.noisyr),
    present_in_main_processed_noisyR_qnorm = used.id %in% rownames(processed.matrix),
    notes = multi.note
  )
}) %>% dplyr::bind_rows()

status.list <- status.list %>%
  dplyr::mutate(
    removed_by_noisyR = present_in_qnorm_no_noisyR & !present_in_main_processed_noisyR_qnorm
  )

print(status.list)

cat("\nSummary: ", sum(status.list$present_in_qnorm_no_noisyR), "/28 present in qnorm-no-noisyR matrix | ",
    sum(status.list$removed_by_noisyR), " of those were removed by noisyR from the main matrix\n", sep = "")

write.csv(status.list, file.path(out_dir, "Luca_selected_genes_present_absent.csv"), row.names = FALSE)

present.genes <- status.list %>% dplyr::filter(present_in_qnorm_no_noisyR)

if (nrow(present.genes) == 0) {
  stop("None of Luca's 28 genes are present in the qnorm-without-noisyR matrix - check gene symbols / data.")
}



expr.subset <- qnorm.no.noisyr[present.genes$ensembl_id, , drop = FALSE]
rownames(expr.subset) <- present.genes$gene_queried # preserve Luca's original spelling/order for display

row.category.subset <- gene.category[present.genes$gene_queried]
row_split <- factor(row.category.subset, levels = unique(gene.category))
print(table(row_split))

write.csv(
  data.frame(gene = rownames(expr.subset), round(expr.subset, 2), check.names = FALSE),
  file.path(out_dir, "Luca_selected_genes_qnorm_absolute_values.csv"),
  row.names = FALSE
)

zmat <- t(scale(t(expr.subset)))
zmat[zmat > 3] <- 3
zmat[zmat < -3] <- -3

write.csv(
  data.frame(gene = rownames(zmat), category = row.category.subset, round(zmat, 4), check.names = FALSE),
  file.path(out_dir, "Heatmap_Luca_selected_genes_zscore_values.csv"),
  row.names = FALSE
)

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

meta_f <- as.data.frame(
  lapply(meta, function(x) if (!is.factor(x)) factor(x, levels = unique(x)) else x),
  stringsAsFactors = FALSE
)
colnames(meta_f) <- colnames(meta)
meta_ordered <- meta_f %>% dplyr::arrange(dplyr::across(dplyr::all_of(items)))

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

breaks  <- seq(-3, 3, 6 / 9)
colours <- rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu"))

luca.heatmap <- ComplexHeatmap::Heatmap(
  matrix              = zmat,
  name                = "Scale",
  col                 = circlize::colorRamp2(breaks = breaks, colors = colours),
  top_annotation      = ha,
  cluster_rows        = CLUSTER.ROWS,
  cluster_columns     = FALSE,
  row_split           = row_split,
  row_title_rot       = 0,
  row_title_gp        = grid::gpar(fontsize = 8),
  row_names_side      = "left",
  show_row_names      = TRUE,
  show_column_names   = (ncol(zmat) <= 20),
  row_names_gp        = grid::gpar(fontsize = 9),
  column_names_gp     = grid::gpar(fontsize = 10)
)

fig.height <- max(5, 0.24 * nrow(zmat) + 2.5) # extra room for the 7 category row-title labels
fig.width  <- 7.5

pdf(file.path(out_dir, "Heatmap_Luca_selected_genes_Het_vs_WT.pdf"), width = fig.width, height = fig.height)
ComplexHeatmap::draw(luca.heatmap)
dev.off()

png(file.path(out_dir, "Heatmap_Luca_selected_genes_Het_vs_WT.png"),
    width = fig.width, height = fig.height, units = "in", res = 300)
ComplexHeatmap::draw(luca.heatmap)
dev.off()

cat("\n Saved (", nrow(zmat), " of 28 genes plotted, height ", round(fig.height, 1), "in) ===\n", sep = "")
print(list.files(out_dir, pattern = "Luca_selected|Heatmap_Luca", full.names = TRUE))

cat("\nDone.\n")
