suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

app_dir <- "C:/Users/Administrator/Documents/R/bulk_mRNA/NatImm_bulkAnalyseR_Macrophage_AEAE_HetVsWT"
out_dir <- "C:/Users/Administrator/Documents/R/bulk_mRNA/NatImm_Macrophage_AEAE_HetVsWT_supplementary_outputs/01_QC"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

n.abundant     <- 500
annotation.col <- "condition"

col.Het <- "#7FC97F"
col.WT  <- "#BEAED4"

load(file.path(app_dir, "expression_matrix.rda"))
load(file.path(app_dir, "metadata.rda"))

expr <- expression.matrix[[1]]
meta <- metadata[[1]]

stopifnot(identical(colnames(expr), meta$sample))
stopifnot(annotation.col %in% colnames(meta))

cat("Expression matrix dim:", dim(expr), "\n")
print(table(meta[[annotation.col]]))

n.abundant <- min(n.abundant, nrow(expr))
top.idx    <- utils::tail(order(rowSums(expr)), n.abundant)
pca.input  <- t(expr[top.idx, , drop = FALSE])

gene.sd <- apply(pca.input, 2, sd, na.rm = TRUE)
keep    <- gene.sd > 0 & !is.na(gene.sd)
cat("\nGenes before zero-variance filter:", ncol(pca.input), "\n")
pca.input <- pca.input[, keep, drop = FALSE]
cat("Genes after zero-variance filter:", ncol(pca.input), "\n")

pca.fit  <- stats::prcomp(pca.input, center = TRUE, scale. = TRUE)
var.expl <- 100 * (pca.fit$sdev^2 / sum(pca.fit$sdev^2))
pc1.pct  <- round(var.expl[1], 1)
pc2.pct  <- round(var.expl[2], 1)
cat("\nPC1:", pc1.pct, "% - PC2:", pc2.pct, "%\n")

pca.df <- data.frame(
  sample = rownames(pca.fit$x),
  PC1    = pca.fit$x[, 1],
  PC2    = pca.fit$x[, 2],
  stringsAsFactors = FALSE
) %>%
  dplyr::left_join(meta, by = "sample") %>%
  dplyr::mutate(
    condition = factor(condition, levels = c("Het", "WT")),
    mapping_qc_flag = factor(
      mapping_qc_flag,
      levels = c("pass_or_moderate", "low_mapping_or_low_assignment")
    )
  )

print(pca.df[, c("sample", "condition", "mapping_qc_flag", "PC1", "PC2")])

padding <- 1.10

pc1.lim <- max(abs(pca.df$PC1)) * padding
pc2.lim <- max(abs(pca.df$PC2)) * padding

cat("PC1: ±", round(pc1.lim, 2), "\n")
cat("PC2: ±", round(pc2.lim, 2), "\n")

write.csv(
  pca.df[, c("sample", "condition", "mapping_qc_flag", "PC1", "PC2")],
  file.path(out_dir, "PCA_Het_vs_WT_coordinates.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    PC = paste0("PC", seq_along(var.expl)),
    variance_explained_percent = round(var.expl, 4)
  ),
  file.path(out_dir, "PCA_Het_vs_WT_variance_explained.csv"),
  row.names = FALSE
)

writeLines(
  c(
    "PCA Het vs WT — summary (v2 symmetric axes)",
    "",
    "Samples: Het = S56, S54, S58 ; WT = S60, S50",
    paste0("Genes used: top ", n.abundant, " by rowSums, ",
           ncol(pca.input), " after zero-variance filter"),
    "PCA: prcomp(center = TRUE, scale. = TRUE)",
    paste0("PC1: ", pc1.pct, "% ; PC2: ", pc2.pct, "%"),
    "",
    "Axes: symmetric around zero, independently per axis.",
    paste0("  PC1: ±", round(pc1.lim, 2),
           "  (driven by S50 at PC1 = ", round(max(pca.df$PC1), 2), ")"),
    paste0("  PC2: ±", round(pc2.lim, 2),
           "  (driven by S60 at PC2 = ", round(max(pca.df$PC2), 2), ")"),
    "",
    "Groups differentiated by colour only.",
    "S50 flagged as low_mapping_or_low_assignment (shown as triangle)."
  ),
  file.path(out_dir, "PCA_Het_vs_WT_summary.txt")
)

x.lab <- paste0("PC1 (proportion of variance = ", pc1.pct, "%)")
y.lab <- paste0("PC2 (proportion of variance = ", pc2.pct, "%)")

base_plot <- ggplot(pca.df, aes(x = PC1, y = PC2)) +
  geom_point(
    aes(colour = condition, shape = mapping_qc_flag),
    size = 4, stroke = 1.1
  ) +
  scale_colour_manual(
    name   = "Condition",
    values = c("Het" = col.Het, "WT" = col.WT)
  ) +
  scale_shape_manual(
    name   = "Mapping QC",
    values = c(
      "pass_or_moderate"              = 16,
      "low_mapping_or_low_assignment" = 17
    )
  ) +
  scale_x_continuous(limits = c(-pc1.lim, pc1.lim)) +
  scale_y_continuous(limits = c(-pc2.lim, pc2.lim)) +
  labs(x = x.lab, y = y.lab) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position  = "right",
    plot.title       = element_text(face = "bold")
  )

label_layer <- ggrepel::geom_label_repel(
  data    = pca.df,
  mapping = aes(x = PC1, y = PC2, colour = condition, label = sample),
  size               = 3.6,
  show.legend        = FALSE,
  max.overlaps       = Inf,
  box.padding        = 0.35,
  point.padding      = 0.3,
  min.segment.length = 0,
  label.size         = 0.2,
  xlim = c(-pc1.lim, pc1.lim),
  ylim = c(-pc2.lim, pc2.lim)
)

title.txt    <- "PCA \u2014 Macrophage A-EAE: Het vs WT"
subtitle.txt <- "Groups differentiated by colour; no ellipses (WT n=2)"

save_variant <- function(p, stem, w = 7, h = 5.5) {
  ggplot2::ggsave(file.path(out_dir, paste0(stem, ".pdf")),
                  plot = p, width = w, height = h)
  ggplot2::ggsave(file.path(out_dir, paste0(stem, ".png")),
                  plot = p, width = w, height = h, dpi = 300)
  cat("Saved:", stem, "\n")
}

save_variant(
  base_plot + label_layer + labs(title = title.txt, subtitle = subtitle.txt),
  "PCA_Het_vs_WT_symmetric_title_labels"
)

save_variant(
  base_plot + labs(title = title.txt, subtitle = subtitle.txt),
  "PCA_Het_vs_WT_symmetric_title_nolabels"
)

save_variant(
  base_plot + label_layer,
  "PCA_Het_vs_WT_symmetric_notitle_labels"
)

save_variant(
  base_plot,
  "PCA_Het_vs_WT_symmetric_notitle_nolabels"
)

print(list.files(out_dir, pattern = "PCA_Het_vs_WT_symmetric", full.names = TRUE))
cat("\nDone.\n")
