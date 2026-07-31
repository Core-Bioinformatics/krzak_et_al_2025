library(Seurat)
library(CellChat)
library(dplyr)
library(readr)
library(circlize)
library(future)

future::plan("sequential")
options(future.globals.maxSize = 8 * 1024^3)
set.seed(42)

object_file <- "/iss-scratch/CoreBioinformatics/rk720/nichenet_stuff/sucnr1_nichenet_object.rds"
output_dir <- "/iss-scratch/CoreBioinformatics/rk720/nichenet_stuff/cellchat_eae_ctrl_vs_ko_v2"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

condition_col <- "eae_condition"
celltype_col <- "celltype_markers"
sample_col <- "orig.ident"
conditions <- c(WT = "EAE_ctrl", KO = "EAE_ko")
min_cells <- 20

obj <- readRDS(object_file)

obj@meta.data %>%
  dplyr::filter(eae_condition %in% c("EAE_ctrl", "EAE_ko"), celltype_markers == "Neutrophils") %>%
  dplyr::count(eae_condition, orig.ident)
DefaultAssay(obj) <- "RNA"

if (inherits(obj[["RNA"]], "Assay5")) obj <- SeuratObject::JoinLayers(obj, assay = "RNA")

obj <- subset(obj, cells = colnames(obj)[
  obj@meta.data[[condition_col]] %in% conditions &
  obj@meta.data[[celltype_col]] != "Unassigned" &
  !is.na(obj@meta.data[[celltype_col]])
])

obj <- NormalizeData(obj, assay = "RNA", verbose = FALSE)
obj$cellchat_group <- as.character(obj@meta.data[[celltype_col]])
obj$samples <- factor(obj@meta.data[[sample_col]])

cell_counts <- table(obj$cellchat_group, obj@meta.data[[condition_col]])

cell_counts
valid_groups <- rownames(cell_counts)[
  apply(cell_counts[, conditions, drop = FALSE] >= min_cells, 1, all)
]

obj <- subset(obj, cells = colnames(obj)[obj$cellchat_group %in% valid_groups])
obj$cellchat_group <- factor(obj$cellchat_group, levels = valid_groups)

table(obj$celltype_markers)

write_csv(
  as.data.frame(cell_counts) %>% rename(celltype = Var1, condition = Var2, n_cells = Freq),
  file.path(output_dir, "cell_counts.csv")
)

message("Cell types retained: ", paste(valid_groups, collapse = ", "))


# CellChat database
CellChatDB.use <- subsetDB(CellChatDB.mouse, search = "Secreted Signaling", key = "annotation")

# Run one condition
run_cellchat <- function(condition_value, condition_label) {
  message("\nRunning CellChat for ", condition_label)

  condition_obj <- subset(obj, cells = colnames(obj)[obj@meta.data[[condition_col]] == condition_value])
  condition_obj$cellchat_group <- droplevels(condition_obj$cellchat_group)

  cellchat <- createCellChat(condition_obj, group.by = "cellchat_group", assay = "RNA")
  cellchat@DB <- CellChatDB.use

  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, population.size = FALSE, nboot = 100)
  cellchat <- filterCommunication(cellchat, min.cells = min_cells)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)

  saveRDS(cellchat, file.path(output_dir, paste0("cellchat_", condition_label, ".rds")))
  cellchat
}

cellchat_wt <- run_cellchat(conditions[["WT"]], "WT")
cellchat_ko <- run_cellchat(conditions[["KO"]], "KO")


cellchat_wt <- readRDS("/iss-scratch/CoreBioinformatics/rk720/nichenet_stuff/cellchat_eae_ctrl_vs_ko_v2/cellchat_WT.rds")
cellchat_ko <- readRDS("/iss-scratch/CoreBioinformatics/rk720/nichenet_stuff/cellchat_eae_ctrl_vs_ko_v2/cellchat_KO.rds")

cellchat_list <- list(WT = cellchat_wt, KO = cellchat_ko)

# Export aggregated networks
for (condition in names(cellchat_list)) {
  write_csv(
    as.data.frame(cellchat_list[[condition]]@net$count) %>% tibble::rownames_to_column("sender"),
    file.path(output_dir, paste0(condition, "_interaction_counts.csv"))
  )

  write_csv(
    as.data.frame(cellchat_list[[condition]]@net$weight) %>% tibble::rownames_to_column("sender"),
    file.path(output_dir, paste0(condition, "_interaction_strengths.csv"))
  )
}


# Prepare matched WT and KO matrices
pad_network <- function(net, sectors) {
  padded <- matrix(0, nrow = length(sectors), ncol = length(sectors), dimnames = list(sectors, sectors))
  padded[rownames(net), colnames(net)] <- net
  padded
}

all_sectors <- unique(c(
  rownames(cellchat_wt@net$count), colnames(cellchat_wt@net$count),
  rownames(cellchat_ko@net$count), colnames(cellchat_ko@net$count)
))

count_networks <- list(
  WT = pad_network(cellchat_wt@net$count, all_sectors),
  KO = pad_network(cellchat_ko@net$count, all_sectors)
)
count_networks
weight_networks <- list(
  WT = pad_network(cellchat_wt@net$weight, all_sectors),
  KO = pad_network(cellchat_ko@net$weight, all_sectors)
)

cell_colours <- CellChat::scPalette(length(all_sectors))
names(cell_colours) <- all_sectors


#### Plots ####

# Circos: number of interactions
pdf(file.path(output_dir, "cellchat_circos_interaction_counts.pdf"), width = 16, height = 8)
par(mfrow = c(1, 2), xpd = TRUE)

for (condition in names(count_networks)) {
  netVisual_chord_cell(
    cellchat_list[[condition]],
    net = count_networks[[condition]],
    color.use = cell_colours,
    cell.order = all_sectors,
    scale = FALSE,
    directional = 1,
    remove.isolate = FALSE,
    transparency = 0.4,
    lab.cex = 0.8,
    title.name = paste0(condition, ": number of interactions")
  )
}

dev.off()

# Circos: interaction strength
pdf(file.path(output_dir, "cellchat_circos_interaction_strength.pdf"), width = 16, height = 8)
par(mfrow = c(1, 2), xpd = TRUE)

for (condition in names(weight_networks)) {
  netVisual_chord_cell(
    cellchat_list[[condition]],
    net = weight_networks[[condition]],
    color.use = cell_colours,
    cell.order = all_sectors,
    scale = FALSE,
    directional = 1,
    remove.isolate = FALSE,
    transparency = 0.4,
    lab.cex = 0.8,
    title.name = paste0(condition, ": interaction strength")
  )
}

dev.off()

# Circle plot: number of interactions
max_count <- max(unlist(count_networks), na.rm = TRUE)

pdf(file.path(output_dir, "cellchat_circle_interaction_counts.pdf"), width = 16, height = 8)
par(mfrow = c(1, 2), xpd = TRUE, mar = c(1, 1, 4, 1))

for (condition in names(count_networks)) {
  netVisual_circle(
    count_networks[[condition]],
    color.use = cell_colours,
    vertex.weight = rep(1, length(all_sectors)),
    vertex.weight.max = 1,
    vertex.size.max = 12,
    weight.scale = TRUE,
    edge.weight.max = max_count,
    edge.width.max = 10,
    label.edge = FALSE,
    vertex.label.cex = 0.8,
    title.name = paste0(condition, ": number of interactions")
  )
}

dev.off()

# Circle plot: interaction strength
max_weight <- max(unlist(weight_networks), na.rm = TRUE)

pdf(file.path(output_dir, "cellchat_circle_interaction_strength.pdf"), width = 16, height = 8)
par(mfrow = c(1, 2), xpd = TRUE, mar = c(1, 1, 4, 1))

for (condition in names(weight_networks)) {
  netVisual_circle(
    weight_networks[[condition]],
    color.use = cell_colours,
    vertex.weight = rep(1, length(all_sectors)),
    vertex.weight.max = 1,
    vertex.size.max = 12,
    weight.scale = TRUE,
    edge.weight.max = max_weight,
    edge.width.max = 10,
    label.edge = FALSE,
    vertex.label.cex = 0.8,
    title.name = paste0(condition, ": interaction strength")
  )
}

dev.off()


# Check contribution in CXCL signalling
pathway <- "CXCL"

p_wt <- netAnalysis_contribution(cellchat_wt, signaling = pathway, title = "WT: CXCL signalling")
p_ko <- netAnalysis_contribution(cellchat_ko, signaling = pathway, title = "KO: CXCL signalling")

ggsave(
  file.path(output_dir, "CXCL_ligand_receptor_contributions.pdf"),
  p_wt + p_ko,
  width = 10,
  height = 4
)

cxcl_wt <- netAnalysis_contribution(
  cellchat_wt,
  signaling = pathway,
  title = "WT: CXCL signalling",
  return.data = TRUE
)

cxcl_ko <- netAnalysis_contribution(
  cellchat_ko,
  signaling = pathway,
  title = "KO: CXCL signalling",
  return.data = TRUE
)

cxcl_contributions <- bind_rows(
  cxcl_wt$LR.contribution %>% mutate(condition = "WT"),
  cxcl_ko$LR.contribution %>% mutate(condition = "KO")
) %>%
  select(condition, everything()) %>%
  group_by(condition) %>%
  arrange(desc(contribution), .by_group = TRUE) %>%
  mutate(rank = row_number()) %>%
  ungroup() %>%
  relocate(condition, rank)

View(cxcl_contributions)

write_csv(
  cxcl_contributions,
  file.path(output_dir, "CXCL_ligand_receptor_contributions.csv")
)



# Top interactions plots
keep_top_edges <- function(networks, n = 10) {
  combined <- networks$WT + networks$KO

  edge_table <- as.data.frame(as.table(combined)) %>%
    dplyr::rename(sender = Var1, receiver = Var2, combined_count = Freq) %>%
    dplyr::filter(combined_count > 0) %>%
    dplyr::arrange(dplyr::desc(combined_count)) %>%
    dplyr::slice_head(n = n)

  mask <- matrix(
    FALSE,
    nrow = nrow(combined),
    ncol = ncol(combined),
    dimnames = dimnames(combined)
  )

  mask[cbind(
    match(edge_table$sender, rownames(mask)),
    match(edge_table$receiver, colnames(mask))
  )] <- TRUE

  filtered_networks <- lapply(networks, function(net) {
    filtered <- net
    filtered[!mask] <- 0
    filtered
  })

  list(networks = filtered_networks, edges = edge_table)
}

plot_top_networks <- function(networks, n, filename) {
  max_count <- max(unlist(networks), na.rm = TRUE)

  pdf(file.path(output_dir, filename), width = 16, height = 8)
  par(mfrow = c(1, 2), xpd = TRUE, mar = c(1, 1, 4, 1))

  for (condition in names(networks)) {
    netVisual_circle(
      networks[[condition]],
      color.use = cell_colours,
      vertex.weight = rep(1, length(all_sectors)),
      vertex.weight.max = 1,
      vertex.size.max = 10,
      weight.scale = TRUE,
      edge.weight.max = max_count,
      edge.width.max = 15,
      label.edge = FALSE,
      vertex.label.cex = 0.85,
      remove.isolate = FALSE,
      title.name = paste0(condition, ": top ", n, " interactions")
    )
  }

  dev.off()
}

for (n_int in c(15, 20, 25, 30)) {
  top_edges <- keep_top_edges(count_networks, n = n_int)

  plot_top_networks(
    top_edges$networks,
    n = n_int,
    filename = paste0("cellchat_circle_top", n_int, "_edges_interaction_counts.pdf")
  )
}