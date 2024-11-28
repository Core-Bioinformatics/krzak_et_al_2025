library(ggplot2)
library(Seurat)
library(harmony)

metadata = read.csv('metadata_cellranger.csv')

started = F
for (row_ix in seq(nrow(metadata))) {
  print(metadata$Barcode[row_ix])
  sample.so = Read10X_h5(paste('cellranger',
                               metadata$Barcode[row_ix],
                               'outs/filtered_feature_bc_matrix.h5',
                               sep='/'))
  sample.so = CreateSeuratObject(sample.so, project = metadata$Barcode[row_ix])

  sample.so$pool = metadata$Pool[row_ix]
  sample.so$barcode = metadata$Barcode[row_ix]
  sample.so$eae = metadata$eae[row_ix]
  sample.so$rep = metadata$rep[row_ix]
  sample.so$condition = metadata$condition[row_ix]
  sample.so$sample_name = metadata$Sample.name[row_ix]
    
  if (started == F) {
    so = sample.so
    started = T
  } else {
    so = merge(so, y = c(sample.so), project = "Krzak")
  }
}

mt_genes = grep('^mt-', row.names(so), value=TRUE)
rp_genes = grep('^Rp[sl]', row.names(so), value=TRUE)

so$percent_mt = PercentageFeatureSet(so, features=mt_genes)
so$percent_rp = PercentageFeatureSet(so, features=rp_genes)

so_subset = subset(so,(
  nFeature_RNA < 3000 &
  nCount_RNA < 3000 & 
  percent_rp < 40 &
  percent_mt < 30))

remove = read.csv('processed/remove_hb.txt')$x

so_subset = subset(so_subset, cells=colnames(so_subset)[!(colnames(so_subset) %in% remove)])

so_subset@meta.data$rep = as.character(so_subset@meta.data$rep)
so_subset <- SCTransform(so_subset, assay = "RNA")
so_subset <- RunPCA(so_subset, assay = "SCT", verbose = T)
so_subset= RunHarmony(so_subset, group.by.vars = 'rep')
so_subset <- RunUMAP(so_subset, reduction = 'harmony', dims = 1:30, umap.method='uwot')

so_subset <- FindNeighbors(so_subset, reduction = "harmony")
so_subset <- FindClusters(so_subset, resolution = 0.4)

library(ShinyCell)

scConf = createConfig(so_subset)
makeShinyApp(so_subset, scConf, gene.mapping = TRUE,
             shiny.title = "ShinyCell", shiny.dir = "processed/shinyapp/") 
