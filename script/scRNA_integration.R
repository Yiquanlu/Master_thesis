library(Seurat)
library(dplyr)
library(SeuratData)
library(patchwork)
library(ggplot2)
library(RColorBrewer)

ifnb <- readRDS("human_RNA_merge_all.rds")

ifnb <- ifnb[, ifnb$lineage %in% c('AST', 'ExNeu', 'GLIALPROG', 'Glioblast', 'Immune',
'IN', 'MG', 'Neuroblast', 'Neuron', 'Neuronal IPC', 'OL', 'Oligo', 'OPC', 'Radial glia')]

ifnb
head(ifnb)
table(ifnb$dataset)
table(ifnb$age_range)

ifnb.list <- SplitObject(ifnb, split.by = "dataset")

# normalize and identify variable features for each dataset independently
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
features <- SelectIntegrationFeatures(object.list = ifnb.list)
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})

immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features, reduction = "rpca", k.anchor = 20)

immune.combined <- IntegrateData(anchorset = immune.anchors)

DefaultAssay(immune.combined) <- "integrated"

plot <- ElbowPlot(object, ndims = 40, reduction = "pca")

ggsave('elbowplot.png', plot,  width = 8, height = 16, dpi = 300)

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)

immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.3)

saveRDS(immune.combined, "human_RNA_all.rds")

custom_palette <- c(
  "#1f77b4", "#ff7f0e",  "#2ca02c", "#d62728", "#9467bd", "#bcbd22", "#8c564b", 
 "#7f7f7f", "#17becf",  "#aec7e8", "#ffbb78","#98df8a",  "#ff9896", "#9467bd", "#c5b0d5", "#ff5733", "#6a1b9a", "#66bb6a")

p1 <- DimPlot(immune.combined, label = TRUE, reduction = "umap", raster = FALSE, pt.size = 0.01)
p2 <- DimPlot(immune.combined, label = TRUE, reduction = "umap", raster = FALSE, group.by = "age_range", pt.size = 0.01)
p3 <- DimPlot(immune.combined, cols = custom_palette, label = TRUE, reduction = "umap", raster = FALSE, group.by = "lineage", pt.size = 0.01)
p4 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, raster = FALSE, pt.size = 0.01, group.by = "cell_type")

plot_grid <- p1 + p2 + p3 + p4 + plot_layout(ncol = 2, nrow = 2)

ggsave('integration_all_rpca.png', plot_grid,  width = 28, height =22, dpi = 300)
