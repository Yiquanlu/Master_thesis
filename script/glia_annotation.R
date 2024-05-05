library(Signac)
library(Seurat)
library(RColorBrewer)
library(dplyr)
library(ggplot2)

custom_palette <- c(
  "#1f77b4", "#ff7f0e",  "#2ca02c", "#d62728", "#9467bd", "#bcbd22", "#8c564b", 
 "#7f7f7f", "#17becf",  "#aec7e8", "#ffbb78","#98df8a",  "#ff9896", "#9467bd", "#c5b0d5"
)

glia <- readRDS("/cfs/klemming/projects/snic/snic2022-23-547/yiquan/gliogenesis.rds")

glia <- FindNeighbors(object = glia, reduction = 'integrated_lsi', dims = 2:30)

glia <- FindClusters(object = glia, verbose = FALSE, algorithm = 3, resolution = 0.4, graph.name = "ATAC_snn")

marker <- FindAllMarkers(glia, only.pos = TRUE)

write.table(marker, file = "gliogenesis_markers.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

p4 <- DimPlot(object = glia, label = TRUE)

p5 <- DimPlot(object = glia, cols = custom_palette, label = TRUE, group.by ="cell", pt.size= 0.1) 

p6 <- DimPlot(object = glia, label = TRUE, group.by ="cell_type", pt.size= 0.1)

p7 <- FeaturePlot(object = glia, features = "Age", raster = FALSE, pt.size= 0.1) + scale_colour_viridis_c()

p8 <- p5 | p7 | p6

ggsave('cluster_glia.png', p4,  width = 10, height = 8, dpi = 300)
ggsave('glia.png', p8,  width = 35, height = 7, dpi = 300)

metadata <- glia@meta.data

metadata$ATAC_snn_res.0.4 <- as.character(metadata$ATAC_snn_res.0.4)

cell_annotations <- c('0' = 'GPC',
			'1' = 'Astrocyte',
			'2' = 'OPCs',
			'3' = 'Astrocyte',
			'4' = 'OPCs',
			'5' = 'Radial glia',
			'6' = 'GPC',
			'7' = 'OPCs',
			'8' = 'OPCs')

metadata$annotation <- cell_annotations[metadata$ATAC_snn_res.0.4]

glia <- AddMetaData(glia, metadata = metadata)

p5 <- DimPlot(object = glia, cols = custom_palette, label = TRUE, group.by ="annotation", pt.size= 0.1)
ggsave('glia_annotation.png', p5,  width = 10, height = 8, dpi = 300)
