library(Signac)
library(Seurat)
library(RColorBrewer)
library(dplyr)
library(ggplot2)

custom_palette <- c(
  "#1f77b4", "#ff7f0e",  "#2ca02c", "#d62728", "#9467bd", "#bcbd22", "#8c564b", 
 "#7f7f7f", "#17becf",  "#aec7e8", "#ffbb78","#98df8a",  "#ff9896", "#9467bd", "#c5b0d5"
)

Glia <- readRDS("/cfs/klemming/projects/snic/snic2022-23-547/yiquan/out/interneurongenesis.rds")

glia <- RunUMAP(glia, reduction = "integrated_lsi", dims = 2:40)
glia <- FindNeighbors(object = glia, reduction = 'integrated_lsi', dims = 2:40)
glia <- FindClusters(object = glia, verbose = FALSE, algorithm = 3, resolution = 0.4, graph.name = "ATAC_snn")

p4 <- DimPlot(object = glia, label = TRUE)

p5 <- DimPlot(object = glia, cols = custom_palette, label = TRUE, group.by ="cell", pt.size= 0.1) 

p6 <- DimPlot(object = glia, label = TRUE, group.by ="cell_type", pt.size= 0.1)

p7 <- FeaturePlot(object = glia, features = "Age", raster = FALSE, pt.size= 0.1) + scale_colour_viridis_c()

p8 <- p5 | p7 | p6

ggsave('cluster_intern.png', p4,  width = 10, height = 8, dpi = 300)
ggsave('intern.png', p8,  width = 35, height = 7, dpi = 300)

metadata <- glia@meta.data

metadata$ATAC_snn_res.0.4 <- as.character(metadata$ATAC_snn_res.0.4)

cell_annotations <- c('0' = 'CGE interneuron',
			'1' = 'MGE interneuron',
			'2' = 'CGE interneuron',
			'3' = 'MGE interneuron',
			'4' = 'MGE progenitors',
			'5' = 'CGE progenitors',
			'6' = 'Radial glia',
			'7' = 'LGE interneuron',
			'8' = 'CGE interneuron',
                        '9' = 'LGE interneuron',
                        '10' = 'LGE interneuron')

metadata$annotation <- cell_annotations[metadata$ATAC_snn_res.0.4]

glia <- AddMetaData(glia, metadata = metadata)

p5 <- DimPlot(object = glia, cols = custom_palette, label = TRUE, group.by ="annotation", pt.size= 0.1)
ggsave('interneuron_annotation.png', p5,  width = 10, height = 8, dpi = 300)

saveRDS(glia, "interneurongenesis.rds")
