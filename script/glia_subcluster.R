library(Signac)
library(Seurat)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(patchwork)

integrated <- readRDS("/cfs/klemming/projects/snic/snic2022-23-547/yiquan/out/signac_integrate5.rds")

DefaultAssay(integrated) <- "ATAC"

metadata <- integrated@meta.data

metadata$ATAC_snn_res.0.7 <- as.character(metadata$ATAC_snn_res.0.7)

cell_annotations <- c('0' = 'Intermediate progenitor cells',
			'1' = 'Dorsal RG',
			'2' = 'Interneuron',
			'3' = 'Interneuron',
			'4' = 'Excitatory neuron',
			'5' = 'Excitatory neuron',
			'6' = 'Interneuron',
			'7' = 'Ventral RG',
			'8' = 'Excitatory neuron',
			'9' = 'Excitatory neuron',
			'10' = 'Excitatory neuron',
			'11' = 'Unknown',
			'12' = 'Excitatory neuron',
			'13' = 'Astrocyte',
			'14' = 'Interneuron',
			'15' = 'Excitatory neuron',
			'16' = 'Excitatory neuron',
			'17' = 'Excitatory neuron',
			'18' = 'Midbrain RG',
			'19' = 'Interneuron',
			'20' = 'Intermediate progenitor cells',
			'21' = 'Excitatory neuron',
			'22' = 'Telen RG',
			'23' = 'Excitatory neuron', 
			'24' = 'Excitatory neuron',
			'25' = 'Dorsal RG',
			'26' = 'Excitatory neuron',
			'27' = 'OPC/Oligo',
			'28' = 'Telen RG',
			'29' = 'Excitatory neuron',
			'30' = 'Microglia',
			'31' = 'Telen RG',
			'32' = 'Interneuron',
			'33' = 'Pericytes',
			'34' = 'Excitatory neuron', 
			'35' = 'Endothelial cell')

metadata$cell <- cell_annotations[metadata$ATAC_snn_res.0.7]

integrated <- AddMetaData(integrated, metadata = metadata)

integrated$Age <- as.numeric(gsub("pcw", "", integrated$Age))
custom_breaks <- seq(7, 24, by = 6)

custom_palette <- c(
  "#1f77b4", "#ff7f0e",  "#2ca02c", "#d62728", "#9467bd", "#bcbd22", "#8c564b", 
 "#7f7f7f", "#17becf",  "#aec7e8", "#ffbb78","#98df8a",  "#ff9896", "#9467bd", "#c5b0d5"
)

p2 <- DimPlot(object = integrated, cols = custom_palette, label = TRUE, raster = FALSE, group.by ="cell", pt.size= 0.01)
p3 <- FeaturePlot (object = integrated, features = "Age", raster = FALSE, pt.size= 0.01) + scale_colour_viridis_c(breaks=custom_breaks)

p1 <- p2 | p3

ggsave("integrate5.png", p2, width =20, height =18, dpi = 300)

glia  <- integrated[, integrated$cell %in% c("Astrocyte", "OPC/Oligo")]

rm(integrated)

DefaultAssay(glia) <- "ATAC"

glia <- RunUMAP(glia, reduction = "integrated_lsi", dims = 2:30)

glia <- FindNeighbors(object = glia, reduction = 'integrated_lsi', dims = 2:30)

glia <- FindClusters(object = glia, verbose = FALSE, algorithm = 3, resolution = 0.4, graph.name = "ATAC_snn")

p4 <- DimPlot(object = glia, label = TRUE)

p5 <- DimPlot(object = glia, cols = custom_palette, label = TRUE, group.by ="cell", pt.size= 0.1)

p6 <- DimPlot(object = glia, label = TRUE, group.by ="cell_type", pt.size= 0.1)

p7 <- FeaturePlot(object = glia, features = "Age", raster = FALSE, pt.size= 0.1) + scale_colour_viridis_c()

p8 <- p5 | p7 | p6

ggsave('cluster_glia.png', p4,  width = 10, height = 8, dpi = 300)
ggsave('glia.png', p8,  width = 35, height = 7, dpi = 300)
saveRDS(glia, "gliogenesis.rds")

DefaultAssay(glia) <- "RNA"

marker <- FindAllMarkers(glia, only.pos = TRUE)

write.table(marker, file = "gliogenesis_markers.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
