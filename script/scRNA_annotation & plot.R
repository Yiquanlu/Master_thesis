library(Seurat)
library(patchwork)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(dplyr)

integrate <- readRDS("human_RNA_all.rds")

table(integrate$age)
head(integrate)

integrate <- RunUMAP(integrate, reduction = "pca", dims = 1:30, min.dist = 0.5, n.neighbors=50L)
integrtae <- FindNeighbors(integrate, reduction = "pca", dims = 1:30)
integrate <- FindClusters(integrate, resolution = 0.3, graph.name = "integrated_snn")

metadata <- integrate@meta.data

metadata$age <- as.character(metadata$age)

age_annotations <- c('10' = '10',
			'11.5' = '11.5',
			'12' = '12',
			'13' = '13',
			'14' = '14',
			'35 GW+12 d' = '34.7',
			'39 GW+1d' = '37.1',
			'6.69999980926514' = '6.7',
			'6.90000009536743' = '6.9',
			'9.19999980926514' = '9.2',
			'8' = '8',
			'9.5' = '9.5',
			'ga22' = '20',
			'ga24' = '22',
			'ga34' = '32',
			'GW16' = '14',
			'GW17' = '15',
			'GW18' = '16',
			'GW19' = '17',
                        'GW20' = '18',
			'GW20.5' = '18.5',
                        'GW21' = '19',
                        'GW22' = '20',
                        'GW23' = '21',
                        'GW24' = '22',
                        'GW25' = '23',
			'GW25+1d' = '23.1',
                        'GW27+2d' = '25.3',
                        'GW28' = '26',
                        'GW28.5' = '26.5',
			'GW30' = '28',
			'GW32' = '30',
			'GW33' = '31',
			'GW34' = '32',
			'GW37+5d' = '35.7',
			'GW38+3d' = '38.4',
			'GW41' = '39') 

metadata$pcw <- age_annotations[metadata$age]


integrate <- AddMetaData(integrate, metadata = metadata)

integrate$pcw <- as.numeric(integrate$pcw)

custom_palette <- c(
  "#1f77b4", "#ff7f0e",  "#2ca02c", "#d62728", "#9467bd", "#bcbd22", "#8c564b",
 "#7f7f7f", "#17becf",  "#aec7e8", "#ffbb78","#98df8a",  "#ff9896", "#9467bd", "#c5b0d5"
)

p1 <- DimPlot(integrate, label = TRUE, reduction = "umap", raster = FALSE, group.by = "cell_type", pt.size = 0.01)
p2 <- FeaturePlot(integrate, features = "pcw", raster = FALSE, pt.size= 0.1) + scale_colour_viridis_c()
p3 <- DimPlot(integrate, cols = custom_palette, label = TRUE, reduction = "umap", raster = FALSE, group.by = "lineage", pt.size = 0.01)
p4 <- DimPlot(integrate, label = TRUE, reduction = "umap", raster = FALSE, pt.size = 0.01)

p5 <- grid.arrange(p1, p2, p3, p4, nrow = 2)

ggsave('integration_all_0.3.png', p5,  width = 28, height = 22, dpi = 300)
ggsave('time_all.png', p2, width = 30, height =27, dpi =300)

metadata <- integrate@meta.data

metadata$seurat_clusters <- as.character(metadata$seurat_clusters)

cell_annotations <- c('0' = 'Excitatory neuron',
			'1' = 'Dorsal RG',
			'2' = 'Interneuron',
			'3' = 'Ventral RG',
			'4' = 'Excitatory neuron',
			'5' = 'Excitatory neuron',
			'6' = 'Excitatory neuron',
			'7' = 'Glioblast',
			'8' = 'Interneuron',
			'9' = 'Interneuron',
			'10' = 'Excitatory neuron',
			'11' = 'Ventral RG',
			'12' = 'Excitatory neuron',
			'13' = 'Excitatory neuron',
			'14' = 'Interneuron',
			'15' = 'Excitatory neuron',
			'16' = 'Astrocyte',
			'17' = 'Excitatory neuron',
			'18' = 'Excitatory neuron',
                        '19' = 'Intermediate progenitor cell',
			'20' = 'Excitatory neuron',
                        '21' = 'Intermediate progenitor cell',
                        '22' = 'Excitatory neuron',
                        '23' = 'Dorsal RG',
                        '24' = 'Excitatory neuron',
                        '25' = 'OPC/Oligo',
			'26' = 'Dorsal RG',
                        '27' = 'Excitatory neuron',
                        '28' = 'Microglia',
                        '30' = 'Excitatory neuron',
			'29' = 'Cajal-Retzius cell')

metadata$annotation_0.3 <- cell_annotations[metadata$seurat_clusters]

integrate <- AddMetaData(integrate, metadata = metadata)

p7 <- DimPlot(integrate, cols = custom_palette, label = TRUE, reduction = "umap", raster = FALSE, group.by = "annotation_0.3", pt.size = 0.01)

ggsave('annotation_integrate_0.3.png', p7,  width = 30, height = 27, dpi = 300)

saveRDS(integrate, "human_RNA_all_0.3.rds")

all.markers <- FindAllMarkers(object = integrate, only.pos = TRUE)

write.table(all.markers, "markers_0.3.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
