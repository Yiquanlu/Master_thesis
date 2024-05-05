library(Signac)
library(Seurat)
library(ggplot2)

brain <- readRDS("/cfs/klemming/projects/snic/snic2022-23-547/yiquan/out/gliogenesis.rds")

features1 <- c('HES1', 'FOXG1', 'NKX2-1', 'EN1', 'EN2', 'TNC', 'VIM', 'SOX2', 'SLC1A3', 'EGFR', 'DLX1', 'HOPX','SPARCL1', 'SLC32A1', 'GFAP', 'ALDH1L1', 'PDGFRA', 'CSPG4', 'OLIG1', 'OLIG2', 'FGF2', 'SOX10', 'AQP4')

DefaultAssay(brain) <- 'RNA'

p1 <- DotPlot(object = brain, features = features1, dot.scale = 8, group.by = "annotation")

ggsave('dot_marker.png', p1, width = 20, height = 7, dpi = 300)


brain <- readRDS("/cfs/klemming/projects/snic/snic2022-23-547/yiquan/second_tri/human_RNA_all_0.3.rds")

DefaultAssay(brain) <- 'integrated'

features2 <- c('HES1', 'VIM', 'PAX6', 'SLC1A3', 'DLX2', 'GFAP', 'AQP4', 'OLIG1', 'OLIG2', 'SOX10', 'PDGFRA', 'P2RY12', 'EOMES', 'SLC6A1', 'SATB2', 'GAD2', 'LHX6', 'EBF3')

p1 <- DotPlot(object = brain, features = features2, dot.scale = 8, group.by = "annotation_0.3")

ggsave('dot_marker2.png', p1, width = 18, height = 6, dpi = 300)


brain <- readRDS("/cfs/klemming/projects/snic/snic2022-23-547/yiquan/out/signac_integrate5.rds")

DefaultAssay(brain) <- "ATAC"

metadata <- brain@meta.data

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
			'18' = 'Ventral RG',
			'19' = 'Interneuron',
			'20' = 'Intermediate progenitor cells',
			'21' = 'Excitatory neuron',
			'22' = 'Dorsal RG',
			'23' = 'Excitatory neuron', 
			'24' = 'Excitatory neuron',
			'25' = 'Dorsal RG',
			'26' = 'Excitatory neuron',
			'27' = 'OPC/Oligo',
			'28' = 'Dorsal RG',
			'29' = 'Excitatory neuron',
			'30' = 'Microglia',
			'31' = 'Dorsal RG',
			'32' = 'Interneuron',
			'33' = 'Pericytes',
			'34' = 'Excitatory neuron', 
			'35' = 'Endothelial cell')

metadata$cell <- cell_annotations[metadata$ATAC_snn_res.0.7]

brain <- AddMetaData(brain, metadata = metadata)

DefaultAssay(brain) <- 'RNA'

features3 <- c('SOX2', 'HES1', 'VIM', 'PAX6', 'NKX2-1', 'GFAP', 'OLIG2', 'SOX10', 'P2RY12', 'EOMES', 'TBR1', 'NEUROD1', 'SLC17A7', 'SATB2', 'SLC32A1', 'GAD2', 'LHX6', 'SP8')

p1 <- DotPlot(object = brain, features = features3, dot.scale = 8, group.by = "cell")

ggsave('dot_marker3.png', p1, width = 18, height = 6, dpi = 300)
