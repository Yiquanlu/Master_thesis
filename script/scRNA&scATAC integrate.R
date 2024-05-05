library(Seurat)
library(Signac)
library(ggplot2)
library(cowplot)
library(patchwork)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(dplyr)
library(EnsDb.Hsapiens.v86)
library(biovizBase)

filter <- readRDS('/cfs/klemming/projects/supr/snic2022-23-547/yiquan/second_tri/human_RNA_all_0.3.rds')


set.seed(42)  # for reproducibility


custom_palette <- c(
  "#1f77b4", "#ff7f0e",  "#2ca02c", "#d62728", "#9467bd", "#bcbd22", "#8c564b",
 "#7f7f7f", "#17becf",  "#aec7e8", "#ffbb78","#98df8a",  "#ff9896", "#9467bd", "#c5b0d5"
)

meta <- filter@meta.data

DefaultAssay(filter) <- "integrated"

filter

chromatin <- readRDS("/cfs/klemming/projects/supr/snic2022-23-547/yiquan/out/signac_integrate5.rds")

DefaultAssay(chromatin) <- "ATAC"

chromatin <- RunUMAP(object = chromatin, reduction = 'lsi', dims = 2:30)

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"
Annotation(chromatin) <- annotations

gene.activities <- GeneActivity(chromatin, features = VariableFeatures(filter))

chromatin[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)

# normalize gene activities
DefaultAssay(chromatin) <- "ACTIVITY"

chromatin_features <- rownames(chromatin)

chromatin <- NormalizeData(chromatin)

chromatin <- ScaleData(chromatin, features = chromatin_features)

chromatin

transfer.anchors <- FindTransferAnchors(
  reference = filter,
  query = chromatin,
  reduction = 'cca',
  features = chromatin_features,
  reference.assay = "integrated", 
  query.assay = "ACTIVITY"
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = filter$annotation_0.3,
  weight.reduction = chromatin[['integrated_lsi']],
  dims = 2:30
)

chromatin <- AddMetaData(object = chromatin, metadata = predicted.labels)

write.table(predicted.labels, "predicted_labels.tsv", sep = "\t", quote = FALSE)

plot1 <- DimPlot(
  object = filter,
  group.by = 'annotation_0.3',
  label = TRUE,
  raster= FALSE,
  pt.size = 0.01) + ggtitle('scRNA-seq')

plot2 <- DimPlot(
  object = chromatin,
  group.by = 'predicted.id',
  label = TRUE,
  raster= FALSE,
  pt.size = 0.05) + ggtitle('scATAC-seq')

p3 <- plot1 | plot2

ggsave('coplot_atac.png', p3,  width = 35, height = 16, dpi = 300)

genes.use <- VariableFeatures(filter)

refdata <- GetAssayData(filter, assay = "RNA", layer = "data")[genes.use, ]

# refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells

imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = chromatin[["integrated_lsi"]], dims = 2:30)

chromatin[["integrated"]] <- imputation

chromatin

DefaultAssay(chromatin) <- "integrated"
DefaultAssay(filter) <- "integrated"

coembed <- merge(x = filter, y = chromatin)

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets

coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30, min.dist = 0.5, n.neighbors = 50L)

coembed

DefaultAssay(coembed) <- "integrated"

table(coembed$dataset)

metadata <- coembed@meta.data

metadata$dataset <- as.character(metadata$dataset)

modality_annotations <- c('braun' = 'RNA',
			'velme' = 'RNA',
			'DC' = 'ATAC',
			'advance' = 'ATAC',
			'biox' = 'ATAC',
			'cell' = 'ATAC',
			'nature' = 'ATAC')


metadata$mod <- modality_annotations[metadata$dataset]

coembed <- AddMetaData(object = coembed, metadata = metadata)

table(coembed$mod)

p1 <-DimPlot(coembed, cols = custom_palette, group.by = "mod", pt.size = 0.01, raster = FALSE)

ggsave("coembed_sc_atac.png", p1, width = 30, height = 27, dpi = 300)

p2 <- DimPlot(coembed, cols = custom_palette, group.by = "annotation_0.3", raster = FALSE, pt.size= 0.01, label = TRUE)

ggsave("coembed_cell_atac.png", p2, width = 30, height = 27, dpi = 300)

plot <- p1 + p2 + plot_layout(ncol=2)

ggsave("coembed_plot_atac.png", plot, width = 30, height = 14, dpi = 300)

rna <- coembed[, coembed$mod %in% c("RNA")]

atac <- coembed[, coembed$mod %in% c("ATAC")]

rm(coembed)

p1 <- DimPlot(rna, cols = custom_palette, group.by = "annotation_0.3", raster = FALSE, pt.size= 0.01, label = TRUE) + ggtitle('scRNA-seq')

p2 <- DimPlot(atac, cols = custom_palette, group.by = "annotation_0.3", raster = FALSE, pt.size= 0.05, label = TRUE) + ggtitle('scATAC-seq')

p3 <- p1|p2

ggsave("coembed_mod_atac.png", p3, width = 30, height = 14, dpi = 300)

