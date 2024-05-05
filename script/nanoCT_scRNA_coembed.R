library(Signac)
library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(EnsDb.Hsapiens.v86)
library(biovizBase)

custom_palette <- c(
  "#1f77b4", "#ff7f0e",  "#2ca02c", "#d62728", "#9467bd", "#bcbd22", "#8c564b",
 "#7f7f7f", "#17becf",  "#aec7e8", "#ffbb78","#98df8a",  "#ff9896", "#9467bd", "#c5b0d5"
)

filter = readRDS("human_first_updated.rds")

head(filter)

filter 

DefaultAssay(filter) <- "RNA"

filter <- NormalizeData(filter)

filter <- FindVariableFeatures(filter)

filter <- ScaleData(filter)

filter <- RunPCA(filter)

filter <- RunUMAP(filter, dims = 1:30, reduction = "pca")

chromatin <- readRDS("/cfs/klemming/projects/supr/snic2022-23-547/yiquan/nanoCT/output/H3K27ac_bin.rds")

DefaultAssay(chromatin) <- "bins"

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"
Annotation(chromatin) <- annotations

gene.activities <- GeneActivity(chromatin, features = VariableFeatures(filter))

chromatin[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)

# normalize gene activities
DefaultAssay(chromatin) <- "ACTIVITY"

filter

chromatin

transfer.anchors <- FindTransferAnchors(
  reference = filter,
  query = chromatin,
  reduction = 'cca',
  features = VariableFeatures(filter),
  reference.assay = "RNA", 
  query.assay = "ACTIVITY",
  k.filter = NA
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = filter$cell_type,
  weight.reduction = chromatin[['lsi']],
  dims = 2:30
)

chromatin <- AddMetaData(object = chromatin, metadata = predicted.labels)

write.table(predicted.labels, "predicted_labels_H3K27ac.tsv", sep = "\t", quote = FALSE)

chromatin$cell_type <- chromatin$predicted.id

plot1 <- DimPlot(
  object = filter,
  group.by = 'cell_type',
  label = TRUE,
  raster= FALSE,
  repel = TRUE,
  pt.size = 0.01) + ggtitle('scRNA-seq')

plot2 <- DimPlot(
  object = chromatin,
  group.by = 'predicted.id',
  label = TRUE,
  raster= FALSE,
  repel = TRUE,
  pt.size = 0.05) + ggtitle('scCUT&Tag')

p3 <- plot1 |  plot2

ggsave('coplot_H3K27ac_1.png', p3,  width = 28, height = 12, dpi = 300)

genes.use <- VariableFeatures(filter)

refdata <- GetAssayData(filter, assay = "RNA", layer = "data")[genes.use, ]

# refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = chromatin[["lsi"]], dims = 2:30)

chromatin[["RNA"]] <- imputation

DefaultAssay(filter) <- "RNA"

filter$modality <- "RNA"
chromatin$modality <- "CUT&Tag"

coembed <- merge(x = filter, y = chromatin)

DefaultAssay(coembed) <- "RNA"

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)

coembed 

rna <- coembed[, coembed$modality %in% c("RNA")]

ct <- coembed[, coembed$modality %in% c("CUT&Tag")]

rm(coembed)

p1 <- DimPlot(rna, cols = custom_palette, group.by = "cell_type", raster = FALSE, pt.size= 0.01, label = TRUE) + ggtitle('scRNA-seq')

p2 <- DimPlot(ct, cols = custom_palette, group.by = "cell_type", raster = FALSE, pt.size= 0.5, label = TRUE) + ggtitle('nanoCT')

p3 <- p1|p2

ggsave("coembed_H3K27ac_bin_coplot_1.png", p3, width = 30, height = 14, dpi = 300)

p4 <- DimPlot(ct, pt.size= 0.5, label = TRUE)

ggsave("coembed_H3K27ac_bin_cluster_1.png", p4, width = 11, height = 10, dpi = 300)
