library(Signac)
library(Seurat)
library(chromVAR)
library(motifmatchr)
library(JASPAR2020)
library(ggseqlogo)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)

# you can do the same step for interneuron
integrate <- readRDS("/cfs/klemming/projects/supr/snic2022-23-547/yiquan/out/gliogenesis.rds")

DefaultAssay(integrate) <- "ATAC"

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
integrate <- AddMotifs(
  object = integrate,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

integrate <- RunChromVAR(
  object = integrate,
  genome = BSgenome.Hsapiens.UCSC.hg38
)


saveRDS(integrate, "gliogenesis.rds")

head(integrate)

integrate <- RegroupIdents(integrate, metadata = "annotation")

table(integrate$annotation)

levels(integrate)

cell_idents <- Idents(integrate)

table(cell_idents)

DefaultAssay(integrate) <- 'chromvar'

differential.activity1 <- FindMarkers(
  object = integrate,
  ident.1 = 'GPC',
  ident.2 = 'Radial glia',
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

differential.activity2 <- FindMarkers(
  object = integrate,
  ident.1 = 'Astrocyte',
  ident.2 = 'OPCs',
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

differential.activity3 <- FindMarkers(
  object = integrate,
  ident.1 = 'Astrocyte',
  ident.2 = 'GPC',
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

differential.activity4 <- FindMarkers(
  object = integrate,
  ident.1 = 'OPCs',
  ident.2 = 'GPC',
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

write.table(differential.activity1, file = "GPC_RG_chromvar.tsv", sep = "\t", quote = FALSE)
write.table(differential.activity2, file = "astro_opc_chromvar.tsv", sep = "\t", quote = FALSE)
write.table(differential.activity3, file = "astro_GPC_chromvar.tsv", sep = "\t", quote = FALSE)
write.table(differential.activity4, file = "opc_GPC_chromvar.tsv", sep = "\t", quote = FALSE)

DefaultAssay <- 'chromvar'

p1 <- FeaturePlot(
  object = integrate,
  features = c('MA0670.1', 'MA1643.1', 'MA0161.2', 'MA1116.1', 'MA0137.3', 'MA0144.2', 'MA0729.1', 'MA1556.1', 'MA0017.2', 'MA0139.1', 'MA1150.1', 'MA1649.1', 'MA1587.1', 'MA1125.1', 'MA0099.3', 'MA0488.1', 'MA0476.1', 'MA0143.4'),
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1,
  ncol = 3,
  raster = FALSE
)

DefaultAssay <- 'ATAC'

p2 <- MotifPlot(
  object = integrate,
  motifs = c('MA0670.1', 'MA1643.1', 'MA0161.2', 'MA0671.1', 'MA1116.1', 'MA0137.3', 'MA0144.2', 'MA0729.1', 'MA1556.1', 'MA0017.2', 'MA0139.1', 'MA1150.1', 'MA1649.1', 'MA1587.1', 'MA1125.1'),
  assay = 'ATAC'
)

ggsave("glia_motif_enrichment.png", p1, width = 22, height= 36, dpi = 300)
ggsave("glia_motifs.png", p1, width = 22, height= 36, dpi = 300)