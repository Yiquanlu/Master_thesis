library(Seurat)
library(Signac)
library(biovizBase)
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(future)
library(ggplot2)

set.seed(123)

# read in peak sets
peaks.cell <- read.table(file = "/cfs/klemming/projects/supr/snic2022-23-547/yiquan/integrate/cell_peaks.bed",
        col.names = c("chr", "start", "end")
)

peaks.nature <- read.table(
        file = "/cfs/klemming/projects/supr/snic2022-23-547/yiquan/integrate/nature_peaks.bed", 
        col.names = c("chr", "start", "end")
)

peaks.DC <- read.table(
        file = "/cfs/klemming/projects/supr/snic2022-23-547/yiquan/integrate/DC_peaks.bed", 
        col.names = c("chr", "start", "end")
)

peaks.biox <- read.table(
        file = "/cfs/klemming/projects/supr/snic2022-23-547/yiquan/integrate/biox_peaks.bed", 
        col.names = c("chr", "start", "end")
)

peaks.advance <- read.table(
	file = "/cfs/klemming/projects/supr/snic2022-23-547/yiquan/integrate/advance_peaks.bed",
        col.names = c("chr", "start", "end")
)

# convert to genomic ranges
gr.cell <- makeGRangesFromDataFrame(peaks.cell)
gr.nature <- makeGRangesFromDataFrame(peaks.nature)
gr.DC <- makeGRangesFromDataFrame(peaks.DC)
gr.biox <- makeGRangesFromDataFrame(peaks.biox)
gr.advance <- makeGRangesFromDataFrame(peaks.advance)

# Create a unified set of peaks to quantify in each dataset
combined.peaks <- reduce(x = c(gr.cell, gr.nature, gr.DC, gr.biox, gr.advance))

peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]

# load metadata
md.cell <- read.delim(
  file = "/cfs/klemming/projects/supr/snic2022-23-547/yiquan/Cell/data/atac_cell_metadata.txt",
  stringsAsFactors = FALSE,
  header = TRUE,
  row.names = 1
)

md.nature <- read.table(
  file = "/cfs/klemming/projects/supr/snic2022-23-547/yiquan/Nature/sub_meta.tsv",
  stringsAsFactors = FALSE,
  sep = '\t',
  header = TRUE,
  row.names = 1
)

md.DC <- read.delim(
  file = "/cfs/klemming/projects/supr/snic2022-23-547/yiquan/DC/DC_submeta.tsv",
  stringsAsFactors = FALSE,
  header = TRUE,
  row.names = 1,
  sep = '\t'
)

md.biox <- read.delim(
  file = "/cfs/klemming/projects/supr/snic2022-23-547/yiquan/Biox/submeta.tsv",
  stringsAsFactors = FALSE,
  header = TRUE,
  row.names = 1,
  sep = '\t'
)

md.advance <- read.delim(
  file = "/cfs/klemming/projects/supr/snic2022-23-547/yiquan/Advance/metadata.tsv",
  stringsAsFactors = FALSE,
  header = TRUE,
  row.names = 1,
  sep = '\t'
)

#Assign cell type and sample annotations to the nature metadata
md.nature$cell_type <- md.nature$CellType
md.nature$Sample <- md.nature$sample

#Convert specimen to Age(weeks) from nature paper
md.nature$specimen <- as.character(md.nature$specimen)

nature_age_annotations <- c("GW17" = "pcw15", "GW18" = "pcw16", "GW20" = "pcw18", 
                            "GW21" = "pcw19")

md.nature$Age <- nature_age_annotations[as.character(md.nature$specimen)]
md.nature$class <- md.nature$cell_type
md.nature$class <- as.character(md.nature$class)

nature_class_annotations <- c('ulEN' = 'Excitatory neurons',
    'dlEN' = 'Excitatory neurons',
    'AstroOligo' = 'Astrocytes/Oligo',
    'earlyEN' = 'Excitatory neurons',
    'IPC' = 'neuronal IPC',
    'IN_CGE'= 'Inhibitory neurons',
    'RG' = 'Radial glia',
    'IN_MGE' = 'Inhibitory neurons',
    'MGE_Progenitors' = 'IN_Progenitors',
    'EndoMural' = 'Endothelial',
    'Microglia' = 'Microglia'
)

md.nature$class <- nature_class_annotations[as.character(md.nature$class)]

md.nature$method <- "v1"
# Remove the columns 
md.nature <- md.nature %>%
  select(-cluster, -specimen, -CellType, -area, -total, -duplicate, -sample, -chimeric, -unmapped, -lowmapq, -mitochondrial, -passed_filters, -TSS_fragments, -DNase_sensitive_region_fragments, -enhancer_region_fragments, -promoter_region_fragments, -on_target_fragments, -blacklist_region_fragments, -peak_region_fragments, -logUMI, -promoter_ratio)

# Combine Sample.ID and Cell.Barcode to create row names
md.cell$row_names <- paste0(md.cell$Sample.ID, "_", md.cell$Cell.Barcode)

# Set the row names
rownames(md.cell) <- md.cell$row_names

# Define a mapping of cluster names to cell type annotations
cluster_annotations <- c("c0" = "Glutamatergic neurons", "c1" = "Glutamatergic neurons", 
                         "c2" = "Glutamatergic neurons", "c3" = "GABAergic neurons", "c4" = "GABAergic neurons", 
                         "c5" = "Glutamatergic neurons", "c6" = "Glutamatergic neurons", "c7" = "Glutamatergic neurons", 
                         "c8" = "neuronal IPC", "c9" = "Late RG", "c10" = "mGPC", 
                         "c11" = "Early RG", "c12" = "GABAergic neurons", "c13" = "Glutamatergic neurons", 
                         "c14" = "Glutamatergic neurons", "c15" = "OPC/Oligo", "c16" = "GABAergic neurons", 
                         "c17" = "Pericytes", "c18" = "Glutamatergic neurons", "c19" = "Microglia", 
                         "c20" = "GABAergic neurons", "c21" = "Endothelial cells")

# Assign cell type annotations to the clusters
md.cell$cell_type <- cluster_annotations[md.cell$Iterative.LSI.Clusters]

md.cell$Sample <- md.cell$Sample.ID
md.cell$class <- md.cell$cell_type

md.cell$class <- as.character(md.cell$class)

cell_class_annotations <- c('Glutamatergic neurons' = 'Excitatory neurons',
    'GABAergic neurons' ='Inhibitory neurons', 
    'Late RG'= 'Radial glia',
    'Microglia' = 'Microglia',
    'Pericytes' = 'Pericytes',
    'neuronal IPC' = 'neuronal IPC',
    'Early RG' = 'Radial glia',
    'mGPC' = 'Astrocytes/Oligo',
    'OPC/Oligo' = 'OPC/Oligo',
    'Endothelial cells' = 'Endothelial'
)

md.cell$class <- cell_class_annotations[as.character(md.cell$class)]
md.cell$method <- "multiome_atac"

# Remove the columns 
md.cell <- md.cell %>%
  select(-row_names, -Sample.ID, -Tissue.ID, -Sample.Type, -Batch, -Iterative.LSI.Clusters, -Tss.Enrichment, -.CR.cellQC.barcode, -.CR.cellQC.chimeric, -.CR.cellQC.unmapped, -.CR.cellQC.lowmapq, -.CR.cellQC.mitochondrial, -.CR.cellQC.passed_filters, -.CR.cellQC.cell_id, -.CR.cellQC.is__cell_barcode, -.CR.cellQC.TSS_fragments, -.CR.cellQC.DNase_sensitive_region_fragments, -.CR.cellQC.enhancer_region_fragments, -.CR.cellQC.promoter_region_fragments, -.CR.cellQC.on_target_fragments, -.CR.cellQC.blacklist_region_fragments, -.CR.cellQC.peak_region_fragments, -.CR.cellQC.peak_region_cutsites, -TSS.Enrichment.Unsmoothed, -.CR.cellQC.total, -.CR.cellQC.duplicate, -Cell.Barcode)

#Assign cell type annotations to the clusters
md.DC$cell_type <- md.DC$Cluster

#Convert Sample to Age(weeks) from DC paper
md.DC$Sample <- as.character(md.DC$Sample)

DC_age_annotations <- c("HES1" = "pcw8.5", "P18861_1001" = "pcw9.5", "P18861_1002" = "pcw11")

md.DC$Age <- DC_age_annotations[as.character(md.DC$Sample)]
md.DC$class <- md.DC$cell_type
md.DC$class <- as.character(md.DC$class)

DC_class_annotations <- c('Glioblast/Pre-OPC' = 'Radial glia',
    'Radial Glia/Glioblast' = 'Radial glia',
    'Endothelial' = 'Endothelial',
    'OPCs' = 'OPC/Oligo',
    'VLMCs' = 'Pericytes',
    'Cortical Interneurons' = 'Inhibitory neurons', 
    'Radial Glia cycling' = 'Radial glia',
    'GABAergic forebrain' = 'Inhibitory neurons', 
    'Excitatory neurons cortex' = 'Excitatory neurons',
    'Forebrain early neuroblast possibly GABAergic' = 'IN_Progenitors',
    'Neuroblast motorneuron/GABAergic?' = 'IN_Progenitors',
    'Microglia' = 'Microglia'
)

md.DC$class <- DC_class_annotations[as.character(md.DC$class)]

md.DC$method <- "v1"

# Remove the columns 
md.DC <- md.DC %>%
    select(-Cluster, -nCount_RNA, -nFeature_RNA)

md.biox$Age <- as.numeric(as.character(md.biox$Age))

biox_age_annotations <- c("7" = "pcw7", "8" = "pcw8", "9" = "pcw9", "10" = "pcw10",
                         "12" = "pcw12", "12.5" = "pcw12.5", "13" = "pcw13")


md.biox$Age <- biox_age_annotations[as.character(md.biox$Age)]
md.biox$method <- md.biox$Chemistry

md.biox$cell_type <- md.biox$ClusterName
md.biox$class <- md.biox$ClusterName
md.biox$class <- as.character(md.biox$class)

biox_class_annotations <- c('Neur_Int_MGE' = 'Inhibitory neurons', 
    'Neur_Int_LGE_CGE' = 'Inhibitory neurons', 
    'Neur_GABA_1' = 'Inhibitory neurons', 
    'Neur_Tel_Glut_SATB2_1' = 'Excitatory neurons',
    'Neur_Tel_Glut_SATB2_2' = 'Excitatory neurons',
    'Neur_Tel_Glut_NBL_1' = 'neuronal IPC',
    'Neur_Tel_Glut_NBL_2' = 'neuronal IPC',
    'Neur_Tel_Glut_PCNA' = 'neuronal IPC',
    'Rgl_Tel_2' = 'Radial glia',
    'Rgl_Tel_ventral' = 'Radial glia',
    'Rgl_Tel_dorsal' = 'Radial glia',
    'Rgl_Tel_1' = 'Radial glia'
) 

md.biox$class <- biox_class_annotations[as.character(md.biox$class)]
md.biox$Sample <- md.biox$sample 

# Remove the columns 
md.biox <- md.biox %>%
  select(-sample, -Clusters, -CellCycle, -Chemistry, -Class, -ClusterName, -Donor, -FRIP, -FractionTSS, -NGenes, -NPeaks, -SEX, -TotalUMI, -passed_filters, -regions)

md.advance$development_stage <- as.character(md.advance$development_stage)

advance_age_annotations <- c("18th week post-fertilization human stage" = "pcw16", "19th week post-fertilization human stage" = "pcw17",
				"23rd week post-fertilization human stage" = "pcw21", "24th week post-fertilization human stage" = "pcw22")

md.advance$Age <- advance_age_annotations[md.advance$development_stage]
md.advance$class <- md.advance$cell_type
md.advance$class <- as.character(md.advance$class)

advance_class_annotations <- c('vascular associated smooth muscle cell' = 'Pericytes',
    'radial glial cell' = 'Radial glia',
    'pericyte'= 'Pericytes',
    'oligodendrocyte precursor cell' ='OPC/Oligo',
    'neural progenitor cell' = 'neuronal IPC',
    'microglial cell' = 'Microglia',
    'medial ganglionic eminence derived interneuron' = 'Inhibitory neurons',
    'inhibitory interneuron' = 'Inhibitory neurons',
    'glutamatergic neuron' = 'Excitatory neurons',
    'endothelial cell' = 'Endothelial',
    'caudal ganglionic eminence derived interneuron' = 'Inhibitory neurons',
    'astrocyte' = 'Astrocytes/Oligo') 

md.advance$class <- advance_class_annotations[as.character(md.advance$class)]
md.advance$method <- 'multiome_atac'

advance_sample_annotations <- c("18th week post-fertilization human stage" = "sample1", "19th week post-fertilization human stage" = "sample2",
                                "23rd week post-fertilization human stage" = "sample3", "24th week post-fertilization human stage" = "sample4")

md.advance$Sample <- advance_sample_annotations[md.advance$development_stage] 
md.advance <- md.advance %>%
  select(-author_cell_type, -age_group, -donor_id, -nCount_RNA, -nFeature_RNA, -nCount_ATAC, -nFeature_ATAC, -TSS_percentile, -nucleosome_signal, -percent_mt, -assay_ontology_term_id, -cell_type_ontology_term_id, -development_stage_ontology_term_id, -disease_ontology_term_id, -self_reported_ethnicity_ontology_term_id, -organism_ontology_term_id, -sex_ontology_term_id, -tissue_ontology_term_id, -suspension_type, -is_primary_data, -batch, -assay, -disease, -organism, -sex, -tissue, -self_reported_ethnicity, -development_stage)


#create fragment objects
path_cell <- c(
    '/cfs/klemming/projects/supr/snic2022-23-547/yiquan/Cell/data/mfragments/m_hft_w16_p7_r1_fragments.tsv.gz',
    '/cfs/klemming/projects/supr/snic2022-23-547/yiquan/Cell/data/mfragments/m_hft_w16_p7_r2_fragments.tsv.gz',
    '/cfs/klemming/projects/supr/snic2022-23-547/yiquan/Cell/data/mfragments/m_hft_w16_p7_r3_fragments.tsv.gz',
    '/cfs/klemming/projects/supr/snic2022-23-547/yiquan/Cell/data/mfragments/m_hft_w20_p3_r1_fragments.tsv.gz',
    '/cfs/klemming/projects/supr/snic2022-23-547/yiquan/Cell/data/mfragments/m_hft_w20_p3_r2_fragments.tsv.gz',
    '/cfs/klemming/projects/supr/snic2022-23-547/yiquan/Cell/data/mfragments/m_hft_w21_p5_r1_fragments.tsv.gz',
    '/cfs/klemming/projects/supr/snic2022-23-547/yiquan/Cell/data/mfragments/m_hft_w21_p5_r2_fragments.tsv.gz',
    '/cfs/klemming/projects/supr/snic2022-23-547/yiquan/Cell/data/mfragments/m_hft_w21_p5_r3_fragments.tsv.gz',
    '/cfs/klemming/projects/supr/snic2022-23-547/yiquan/Cell/data/mfragments/m_hft_w21_p5_r4_fragments.tsv.gz',
    '/cfs/klemming/projects/supr/snic2022-23-547/yiquan/Cell/data/mfragments/m_hft_w24_p6_r1_fragments.tsv.gz',
    '/cfs/klemming/projects/supr/snic2022-23-547/yiquan/Cell/data/mfragments/m_hft_w24_p6_r2_fragments.tsv.gz',
    '/cfs/klemming/projects/supr/snic2022-23-547/yiquan/Cell/data/mfragments/m_hft_w24_p6_r3_fragments.tsv.gz',
    '/cfs/klemming/projects/supr/snic2022-23-547/yiquan/Cell/data/mfragments/m_hft_w24_p6_r4_fragments.tsv.gz'
)

# Create a list of Fragment objects
frag_cell <- lapply(path_cell, function(path) CreateFragmentObject(path = path, cells = rownames(md.cell), validate.fragments = FALSE))

#create fragment objects
path_nature <- c('/cfs/klemming/projects/supr/snic2022-23-547/yiquan/Nature/fragment/Cortex_GW17_fragments.tsv.gz', 
                 '/cfs/klemming/projects/supr/snic2022-23-547/yiquan/Nature/fragment/Cortex_GW18_fragments.tsv.gz', 
                 '/cfs/klemming/projects/supr/snic2022-23-547/yiquan/Nature/fragment/Cortex_GW21_fragments.tsv.gz', 
                 '/cfs/klemming/projects/supr/snic2022-23-547/yiquan/Nature/fragment/Insula_GW20_fragments.tsv.gz', 
                 '/cfs/klemming/projects/supr/snic2022-23-547/yiquan/Nature/fragment/M1_GW20_fragments.tsv.gz', 
                 '/cfs/klemming/projects/supr/snic2022-23-547/yiquan/Nature/fragment/MGE_GW20_fragments.tsv.gz', 
                 '/cfs/klemming/projects/supr/snic2022-23-547/yiquan/Nature/fragment/Parietal_GW20_fragments.tsv.gz', 
                 '/cfs/klemming/projects/supr/snic2022-23-547/yiquan/Nature/fragment/PFC_GW20_fragments.tsv.gz',
                 '/cfs/klemming/projects/supr/snic2022-23-547/yiquan/Nature/fragment/Somato_GW20_fragments.tsv.gz',
                 '/cfs/klemming/projects/supr/snic2022-23-547/yiquan/Nature/fragment/V1_GW20_fragments.tsv.gz',
                 '/cfs/klemming/projects/supr/snic2022-23-547/yiquan/Nature/fragment/Temporal_GW20_fragments.tsv.gz'
)

# Create a list of Fragment objects
frag_nature <- lapply(path_nature, function(path) CreateFragmentObject(path = path, cells = rownames(md.nature), validate.fragments = FALSE))

path_DC <- c('/cfs/klemming/projects/supr/snic2022-23-547/yiquan/DC/fragments/m_HES1_fragments.tsv.gz',
             '/cfs/klemming/projects/supr/snic2022-23-547/yiquan/DC/fragments/m_P18861_1001_fragments.tsv.gz',
             '/cfs/klemming/projects/supr/snic2022-23-547/yiquan/DC/fragments/m_P18861_1002_fragments.tsv.gz'
)

# Create a list of Fragment objects
frag_DC <- lapply(path_DC, function(path) CreateFragmentObject(path = path, cells = rownames(md.DC),validate.fragments = FALSE))

path_biox <-c(
    '/cfs/klemming/projects/supr/snic2022-23-547/yiquan/Biox/fragment/m_10X291_2.tsv.gz',
    '/cfs/klemming/projects/supr/snic2022-23-547/yiquan/Biox/fragment/m_10X291_3.tsv.gz',
    '/cfs/klemming/projects/supr/snic2022-23-547/yiquan/Biox/fragment/m_10X313_1.tsv.gz',
    '/cfs/klemming/projects/supr/snic2022-23-547/yiquan/Biox/fragment/m_10X313_2.tsv.gz',
    '/cfs/klemming/projects/supr/snic2022-23-547/yiquan/Biox/fragment/m_10X346_2.tsv.gz',
    '/cfs/klemming/projects/supr/snic2022-23-547/yiquan/Biox/fragment/m_10X347_2.tsv.gz',
    '/cfs/klemming/projects/supr/snic2022-23-547/yiquan/Biox/fragment/m_10X366_4.tsv.gz',
    '/cfs/klemming/projects/supr/snic2022-23-547/yiquan/Biox/fragment/m_10X369_1.tsv.gz',
    '/cfs/klemming/projects/supr/snic2022-23-547/yiquan/Biox/fragment/m_10X370_1.tsv.gz',
    '/cfs/klemming/projects/supr/snic2022-23-547/yiquan/Biox/fragment/m_10X402_1.tsv.gz',
    '/cfs/klemming/projects/supr/snic2022-23-547/yiquan/Biox/fragment/m_10X406_7.tsv.gz',
    '/cfs/klemming/projects/supr/snic2022-23-547/yiquan/Biox/fragment/m_10X420_1.tsv.gz'
)

# Create a list of Fragment objects
frag_biox <- lapply(path_biox, function(path) CreateFragmentObject(path = path, cells = rownames(md.biox), validate.fragments = FALSE))

path_advance <- c('/cfs/klemming/projects/supr/snic2022-23-547/yiquan/Advance/m_4.tsv.gz',
             	 '/cfs/klemming/projects/supr/snic2022-23-547/yiquan/Advance/m_8.tsv.gz',
		 '/cfs/klemming/projects/supr/snic2022-23-547/yiquan/Advance/m_11.tsv.gz',
                 '/cfs/klemming/projects/supr/snic2022-23-547/yiquan/Advance/m_16.tsv.gz'
)

# Create a list of Fragment objects
frag_advance <- lapply(path_advance, function(path) CreateFragmentObject(path = path, cells = rownames(md.advance),validate.fragments = FALSE))

#Quantify peaks in each dataset
cell.counts <- FeatureMatrix(
  fragments = frag_cell,
  features = combined.peaks,
  cells = rownames(md.cell)
)
#Create the objects
cell_assay <- CreateChromatinAssay(cell.counts, fragments = frag_cell)
cell <- CreateSeuratObject(cell_assay, assay = "ATAC", meta.data=md.cell)
                
nature.counts <- FeatureMatrix(
  fragments = frag_nature,
  features = combined.peaks,
  cells = rownames(md.nature)
)

nature_assay <- CreateChromatinAssay(nature.counts, fragments = frag_nature)
nature <- CreateSeuratObject(nature_assay, assay = "ATAC", meta.data=md.nature)

DC.counts <- FeatureMatrix(
  fragments = frag_DC,
  features = combined.peaks,
  cells = rownames(md.DC)
)
DC_assay <- CreateChromatinAssay(DC.counts, fragments = frag_DC)
DC <- CreateSeuratObject(DC_assay, assay = "ATAC", meta.data=md.DC)

biox.counts <- FeatureMatrix(
  fragments = frag_biox,
  features = combined.peaks,
  cells = rownames(md.biox)
)

biox_assay <- CreateChromatinAssay(biox.counts, fragments = frag_biox)
biox <- CreateSeuratObject(biox_assay, assay = "ATAC", meta.data=md.biox)
                   
advance.counts <- FeatureMatrix(
  fragments = frag_advance,
  features = combined.peaks,
  cells = rownames(md.advance)
)

advance_assay <- CreateChromatinAssay(advance.counts, fragments = frag_advance)
advance <- CreateSeuratObject(advance_assay, assay = "ATAC", meta.data=md.advance)       

# add information to identify dataset of origin
cell$dataset <- 'cell'
nature$dataset <- 'nature'   
DC$dataset <- 'DC'
biox$dataset <- 'biox'
advance$dataset <- 'advance'

combined <- merge(
  x = cell,
  y = list(nature, DC, biox, advance)
)

nature <- RunTFIDF(nature)
nature <- FindTopFeatures(nature, min.cutoff = 20)
nature <- RunSVD(nature)
nature <- RunUMAP(nature, dims = 2:30, reduction = 'lsi')
                   
cell <- RunTFIDF(cell)
cell <- FindTopFeatures(cell, min.cutoff = 20)
cell <- RunSVD(cell)
cell <- RunUMAP(cell, dims = 2:30, reduction = 'lsi')
                   
biox <- RunTFIDF(biox)
biox <- FindTopFeatures(biox, min.cutoff = 20)
biox <- RunSVD(biox)
biox <- RunUMAP(biox, dims = 2:30, reduction = 'lsi')

DC <- RunTFIDF(DC)
DC <- FindTopFeatures(DC, min.cutoff = 20)
DC <- RunSVD(DC)
DC <- RunUMAP(DC, dims = 2:30, reduction = 'lsi')

advance <- RunTFIDF(advance)
advance <- FindTopFeatures(advance, min.cutoff = 20)
advance <- RunSVD(advance)
advance <- RunUMAP(advance, dims = 2:30, reduction = 'lsi')

combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 20)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, dims = 2:30, reduction = 'lsi')

p1 <- DimPlot(object = combined, label = TRUE, raster = FALSE, group.by = "dataset", pt.size = 0.01)

ggsave('merge5.png', p1, width = 17, height = 15, dpi= 300)

integration.anchors <- FindIntegrationAnchors(
  object.list = list(cell, DC, nature, biox, advanc),
  anchor.features = rownames(cell),
  reduction = "rlsi",
  dims = 2:30,
  k.filter = NA
)

# integrate LSI embeddings
integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = combined[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30
)

# create a new UMAP using the integrated embeddings
integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30)

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to UCSC style since the data was mapped to hg38
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"
Annotation(integrated) <- annotations

integrated$Age <- as.numeric(gsub("pcw", "", integrated2$Age))
custom_breaks <- seq(7, 24, by = 6)

p4 <- DimPlot(object = integrated, label = TRUE, raster = FALSE, group.by ="class", pt.size= 0.01) + NoLegend()
p5 <- FeaturePlot (object = integrated, features = "Age", raster = FALSE, pt.size= 0.01) + scale_colour_viridis_c(breaks=custom_breaks)

p6 <-  p4 | p5

p7 <- DimPlot(object = integrated, label = TRUE, raster = FALSE, group.by = "cell_type", pt.size=0.01) + NoLegend()
p8 <- DimPlot(object = integrated, label = TRUE, raster = FALSE, group.by = "dataset", pt.size = 0.01)

p9 <-  p7 | p8

ggsave('coplot_signac5_1.png', p6, width= 30, height = 14, dpi=300)

ggsave('coplot_signac5_2.png', p9, width = 30, height = 14, dpi= 300)

saveRDS(integrated, '/cfs/klemming/projects/supr/snic2022-23-547/yiquan/out/signac_integrate5.rds')                   
