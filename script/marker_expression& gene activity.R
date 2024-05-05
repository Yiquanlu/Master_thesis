library(Signac)
library(Seurat)
library(Rmagic)
library(viridis)
library(gridExtra)
library(ggplot2)

# marker test for scRNA-seq datasets
integrate <- readRDS("human_RNA_all_0.3.rds")

features <- c('HES1', 'VIM', 'PAX6', 'SLC1A3', 'DLX1', 'GFAP', 'AQP4', 'OLIG1', 'OLIG2', 'SOX10', 'PDGFRA', 'P2RY12', 'EOMES', 'SLC6A1', 'SATB2', 'GAD2', 'LHX6', 'EBF3')

p6 <- FeaturePlot(integrate, features= features, raster = FALSE, pt.size = 0.01, ncol = 3)

ggsave('human_marker_gene.png', p6,  width = 24, height = 42, dpi = 300)


#marker test for glia subcluster in scATAC-seq datasets
brain <- readRDS("/cfs/klemming/projects/snic/snic2022-23-547/yiquan/out/gliogenesis.rds")

DefaultAssay(brain) <- 'RNA'

features2 <- c('SOX2','FOS', 'JUN', 'SLC1A3', 'EGFR', 'DLX1', 'HOPX','SPARCL1', 'GFAP', 'ALDH1L1', 'PDGFRA', 'CSPG4', 'OLIG1', 'OLIG2', 'ETV4', 'SOX10', 'AQP4', 'ENPP6')

plot_list <- lapply(features2, function(feature) {
  FeaturePlot(
    object = brain,
    features = feature,
    pt.size = 0.1,
    max.cutoff = 'q95',
    ncol = 3,
    raster = FALSE
  ) + scale_colour_viridis_c()
})

# Combine plots into a grid
combined_plot <- gridExtra::grid.arrange(grobs = plot_list, ncol = 3)


# Save the combined plot
ggsave('glia_marker.png', combined_plot, width = 18, height = 30, dpi = 300)

brain <- magic(brain)

DefaultAssay(brain) <- 'MAGIC_RNA'

p1 <- VlnPlot(brain, ncol=3, features = c('SOX2','FOS', 'JUN', 'SLC1A3', 'EGFR', 'DLX1', 'HOPX', 'SPARCL1', 'GFAP', 'ALDH1L1', 'PDGFRA', 'CSPG4', 'OLIG1', 'OLIG2', 'ETV4', 'SOX10', 'AQP4', 'ENPP6'))

ggsave('vln_feature.png', p1, width = 22, height = 42, dpi = 300)

plot_list <- lapply(features2, function(feature) {
  FeaturePlot(
    object = brain,
    features = feature,
    pt.size = 0.1,
    max.cutoff = 'q95',
    ncol = 3 
  ) + scale_colour_viridis_c()
})

# Combine plots into a grid
combined_plot <- gridExtra::grid.arrange(grobs = plot_list, ncol = 3)

# Save the combined plot
ggsave('glia_marker_magic.png', combined_plot, width = 18, height = 30, dpi = 300)


#marker test for interneuron subcluster in scATAC-seq datasets
features2 <- c('SOX2', 'DLX1', 'DLX2', 'ASCL1', 'OLIG1', 'OLIG2', 'NES','NKX2-1', 'SLC32A1', 'GAD2', 'LHX6', 'SP8', 'SST', 'PVALB', 'VIP', 'LAMP5', 'SIX3', 'EBF1', 'NR2F1', 'NR2F2', 'NFIX')

brain <- readRDS("/cfs/klemming/projects/snic/snic2022-23-547/yiquan/out/interneurongenesis.rds")

DefaultAssay(brain) <- 'RNA'

# Use lapply to generate separate plots for each feature
plot_list <- lapply(features2, function(feature) {
  FeaturePlot(
    object = brain,
    features = feature,
    pt.size = 0.05,
    max.cutoff = 'q95',
    ncol = 3,
    raster = FALSE
  ) + scale_colour_viridis_c()
})

# Combine plots into a grid
combined_plot <- gridExtra::grid.arrange(grobs = plot_list, ncol = 3)


# Save the combined plot
ggsave('intern_marker.png', combined_plot, width = 22, height = 42, dpi = 300)

brain <- magic(brain)

DefaultAssay(brain) <- 'MAGIC_RNA'

# Use lapply to generate separate plots for each feature
plot_list <- lapply(features2, function(feature) {
  FeaturePlot(
    object = brain,
    features = feature,
    pt.size = 0.05,
    max.cutoff = 'q95',
    ncol = 3,
    raster = FALSE
  ) + scale_colour_viridis_c()
})

# Combine plots into a grid
combined_plot <- gridExtra::grid.arrange(grobs = plot_list, ncol = 3)

# Save the combined plot
ggsave('intern_marker_magic.png', combined_plot, width = 22, height = 42, dpi = 300)
