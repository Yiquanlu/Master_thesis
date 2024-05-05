library(Seurat)
library(dplyr)

second <- readRDS("human_velme.rds")

columns_to_keep <- c("sample", "region", "age", "age_range", "chemistry", "lineage", "cell_type")

metadata_second <- second@meta.data[, columns_to_keep]

counts_second <- second[["RNA"]]@counts

metadata_second$dataset <- "velme"

second <-  CreateSeuratObject(counts = counts_second, project = "velme", assay = "RNA", meta.data = metadata_second)

second

first <- readRDS("/cfs/klemming/projects/supr/snic2022-23-547/yiquan/first_tri/first_trimemster.rds")

columns_to_keep <- c("Class", "CellClass", "Age", "Chemistry", "SampleID", "Subregion")

# Subset the metadata to include only the specified columns
metadata_first <- first@meta.data[, columns_to_keep]

metadata_first$age_range <- "1st trimemster"

metadata_first <- metadata_first %>%
  rename(region = Subregion,
         sample = SampleID,
         age = Age,
         chemistry = Chemistry,
         lineage = Class,
         cell_type = CellClass)

metadata_first$dataset <- "braun"

metadata_first$lineage <- metadata_first$cell_type

counts_first <- first[["RNA"]]@counts

first <-  CreateSeuratObject(counts = counts_first, project = "first.trimester", assay = "RNA", meta.data = metadata_first)

first

overlapping_features <- intersect(rownames(first), rownames(second))

first_subset <- subset(first, features = overlapping_features)
second_subset <- subset(second, features = overlapping_features)

merged_seurat <- merge(x = first_subset, y = second_subset, add.cell.ids = c("braun", "velme"))

merged_seurat

head(merged_seurat)

table(merged_seurat$lineage)
table(merged_seurat$cell_type)

saveRDS(merged_seurat, "human_RNA_merge_all.rds")
