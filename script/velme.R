library(Seurat)
library(biomaRt)
library(dplyr)

second <- readRDS("velme.rds")

second 

meta1 <- read.delim(file = "meta_update.tsv", header = TRUE, row.names = 1, sep = "\t")

second <- AddMetaData(object = second, metadata = meta1)

filter <- second[,  second$age_range %in% c("2nd trimester", "3rd trimester")]

metadata <- filter@meta.data

counts <- filter[["RNA"]]@counts

# Define Ensembl and gene symbol attributes
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
attributes <- c("ensembl_gene_id", "external_gene_name")

# Get a mapping of Ensembl IDs to gene symbols
ensembl_to_symbol <- getBM(attributes = attributes, mart = ensembl)

# Match Ensembl IDs to gene symbols, handling missing values
gene_symbols <- character(length = nrow(counts))
match_indices <- match(rownames(counts), ensembl_to_symbol$ensembl_gene_id)
valid_matches <- !is.na(match_indices)
gene_symbols[valid_matches] <- ensembl_to_symbol$external_gene_name[match_indices[valid_matches]]

# Filter out rows with missing or empty gene symbols
non_empty_rows <- gene_symbols != "" & !is.na(gene_symbols)
counts_filtered <- counts[non_empty_rows, ]
gene_symbols_filtered <- gene_symbols[non_empty_rows]

# Assign gene symbols as row names
rownames(counts_filtered) <- gene_symbols_filtered

# Create a new Seurat object
filter <- CreateSeuratObject(counts = counts_filtered, project = "velme", assay = "RNA", meta.data = metadata)

filter 

# Save the updated Seurat object
saveRDS(filter, "human_velme.rds")
