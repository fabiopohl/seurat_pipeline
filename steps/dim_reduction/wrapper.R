#' Reduce dimensionality using Seurat
#' 
#' @param input Seurat R object
suppressMessages(library(Seurat))

# Get input
seurat_file <- snakemake@input

# Load Seurat object
seurat <- load(seurat_file)

# Run PCA
seurat <- RunPCA(seurat)

# Jackstraw to select number of PCs
seurat <- JackStraw(seurat, num_replicate = 100)
seurat <- ScoreJackStraw(seurat)

# Testing


