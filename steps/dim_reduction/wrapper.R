#' Reduce dimensionality using Seurat
#' 
#' @param input Seurat R object
suppressMessages(library(Seurat))

# Get input
seurat_file <- snakemake@input

# Load Seurat object
seurat <- readRDS(seurat_file)

# Run PCA
seurat <- RunPCA(seurat)

# Jackstraw to select number of PCs
seurat <- JackStraw(seurat, num.replicate = 100)
seurat <- ScoreJackStraw(seurat, dim = 1:20)
sig_pcs <- sum(seurat@reductions$pca@jackstraw$overall.p.values < 0.01)

# Run UMAP
seurat <- RunUMAP(seurat, dim = 1:sig_pcs)

# Testing
seurat_file <- "/data/proj/GCB_FBP/repos/seurat_pipeline/steps/dim_reduction/test/seurat_obj.rds"

