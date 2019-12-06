#' Reduce dimensionality using Seurat
#' 
#' @param input Seurat R object
#' @return Seurat object with PCA and UMAP slots
suppressMessages(library(Seurat))

# Get input
set.seed(666)
seurat_file <- snakemake@input[[1]]

# Load Seurat object
seurat <- readRDS(seurat_file)

# Run PCA
seurat <- RunPCA(seurat)

# Jackstraw to select number of PCs
seurat <- JackStraw(seurat, num.replicate = 100)
seurat <- ScoreJackStraw(seurat, dim = 1:20)
sig_pcs <- sum(seurat@reductions$pca@jackstraw$overall.p.values < 0.01)

# Run UMAP
seurat <- RunUMAP(seurat, n.components = 3, dim = 1:sig_pcs)

# Get output
output_file <- snakemake@output[[1]]

# Save seurat object
saveRDS(file = output_file, seurat)

# Testing
#seurat_file <- "/data/proj/GCB_FBP/repos/seurat_pipeline/steps/dim_reduction/test/seurat_obj.rds"
