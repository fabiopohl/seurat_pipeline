#' Clustering using Seurat
#'
#' @param input Seurat Object with reduced dimensions
#' @return Seurat Object with clusters information
suppressMessages(library(Seurat))

# Get input
seurat_file <- snakemake@input[[1]]

# Load Seurat object
seurat <- readRDS(seurat_file)
ndims <- length(seurat@commands$RunUMAP.SCT.pca$dims)

# Cluster cells
seurat <- FindNeighbors(seurat, dim = 1:ndims)
seurat <- FindClusters(seurat, resolution = 2)

# Get output
output <- snakemake@output[[1]]

# Save Seurat Object
saveRDS(file = output, seurat)

# Testing
#seurat_file <- "/data/proj/GCB_FBP/repos/seurat_pipeline/steps/clustering/test/seurat_obj_red.rds"
