#' Reduce dimensionality using Seurat
#' 
suppressMessages(library(Seurat))

# Get input
seurat_file <- snakemake@input

# Load Seurat object
seurat <- load(seurat_file)


