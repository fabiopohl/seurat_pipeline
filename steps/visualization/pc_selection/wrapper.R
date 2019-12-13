#' PC selection plots - Jackstraw and elbow
#'
#' @param seurat Seurat object .rds
#' @return .eps plots
suppressMessages(library(Seurat))

# Get input 
seurat_file <- snakemake@input[["seurat"]]

# Load Seurat Object
seurat <- readRDS(seurat_file)

# Get output
output_jack <- snakemake@output[["jackstraw"]]
output_elbow <- snakemake@output[["elbow"]]

# Plot Jackstraw 
setEPS()
postscript(output_jack, width = 8.0, height = 8.0)
JackStrawPlot(seurat, dim = 1:15)
temp <- dev.off()

# Plot elbow
setEPS()
postscript(output_elbow, width = 8.0, height = 8.0)
ElbowPlot(seurat)
temp <- dev.off()

# Testing
#seurat_file <- "/data/proj/GCB_FBP/repos/seurat_pipeline/steps/visualization/pc_selection/test/seurat_obj_red.rds"
