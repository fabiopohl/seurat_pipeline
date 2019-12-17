#' Find marker genes using Seurat and compare to markers gene list
#'
#' @param seurat Seurat object with clustering information
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))

# Get input
seurat_file <- snakemake@input[["seurat"]]
markers_file <- snakemake@input[["markers"]]

# Load data
seurat <- readRDS(seurat_file)
markers_list <- read.csv(file = markers_file, stringsAsFactors = FALSE)

# Find markers
seurat_markers <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Compare to markers list
celltypes <- unique(markers_list$CellType)

cells <- lapply(celltypes, function(cell) {
                    markers <- filter(markers_list, CellType == cell)$Gene
                    marker_data <- filter(seurat_markers, gene %in% markers)
                    marker_data$cell_type <- rep(cell, nrow(marker_data))
                    marker_data
})
cell_dt <- do.call("rbind", cells)
cell_dt <- arrange(cell_dt, cluster, cell_type, pct.1, avg_logFC)

# Get output
seurat_markers_output <- snakemake@output[["seurat_output"]]
known_markers <- snakemake@output[["known_markers"]]

# Save output
write.table(file = seurat_markers_output, seurat_markers, row.names = FALSE)
write.table(file = known_markers, cell_dt, row.names = FALSE)

# Testing
#seurat_file <- "/data/proj/GCB_FBP/repos/seurat_pipeline/steps/findmarkers/test/seurat_obj_red_cluster.rds"
#markers_file <- "/data/proj/GCB_FBP/repos/seurat_pipeline/steps/findmarkers/test/MARKERS_gwen.csv"
