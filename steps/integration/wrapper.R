#' Integrate datasets using anchor points
#'
#' @param seurat List of seurat objects
#' @param samples List of samples to integrate
#' @return Single seurat object with integrated data
suppressMessages(library(Seurat))
#library(Seurat, lib = "/home/fabio/R/x86_64-redhat-linux-gnu-library/3.6/")
options(future.globals.maxSize =  10485760000)

# Get input
seurat_files <- snakemake@input[["seurat"]]
samples <- snakemake@params[["samples"]]

# Load seurat files
seurat_list <- lapply(seurat_files, readRDS)

# Add sample prefix to barcodes
seurat_list <- lapply(1:length(seurat_list), function(idx) {
                          RenameCells(seurat_list[[idx]], new.names=paste0(samples[idx], "_", colnames(seurat_list[[idx]])))
})

# Redoing 
# Get counts data
new_seurat <- lapply(1:length(seurat_list), function(i) {
                         CreateSeuratObject(counts = (GetAssayData(object = seurat_list[[i]], slot = "counts")))
})

# Normalise
for (i in 1:length(new_seurat)) {
    new_seurat[[i]] <- SCTransform(new_seurat[[i]], verbose = FALSE)
}
seurat_list <- new_seurat

# Integrate datasets
seurat_features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3000)
seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = seurat_features)
seurat_anchors <- FindIntegrationAnchors(object.list = seurat_list, normalization.method = "SCT",
                                         anchor.features = seurat_features, verbose = FALSE)
seurat_integrated <- IntegrateData(anchorset = seurat_anchors, normalization.method = "SCT")

# Get output
output_file <- snakemake@output[[1]]

# Save object
saveRDS(file = output_file, seurat_integrated)

# Testing
#seurat_files <- c("test/seurat_obj_red_cluster_1.rds",  "test/seurat_obj_red_cluster_2.rds")
#samples <- c("sample1", "sample2")
#output_file <- "test/seurat_obj_integrated.rds"
