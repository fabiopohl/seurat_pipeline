#' Plot expression of features of interest in UMAP
#'
#' @param seurat Seurat object after UMAP
suppressMessages(library(Seurat))

# Get input
seurat_file <- snakemake@input[["seurat"]]
features_file <- snakemake@input[["features"]]

# Get output 
umap <- snakemake@ouput[["umap"]]

# Load Seurat Object
seurat <- readRDS(seurat_file)

# Load features to plot
feats <- read.table(features_file, stringsAsFeatures = FALSE)

# Plot features
setEPS()
postscript("~/feats.eps")
FeaturePlot(seurat, features = feats)
temp <- dev.off()

# Testing
seurat_file <- "/data/proj/GCB_FBP/repos/seurat_pipeline/steps/visualization/expression/test/seurat_obj_red_cluster.rds"
feats <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "PTPRZ1", "MBP")
