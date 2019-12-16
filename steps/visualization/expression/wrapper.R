#' Plot expression of features of interest in UMAP
#'
#' @param seurat Seurat object after UMAP
#' @return Feature plots 
suppressMessages(library(Seurat))
suppressMessages(library(plotly))
suppressMessages(library(RColorBrewer))

# Get input
seurat_file <- snakemake@input[["seurat"]]

# Get feature
feat <- snakemake@params[["feat"]]

# Get output 
umap_file <- snakemake@output[["umap"]]
html_file <- snakemake@output[["html"]]

# Load Seurat Object
seurat <- readRDS(seurat_file)
DefaultAssay(seurat) <- "RNA"

# Plot features
setEPS()
postscript(umap_file)
FeaturePlot(seurat, features = feat)
temp <- dev.off()

# Plot 3D UMAP
cols <- brewer.pal(9, "OrRd")
pal <- colorRampPalette(cols)

metadata <- colnames(seurat[[]])
umap_emb <- Embeddings(seurat, reduction = "umap")[, 1:3]

if (feat %in% metadata) {
    color_feat <- seurat[[feat]][[1]]
    p <- plot_ly(type = "scatter3d", x = umap_emb[, 1], y = umap_emb[, 2], z = umap_emb[, 3], mode = "markers", color = color_feat, colors = pal(20), size = 1)
    htmlwidgets::saveWidget(as_widget(p), html_file)
} else {
    color_counts <- as.vector(seurat[["RNA"]][feat, ])
    p <- plot_ly(type = "scatter3d", x = umap_emb[, 1], y = umap_emb[, 2], z = umap_emb[, 3], mode = "markers", color = color_counts, colors = pal(20), size = 1)
    htmlwidgets::saveWidget(as_widget(p), html_file)
}

# Testing
#seurat_file <- "/data/proj/GCB_FBP/repos/seurat_pipeline/steps/visualization/expression/test/seurat_obj_red_cluster.rds"
#umap_file <-  "/data/proj/GCB_FBP/repos/seurat_pipeline/steps/visualization/expression/test/umap_ncount.eps"
#html_file <-  "/data/proj/GCB_FBP/repos/seurat_pipeline/steps/visualization/expression/test/3dumap_ncount.html"
#feats <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "PTPRZ1", "MBP")
#features_file <- "/data/proj/GCB_FBP/repos/seurat_pipeline/steps/visualization/expression/test/features_example.txt"
#feat <- "nCount_RNA"
#feat <- "PTPRZ1"
