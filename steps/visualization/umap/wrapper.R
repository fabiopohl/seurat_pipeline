#' UMAP for clustering
#' 
#' @param seurat Seurat object after UMAP and clustering
#' @return UMAP plots
suppressMessages(library(Seurat))
suppressMessages(library(plotly))

# Get input
seurat_file <- snakemake@input[["seurat"]]

# Get output
umap_file <- snakemake@output[["umap"]]
html_file <- snakemake@output[["html"]]
wd <- getwd()

# Load seurat object
seurat <- readRDS(seurat_file)

# Plot 2D UMAP
setEPS()
postscript(umap_file)
DimPlot(seurat, reduction = "umap")
temp <- dev.off()

# Plot 3D UMAP
umap_emb <- Embeddings(seurat, reduction = "umap")[, 1:3]
p <- plot_ly(type = "scatter3d", x = umap_emb[, 1], y = umap_emb[, 2], z = umap_emb[, 3], mode = "markers", color = seurat$seurat_clusters, text = seurat$seurat_clusters, size = 1)

htmlwidgets::saveWidget(as_widget(p), paste0(wd, "/", html_file))

# Testing
#seurat_file <- "/data/proj/GCB_FBP/repos/seurat_pipeline/steps/visualization/umap/test/seurat_obj_red_cluster.rds"
#umap_file <- "/data/proj/GCB_FBP/repos/seurat_pipeline/steps/visualization/umap/test/2d_umap.eps"
#html_file <- "/data/proj/GCB_FBP/repos/seurat_pipeline/steps/visualization/umap/test/3d_umap.html"
