#' Generate automatic report using Rmarkdown
#'
#' @param input Seurat object
#' @return Interactive report in HTML
suppressMessages(library(Seurat))
suppressMessages(library(plotly))
suppressMessages(library(rmarkdown))
suppressMessages(library(htmlwidgets))
suppressMessages(library(RColorBrewer))

library(knitr)
library(htmltools)

test <- "world"
test2 <- "hello"

# Get input
seurat_file <- snakemake@input[[1]]

# Load Seurat object
seurat <- readRDS(file = seurat_file)

# Plot 3D UMAP
umap_emb <- Embeddings(seurat, reduction = "umap")[, 1:3]
umap_cluster <- plot_ly(type = "scatter3d", x = umap_emb[, 1], y = umap_emb[, 2], z = umap_emb[, 3], mode = "markers", color = seurat$seurat_clusters, text = seurat$seurat_clusters, size = 1)

# Total UMI
feat <- "nCount_RNA"
cols <- brewer.pal(8, "OrRd")
pal <- colorRampPalette(cols)
color_feat <- seurat[[feat]][[1]]
umap_ntotal <- plot_ly(type = "scatter3d", x = umap_emb[, 1], y = umap_emb[, 2], z = umap_emb[, 3], mode = "markers", color = color_feat, size = 1)

# Total features
feat <- "nFeature_RNA"
cols <- brewer.pal(8, "OrRd")
pal <- colorRampPalette(cols)
color_feat <- seurat[[feat]][[1]]
umap_nfeats <- plot_ly(type = "scatter3d", x = umap_emb[, 1], y = umap_emb[, 2], z = umap_emb[, 3], mode = "markers", color = color_feat, colors = pal(20), size = 1)
umap_nfeats <- plot_ly(type = "scatter3d", x = umap_emb[, 1], y = umap_emb[, 2], z = umap_emb[, 3], mode = "markers", color = color_feat, size = 1)

# OPC markers
# PTPRZ1
feat <- "PTPRZ1"
color_counts <- as.vector(seurat[["RNA"]][feat, ])
ptp <- plot_ly(type = "scatter3d", x = umap_emb[, 1], y = umap_emb[, 2], z = umap_emb[, 3], mode = "markers", color = color_counts, colors = pal(20), size = 1)

# PDGFRA
feat <- "PDGFRA"
color_counts <- as.vector(seurat[["RNA"]][feat, ])
pdg <- plot_ly(type = "scatter3d", x = umap_emb[, 1], y = umap_emb[, 2], z = umap_emb[, 3], mode = "markers", color = color_counts, colors = pal(20), size = 1)

render("example.Rmd", output_file = "example_1.html")
# Test
seurat_file <- "test/seurat_obj_red_cluster.rds"
