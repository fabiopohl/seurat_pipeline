#' Visualizations for sample quality control
#'
#' @param mat Raw sparse matrix
#' @param genes Genes table for matrix
#' @return Plots for quality control
suppressMessages(library(Matrix))
suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(DropletUtils))
suppressMessages(library(cowplot))

# Get input 
mat_file <- snakemake@input[["mat"]]
genes_file <- snakemake@input[["genes"]]

# Load data
mat <- readMM(mat_file)
genes <- read.table(file = genes_file, stringsAsFactors = FALSE)

mat <- t(mat)
colnames(mat) <- paste0("cell_", 1:ncol(mat))
rownames(mat) <- genes[[1]]

# Create seurat object
seurat <- CreateSeuratObject(counts = mat, min.cells = 3, min.features = 300)
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")

# Get output 
vln_file <- snakemake@output[["vln"]]
scatter_file <- snakemake@output[["scatter"]]
knee_file <- snakemake@output[["knee"]]

# Knee plot
bcrank <- barcodeRanks(mat)
uniq <- !duplicated(bcrank$rank)

setEPS() 
postscript(knee_file, width = 4.0, height = 4.0)
ggplot(as.data.frame(bcrank[uniq, ]), aes(rank, total)) +
    geom_point() +
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    labs(x = "Cell rank", y = "Total UMI count")
dev.off()

# Plot features
setEPS()
postscript(vln_file, width = 8.0, height = 4.0)
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
dev.off()

setEPS()
postscript(scatter_file, width = 8.0, height = 4.0)
p1 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
p2 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(p1, p2))
dev.off()

# Testing
mat_file <- "/data/proj/GCB_FBP/repos/seurat_pipeline/steps/visualization/qc/test/unspliced_spliced_collapsed.mtx"
genes_file <- "/data/proj/GCB_FBP/repos/seurat_pipeline/steps/visualization/qc/test/unspliced_spliced_collapsed.genes.txt"
