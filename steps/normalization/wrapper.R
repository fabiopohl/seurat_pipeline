#' Filter low quality cells and normalise data (SCTransform)
#'
#' @param mat Sparse matrix
#' @param genes Genes table
#' @param barcodes Barcodes table
#' @return Seurat R object
suppressMessages(library(Seurat))
suppressMessages(library(Matrix))

# Get input
mat_file <- snakemake@input[["mat"]]
genes_file <- snakemake@input[["genes"]]
barcodes_file <- snakemake@input[["barcodes"]]

# Get params
min_features <- snakemake@params[["min_features"]]
max_features <- snakemake@params[["max_features"]]
pct_mt <- snakemake@params[["pct_mt"]]
min_total <- snakemake@params[["min_total"]]
vars_regress <- snakemake@params[["vars_regress"]]

# Load data
mat <- readMM(mat_file)
genes <- read.table(file = genes_file, stringsAsFactors = FALSE)
barcodes <- read.table(file = barcodes_file, stringsAsFactors = FALSE)

# Adjust matrix
mat <- t(mat)
colnames(mat) <- barcodes[[1]]
rownames(mat) <- genes[[1]]

# Filter low total reads barcodes
#ntotal <- colSums(mat)
#mat_fil <- mat[, which(ntotal > min_total)]

# Create Seurat Object
seurat <- CreateSeuratObject(counts = mat, min.cells = 10, min.features = 500)
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")

# Filter low quality cells
seurat <- subset(seurat, subset = nFeature_RNA > min_features & nFeature_RNA < max_features & percent.mt < pct_mt& nCount_RNA > min_total)

# Normalise and scale data
if (vars_regress == "NULL") {
    vars_regress <- NULL
}
seurat <- SCTransform(seurat, vars.to.regress = vars_regress)

# Get output
output <- snakemake@output[[1]]

# Save Seurat obj
saveRDS(file = output, seurat)
