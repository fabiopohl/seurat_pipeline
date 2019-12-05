#' Remove empty cells (EmptyDrops)
#'
#' @param mat Sparse matrix 
#' @param barcodes Barcodes table
#' @param genes Gene names table
#' @param limit Total number of reads to consider empty cell (Recommended: 100)
#' @param cutoff FDR cutoff for defining empty cell
#' @return Filtered matrix
suppressMessages(library(DropletUtils))
suppressMessages(library(Matrix))

filter_cells <- function(mat, limit, cutoff){

    # Remove empty droplets
    e.out <- emptyDrops(mat, lower = limit)
    summary(e.out$FDR <= cutoff)
    
    # Filter empty droplets
    mat_noempty <- mat[, which(e.out$FDR <= cutoff)]
    return(mat_noempty)
    
}

# Get input
mat_file <- snakemake@input[["mat"]]
barcodes_file <- snakemake@input[["barcodes"]]
genes_file <- snakemake@input[["genes"]]

# Get params
limit <- snakemake@params[["limit"]]
cutoff <- snakemake@params[["cutoff"]]

# Prepare matrix
matr <- readMM(mat_file)
barcodes <- read.table(file = barcodes_file, stringsAsFactors = FALSE)
genes <- read.table(file = genes_file, stringsAsFactors = FALSE)

matr <- t(matr)
colnames(matr) <- barcodes[[1]]
rownames(matr) <- genes[[1]]

# Run function
mat_filt <- filter_cells(matr, limit, cutoff)

# Get output
mat_output <- snakemake@output[["mat_output"]]
barcodes_output <- snakemake@output[["barcodes_output"]]
genes_output <- snakemake@output[["genes_output"]]

# Save output
write.table(file = barcodes_output, data.frame(colnames(mat_filt)))
write.table(file = genes_output, data.frame(rownames(mat_filt)))

colnames(mat_filt) <- NULL
rownames(mat_filt) <- NULL
mat_filt <- t(mat_filt)
writeMM(file = mat_output, mat_filt)

# Testing
#mat_file = "/data/proj/GCB_FBP/czi/data/czi_kb_lidx_191108/output_11723WAPool01-S__18_13782/counts_unfiltered/spliced_unspliced.mtx"
#barcodes_file = "/data/proj/GCB_FBP/czi/data/czi_kb_lidx_191108/output_11723WAPool01-S__18_13782/counts_unfiltered/spliced_unspliced_barcodes.txt"
#genes_file = "/data/proj/GCB_FBP/czi/data/czi_kb_lidx_191108/output_11723WAPool01-S__18_13782/counts_unfiltered/spliced_unspliced_genes.txt"
#limit <- 100
#cutoff <- 0.001
