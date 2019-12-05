#' Merge unspliced and spliced matrix from Kallisto KB output
#'
#' @param mat_spliced Full path to spliced matrix from KB
#' @param barcodes_spliced Full path to spliced barcodes table
#' @param genes_spliced Full path to spliced genes table
#' @param mat_unspliced Full path to unspliced matrix from KB
#' @param barcodes_unspliced Full path to unspliced barcodes table
#' @param genes_unspliced Full path to unspliced genes table
#' @return Unspliced + spliced matrix
library(Matrix)
unspliced_plus_spliced <- function(mat_spliced, barcodes_spliced, genes_spliced,
                                   mat_unspliced, barcodes_unspliced, genes_unspliced,
                                   mat_output, barcodes_output, genes_output) {

    # Unspliced
    genes_unspliced <- read.table(file = genes_unspliced, header = FALSE, stringsAsFactors = FALSE)
    barcodes_unspliced <- read.table(file = barcodes_unspliced, header = FALSE, stringsAsFactors = FALSE)
    mat_unspliced <- readMM(mat_unspliced)
    rownames(mat_unspliced) <- barcodes_unspliced[[1]]
    
    # Spliced
    genes_spliced <- read.table(file = genes_spliced, header = FALSE, stringsAsFactors = FALSE)
    barcodes_spliced <- read.table(file = barcodes_spliced, header = FALSE, stringsAsFactors = FALSE)
    mat_spliced <- readMM(mat_spliced)
    rownames(mat_spliced) <- barcodes_spliced[[1]]
    
    # Intersect by barcodes
    barcodes_common <- intersect(barcodes_unspliced[[1]], barcodes_spliced[[1]])
    mat_unspliced_common <- mat_unspliced[barcodes_common, ]
    mat_spliced_common <- mat_spliced[barcodes_common, ]
    
    # Sum matrixs
    mat_sum <- mat_unspliced_common + mat_spliced_common
    rownames(mat_sum) <- NULL
    
    # Save matrix and new barcode table
    writeMM(file = mat_output, mat_sum)
    write.table(file = barcodes_output, data.frame(barcodes_common))
    write.table(file = genes_output, genes_unspliced)
}

# Get parameters
mat_spliced <- snakemake@input[["mat_spliced"]]
barcodes_spliced <- snakemake@input[["barcodes_spliced"]]
genes_spliced <- snakemake@input[["genes_spliced"]]

mat_unspliced <- snakemake@input[["mat_unspliced"]]
barcodes_unspliced <- snakemake@input[["barcodes_unspliced"]]
genes_unspliced <- snakemake@input[["genes_unspliced"]]

mat_output <- snakemake@output[["mat_output"]]
barcodes_output <- snakemake@output[["barcodes_output"]]
genes_output <- snakemake@output[["genes_output"]]

# Run function
unspliced_plus_spliced(mat_spliced, barcodes_spliced, genes_spliced,
                       mat_unspliced, barcodes_unspliced, genes_unspliced,
                       mat_output, barcodes_output, genes_output)
