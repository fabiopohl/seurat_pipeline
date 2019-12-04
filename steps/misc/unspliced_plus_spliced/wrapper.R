#' Merge unspliced and spliced matrix from Kallisto KB output
#'
#' @param sample_path Full path to sample
#' @param save_mat Save matrix (Default = TRUE)
#' @return Named sparse matrix
library(Matrix)
unspliced_plus_spliced <- function(sample_path, save_mat = TRUE) {
    # Read data
    # Unspliced
    genes_unspliced <- read.table(file = paste0(sample_path, "counts_unfiltered/unspliced.genes.txt"),
                                  header = FALSE, stringsAsFactors = FALSE)
    barcodes_unspliced <- read.table(file = paste0(sample_path, "counts_unfiltered/unspliced.barcodes.txt"),
                                     header = FALSE, stringsAsFactors = FALSE)
    mat_unspliced <- readMM(paste0(sample_path, "counts_unfiltered/unspliced.mtx"))
    rownames(mat_unspliced) <- barcodes_unspliced[[1]]
    
    # Spliced
    genes_spliced <- read.table(file = paste0(sample_path, "counts_unfiltered/spliced.genes.txt"),
                                              header = FALSE, stringsAsFactors = FALSE)
    barcodes_spliced <- read.table(file = paste0(sample_path, "counts_unfiltered/spliced.barcodes.txt"),
                                   header = FALSE, stringsAsFactors = FALSE)
    mat_spliced <- readMM(paste0(sample_path, "counts_unfiltered/spliced.mtx"))
    rownames(mat_spliced) <- barcodes_spliced[[1]]
    
    # Intersect by barcodes
    barcodes_common <- intersect(barcodes_unspliced[[1]], barcodes_spliced[[1]])
    mat_unspliced_common <- mat_unspliced[barcodes_common, ]
    mat_spliced_common <- mat_spliced[barcodes_common, ]
    
    # Sum matrixs
    mat_sum <- mat_unspliced_common + mat_spliced_common
    rownames(mat_sum) <- NULL
    
    # Save matrix and new barcode table
    if (save_mat == TRUE) {

        writeMM(file = paste0(sample_path, "counts_unfiltered/spliced_unspliced.mtx"), mat_sum)
        write.table(file = paste0(sample_path, "counts_unfiltered/spliced_unspliced_barcodes.txt"), data.frame(barcodes_common))
        write.table(file = paste0(sample_path, "counts_unfiltered/spliced_unspliced_genes.txt"), genes_unspliced)

    }

    rownames(mat_sum) <- barcodes_common
    colnames(mat_sum) <- genes_unspliced[[1]]
    mat_sum <- t(mat_sum)
    return(mat_sum)
}

# Get parameters
sample_path <- snakemake@input[["sample"]]
