#' Collapse genes ID to genes symbol - select highest expressed transcript
#'
#' @param mat Full path to sparse matrix R Object
#' @param id_name Full path to CSV file with Gene ID mapping to Gene symbols
#' @return Collapsed matrix
library(dplyr)
library(Matrix)

collapse_genes <- function(mat, id_name) {
    
    t2g <- read.table(file = id_name, sep = ",",
                      header = TRUE, stringsAsFactors = FALSE)
    genes <- rownames(mat)
    genes_keep <- intersect(genes, t2g$Gene.stable.ID.version)
    t2g_keep <- filter(t2g, Gene.stable.ID.version %in% genes_keep, !duplicated(Gene.stable.ID.version))
    mat_keep <- mat[genes_keep, ]
    matmn <- data.frame(idx = 1:nrow(mat_keep), mn = rowMeans(mat_keep))
    matmn$geneid_version <- rownames(matmn)
    matmn_t2g <- inner_join(t2g_keep, matmn, by = c("Gene.stable.ID.version" = "geneid_version"))
    to_collapse <- matmn_t2g %>%
        group_by(Gene.name) %>%
        arrange(mn) %>%
        mutate(to_filter = !duplicated(Gene.name)) %>%
        ungroup() %>%
        arrange(idx)
    
    mat_final <- mat_keep[to_collapse$to_filter, ]
    rownames(mat_final) <- to_collapse[to_collapse$to_filter, "Gene.name"][[1]]
    return(mat_final)

}

# Get parameters
mat_path <- snakemake@input[[mat]]

# Test
#genes <- read.table(file = "/Users/fabiopohl/proj/czi/data/czi_kb_lidx_191108/output_11723WAPool01-S__18_13782/counts_unfiltered/unspliced.genes.txt", stringsAsFactors = FALSE)
#barcodes <- read.table(file = "/Users/fabiopohl/proj/czi/data/czi_kb_lidx_191108/output_11723WAPool01-S__18_13782/counts_unfiltered/spliced_unspliced_barcodes.txt", stringsAsFactors = FALSE)
#mat <- readMM("/Users/fabiopohl/proj/czi/data/czi_kb_lidx_191108/output_11723WAPool01-S__18_13782/counts_unfiltered/spliced_unspliced.mtx")
#t2g <- read.table(file = "/Users/fabiopohl/references/geneid_genename.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
#
## Transpose matrix
#mat <- t(mat)
#rownames(mat) <- genes[[1]]
#colnames(mat) <- barcodes[[1]]
#
#test <- collapse_genes(mat, t2g)
