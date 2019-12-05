#' Collapse genes ID to genes symbol - select highest expressed transcript
#'
#' @param mat Full path to sparse matrix R Object
#' @param barcodes Full path to barcodes table
#' @param genes Full path to genes table
#' @param id_to_genes Full path to CSV file with Gene ID mapping to Gene symbols
#' @return Collapsed matrix
library(dplyr)
library(Matrix)

# Get input
mat_file <- snakemake@input[["mat"]]
barcodes_file <- snakemake@input[["barcodes"]]
genes_file <- snakemake@input[["genes"]]

# Params
id_to_genes <- snakemake@params[["id_to_genes"]]

# Read files
t2g <- read.table(file = id_to_genes, sep = ",",
                  header = TRUE, stringsAsFactors = FALSE)
genes <- read.table(file = genes_file, stringsAsFactors = FALSE)
barcodes <- read.table(file = barcodes_file, stringsAsFactors = FALSE)
mat <- readMM(mat_file)

# Transpose matrix and add names
mat <- t(mat)
rownames(mat) <- genes[[1]]
colnames(mat) <- barcodes[[1]]

# Intersect genes
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

# Get output files
mat_output <- snakemake@output[["mat_output"]]
barcodes_output <- snakemake@output[["barcodes_output"]]
genes_output <- snakemake@output[["genes_output"]]

# Write output
write.table(file = genes_output, data.frame(rownames(mat_final)))
write.table(file = barcodes_output, data.frame(colnames(mat_final)))

colnames(mat_final) <- NULL
rownames(mat_final) <- NULL
writeMM(file = mat_output, mat_final)

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
