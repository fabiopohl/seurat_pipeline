#' Pseudo align FASTQ files to referece genome using Kallisto KB wrapper
#'
#' @param r1 Full path to reads 1 FASTQ
#' @param r2 Full path to reads 2 FASTQ
#' @param reference Full path to reference genome index
#' @param t2g Full path to transcripts to genes table
#' @param c1 Full path to cdna to capture
#' @param c2 Full path to introns to capture
#' @param whitelist Full path to 10x barcode whitelist
#' @param barcode_version Barcodes version
#' @param output_folder Full path to output folder

# Get inputs and parameters
r1_file <- snakemake@input[["r1"]]
r2_file <- snakemake@input[["r2"]]
rpairs <- paste(r1_file, r2_file, sep = " ")
allreads <- paste(rpairs, collapse = " ")
params <- snakemake@params
params <- params[names(params) != ""]

# Run kallisto kb
cmd <- paste(list("count", 
                  "-i", params[["reference"]],
                  "-w", params[["whitelist"]],
                  "-m 40G -t 1",
                  "-o", params[["output_folder"]],
                  "-g", params[["t2g"]],
                  "-x", params[["barcode_version"]],
                  "-c1", params[["c1"]],
                  "-c2", params[["c2"]],
                  "--lamanno", allreads), collapse = " ")

system2("/home/fabio/anaconda3/bin/kb", cmd)

# Testing
#r1 = "/data/proj/GCB_FBP/repos/seurat_pipeline/steps/align/kallisto_kb/test/sample_r1.fastq"
#r2 = "/data/proj/GCB_FBP/repos/seurat_pipeline/steps/align/kallisto_kb/test/sample_r2.fastq"
#allreads = paste(r1, r2, sep = " ")
#
#params <- list(
#reference = "/data/proj/GCB_FBP/references/index/kb/linnarson/index.idx",
#t2g = "/data/proj/GCB_FBP/references/index/kb/linnarson/transcripts_to_genes.txt",
#c1 = "/data/proj/GCB_FBP/references/index/kb/linnarson/cdna_transcripts_to_capture.txt",
#c2 = "/data/proj/GCB_FBP/references/index/kb/linnarson/intron_transcripts_to_capture.txt",
#whitelist = "/data/proj/GCB_FBP/references/whitelist/10xv3_whitelist.txt",
#barcode_version = "0,0,16:0,16,27:1,0,0",
#output_folder = "/data/proj/GCB_FBP/repos/seurat_pipeline/steps/align/kallisto_kb/test/")
#
#cmd <- paste(list("count", 
#                  "-i", params[["reference"]],
#                  "-w", params[["whitelist"]],
#                  "-m 40G -t 1",
#                  "-o", params[["output_folder"]],
#                  "-g", params[["t2g"]],
#                  "-x", params[["barcode_version"]],
#                  "-c1", params[["c1"]],
#                  "-c2", params[["c2"]],
#                  "--lamanno", allreads), collapse = " ")
#
#system2("/home/fabio/anaconda3/bin/kb", cmd)
