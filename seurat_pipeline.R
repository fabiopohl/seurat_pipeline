# Seurat pipeline for fetal samples

# Get sample names
samples <- read.table(file = "../fetal_samples.txt", stringsAsFactors = FALSE)[[1]]
samples_dir <- "/data/proj/GCB_FBP/czi/data/czi_kb_lidx_191108/output_"

# Sum matrixs
source("./unspliced_plus_spliced.R")
mat_list <- lapply(paste0(samples_dir, samples, "/"), unspliced_plus_spliced)

# Collapse genes
source("./collapse_genes.R")
collapsed_mats <- lapply(mat_list, collapse_genes)

# Remove empty droplets, create Seurat object and normalise (SCTranform)
source("./filter_normalise.R")
seurat_list <- lapply(1:length(collapsed_mats), function(i) {
                          seurat_create_normalise(collapsed_mats[[i]], 100, 0.001, samples[[i]])
})

# Integrate datasets
source("./integrate_datasets.R")
seurat_integrated <- seurat_integrate(seurat_list)
