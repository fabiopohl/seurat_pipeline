# Test functions for Seurat pipeline

# Sum matrix
source("./unspliced_plus_spliced.R")

mat <- unspliced_plus_spliced("~/proj/steve/invitro_preimplantation/data/AB25/")

# Collapse matrix
source("./collapse_genes.R")

mat_collapsed <- collapse_genes(mat)

# Remove empty droplets, create Seurat obj and normalise (SCTransform)
source("./filter_normalise.R")

seurat <- seurat_create_normalise(mat_collapsed, 100, 0.001, "ab25")
