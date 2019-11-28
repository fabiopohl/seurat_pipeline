# Function to remove empty cells (EmptyDrops) and normalise data (SCTransform)
library(Seurat)
library(DropletUtils)

seurat_create_normalise <- function(mat, limit, cutoff, sample_name, regress = NULL) {
    
    # Remove empty droplets
    e.out <- emptyDrops(mat, lower = limit)
    #summary(e.out$FDR <= cutoff)
    
    # Filter empty droplets
    mat_noempty <- mat[, which(e.out$FDR <= cutoff)]
    colnames(mat_noempty) <- paste0(sample_name, colnames(mat_noempty))
    
    # Create Seurat Object
    seurat <- CreateSeuratObject(counts = mat_noempty, project = sample_name, min.cells = 3, min.features = 200)
    seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
    
    # Filter low quality cells
    seurat <- subset(seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5 & nCount_RNA > 300)
    
    # Normalise and scale data
    seurat <- SCTransform(seurat, vars.to.regress = regress)
    return(seurat)

}
