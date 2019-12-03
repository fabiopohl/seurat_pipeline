# Function to integrate datasets

seurat_integrate <- function(seurat_list) {
    seurat_features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3000)
    seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = seurat_features)
    seurat_anchors <- FindIntegrationAnchors(object.list = seurat_list, normalization.method = "SCT",
                                             anchor.features = seurat_features, verbose = FALSE)
    seurat_integrated <- IntegrateData(anchorset = seurat_anchors, normalization.method = "SCT")
}    
