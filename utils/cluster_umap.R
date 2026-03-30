#cluster_umap.R -------------------
#Provides clustering functionality for spatial transcriptomics data.
#
#Includes functionality for:
# - finding neighborhoods
# - defining clusters
# - projecting clusters into 2 dimensions (UMAP)
#
#Notes:
# - defines reusable function(s) only -- no top level execution.

library(Seurat)

cluster_umap <- function(object, numPCA, clusterResolution) {
  object <- FindNeighbors(object, dims = 1:numPCA)
  object <- FindClusters(object, resolution = clusterResolution)
  object <- RunUMAP(object, dims = 1:numPCA)
  object
}