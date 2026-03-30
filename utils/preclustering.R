#preclustering.R -------------------
#Utility functions to prepare seurat object for clustering and downstream analysis.
#
#Includes functionality for:
# - normalizing data (via log transform)
# - finding variable features
# - scaling data (set mean and variance to (0,1) respectively)
# - principal component analysis (PCA)
#
#Outputs the updated seurat object
#Notes:
# - defines reusable function(s) only -- no top level execution.

library(Seurat)
library(dplyr)

here::i_am("repo/utils/preclustering.R")

precluster <- function(object, numFeatures) {
  object <- NormalizeData(object)
  object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = numFeatures)
  object <- ScaleData(object, features = VariableFeatures(object))
  object <- RunPCA(object, features = VariableFeatures(object))
  
  object
  
}

