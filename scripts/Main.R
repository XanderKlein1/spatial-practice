#Load prerequisite libraries
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

#Load in the data
intestine <- readRDS("./data/intestine.rds")
DefaultAssay(intestine) <- "Spatial.008um"

#Normalize the data with a log transform.
intestine <- NormalizeData(intestine)

#Now we use the normalized data to find the genes with the highest variability.
intestine <- FindVariableFeatures(intestine, selection.method = "vst", nfeatures = 2000)

#Next we must scale the data by setting the mean and variance of the log transformed data to 0 and 1, respectively.
#This ensures that all genes are weighted equally in downstream analysis, and allows us to accurately identify lowly
#expressed DEGs.

#This is a computationally expensive algorithm, so we will have to run on the HPC and then bring the RDS back here.
intestine <- ScaleData(intestine, features = VariableFeatures(intestine))
intestine <- RunPCA(intestine, features = VariableFeatures(intestine))
saveRDS(intestine, file = "./outputs/intestinePCA.rds")

#We will use the first 18 PCs for downstream clustering (determined via elbow plot)
intestine <- FindNeighbors(intestine, dims = 1:18)
intestine <- FindClusters(intestine, resolution = 0.8)
saveRDS(intestine, file = "./outputs/intestineClustered.rds")

#Now we will run the UMAP algorithm to project our clusters into two dimensions.
intestine <- RunUMAP(intestine, dims = 1:18)
saveRDS(intestine, file = "./outputs/intestineUMAP.rds")