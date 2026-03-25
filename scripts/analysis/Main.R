#Load prerequisite libraries
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(here)

#Load in the data
localdir <- here("data")
intestine <- Load10X_Spatial(data.dir = localdir, bin.size = c(8,16))

#We will do all analysis on the 8 um bin data
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

#We will use the first 18 PCs for downstream clustering (determined via elbow plot)
intestine <- FindNeighbors(intestine, dims = 1:18)
intestine <- FindClusters(intestine, resolution = 0.8)

#Now we will run the UMAP algorithm to project our clusters into two dimensions.
intestine <- RunUMAP(intestine, dims = 1:18)
saveRDS(intestine, file = here("data", "intestine_analysis.rds"))
