
#Load prerequisite libraries
library(future)
plan(sequential)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

#Increase the maximum size of global variables (our Seurat object) so it can work with SCTransform.
options(future.globals.maxSize = 3000 * 1024^2)

intestine <- readRDS("./data/intestine.rds")
DefaultAssay(intestine) <- "Spatial.016um"
intestine <- SCTransform(intestine, assay = "Spatial.016um", method = 'glmGamPoi', verbose = FALSE, ncells=10000, conserve.memory=TRUE)

saveRDS(intestine, file = "./outputs/intestineSCT16um.rds")