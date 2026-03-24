#Load prerequisite libraries
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

#Read in the output from the Main script.
intestine <- readRDS("./outputs/intestineUMAP.rds")

#QC plot
qc_p <- SpatialFeaturePlot(intestine, features = c("nCount_Spatial.008um"),
                           ncol = 2)
ggsave(filename = "qc_plot.pdf",
       plot = qc_p,
       width = 7,
       height = 5,
       dpi = 300)

#QC plot of the normalized and scaled expression data (what actually gets used downstream)

#Generate an elbow plot to visualize which dimensions from PCA are important.
elb_p <- ElbowPlot(intestine)
ggsave(filename = "elbow_plot.pdf",
       plot = elb_p,
       width = 7,
       height = 5,
       dpi = 300)

#Create UMAP plot to visualize the clustering.
umap_p <- DimPlot(intestine, reduction = "umap")
ggsave(filename = "umap_plot.pdf",
       plot = umap_p,
       width = 7,
       height = 5,
       dpi = 300)

#Overlay UMAP onto tissue image.
umap_spatial <- SpatialDimPlot(intestine, label = TRUE, label.size = 3)
ggsave(filename = "umap_spatial.pdf",
       plot = umap_spatial,
       width = 7,
       height = 5,
       dpi = 300)