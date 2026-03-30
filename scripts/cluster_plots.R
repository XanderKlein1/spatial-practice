#cluster_plots.R
#Generates 4 individual plots and 1 multi-plot image as pdf files:
# - Raw tissue image
# - UMAP projection of clusters
# - Clusters overlayed onto tissue image (with cluster labels)
# - Clusters overlayed onto tissue image (without cluster labels)


#Load prerequisite libraries
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(here)

#Read in the output from the Main script.
here::i_am("repo/scripts/cluster_plots.R")
localdir <- here("rds_objects", "intestine_analysis.rds")
intestine <- readRDS(localdir)
savedir <- here("repo", "figures")

#Plot the original tissue image
raw_image_p <- SpatialDimPlot(intestine, images = "slice1.008um", alpha = c(0,0), image.alpha = 1) + NoLegend()
ggsave(filename = here(savedir, "raw_image.pdf"),
       plot = raw_image_p,
       width = 7,
       height = 5,
       dpi = 300)

#Plot the UMAP projection of clusters
umap_p <- DimPlot(intestine, reduction = "umap")
ggsave(filename = here(savedir, "umap_plot.pdf"),
       plot = umap_p,
       width = 7,
       height = 5,
       dpi = 300)

#Plot the clusters overlayed onto the tissue image (with cluster labels)
cluster_overlay_p <- SpatialDimPlot(intestine, label = TRUE, label.size = 3)
ggsave(filename = here(savedir, "cluster_overlay.pdf"),
       plot = cluster_overlay_p,
       width = 7,
       height = 5,
       dpi = 300)

#Plot the clusters overlayed onto the tissue image (without cluster labels)
cluster_overlay_nolabel_p <- SpatialDimPlot(intestine, label = FALSE) + NoLegend()
ggsave(filename = here(savedir, "cluster_overlay_nolabel.pdf"),
       plot = cluster_overlay_nolabel_p,
       width = 7,
       height = 5,
       dpi = 300)

plots <- (raw_image_p | umap_p) / (cluster_overlay_p | cluster_overlay_nolabel_p)
ggsave(filename = here(savedir, "plots.pdf"),
       plot = plots,
       width = 9,
       height = 7,
       dpi = 300)
