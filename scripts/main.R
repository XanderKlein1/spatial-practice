#main.R
#Runs the main analysis workflow
#
#Notes:
#Top-level execution.

#Load prerequisite libraries
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(here)

#Load in the data
here::i_am("repo/scripts/main.R")
source_utils <- list.files(path = here("repo/utils"), pattern = "\\.R$", full.names = TRUE)
lapply(source_utils, source)

localdir <- here("rds_objects", "intestine.rds")
intestine <- readRDS(localdir)

#We will do all analysis on the 8 um bin data
DefaultAssay(intestine) <- "Spatial.008um"

#Preprocess the data:
intestine <- run_precluster(intestine, 2000)

#Looking at the elbow plot, we can determine principal components we should use for downstream analysis.
#In this case, it looks like the variance tapers off around 18 so we will use 18 PCs.
elb_p <- ElbowPlot(intestine)
ggsave(filename = here("repo", "figures", "elbow_plot.pdf"),
       plot = elb_p,
       width = 7,
       height = 5,
       dpi = 300)

#Cluster the data and project with UMAP:
intestine <- cluster_umap(intestine, 18, 0.8)
saveRDS(intestine, file = here("rds_objects", "intestine_analysis.rds"))

interaction_heatmap <- run_interaction_analysis(intestine, 4, 250, 100)
ggsave(filename = here("repo", "figures", "interaction_heatmap.pdf"),
       plot = interaction_heatmap,
       width = 7,
       height = 5,
       dpi = 300)

cluster_heatmap <- create_cluster_heatmap(intestine, n_genes = 5)
