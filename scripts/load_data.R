#load_data.R
#Loads in the Visium HD data and saves as a .rds for downstream analysis.
#
#Dataset is pulled from: https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-mouse-intestine
#
#Notes:
# - only needs to be run once, to generate the initial .rds file

library(here)
library(Seurat)

here::i_am("repo/scripts/load_data.R")
localdir <- here("data")
intestine <- Load10X_Spatial(data.dir = localdir, bin.size = c(8,16))
saveRDS(intestine, file = here("rds_objects", "intestine.rds"))
