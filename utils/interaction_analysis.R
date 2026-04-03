#interaction_analysis.R
#Takes a pre-clustered seurat object and runs a cluster interaction analysis on it.
#
#Based on the analysis used in:
#Hildebrandt, F., Iturritza, M.U., Zwicker, C. et al. Host-pathogen interactions in the 
#Plasmodium-infected mouse liver at spatial and single-cell resolution. 
#Nat Commun 15, 7105 (2024). https://doi.org/10.1038/s41467-024-51418-2
#
#Algorithm Overview:
#Identify the k nearest neighbors for each cell.
#Count the interaction frequency between clusters (# of times cells from one cluster appear in another cluster's neighborhood).
#Normalize neighbor counts against cluster size.
#Generate a null model by running many random permutations of cluster labellings.
#Perform binomial tests to identify significant differences between the observed cluster interactions and the null model.

#Notes:
#100 permutations is generally enough for a good null model.
#A good starting size for the neighborhood is 4.
#distance should be determined based on tissue size.
#
#This script defines reusable function(s) only -- no top level execution.
#
#Output: a Heatmap plot of the under/over represented cluster interactions.

#Load prerequisite libraries
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(FNN)
library(pheatmap)


run_interaction_analysis <- function(object, k, max_dist, n_perm) {
  coords <- GetTissueCoordinates(intestine)
  
  #Find 4 nearest neighbors for each point.
  nn <- get.knnx(data = coords[, c("x", "y")],
                 query = coords[, c("x", "y")],
                 k = k+1)
  nn_idx <- nn$nn.index[, -1]
  nn_dist <- nn$nn.dist[, -1]
  
  #Set values that exceed the maximum distance to NA
  nn_idx[nn_dist > max_dist] <- NA
  
  #Create a dataframe to connect each cell to its corresponding cluster.
  cluster_id <- as.character(Idents(intestine))
  spot_id <- names(Idents(intestine))
  spotxcluster_df <- data.frame(spot_id, cluster_id)
  
  #Merge the dataframe with the x,y coordinates of the spots. First we have to add a spot_id field to the coords df so it can merge.
  coords$spot_id <- rownames(coords)
  df_merged <- left_join(coords, spotxcluster_df, by = "spot_id")
  
  #Now we can calculate the interaction effect between clusters
  #First we generate an empty matrix to store interaction frequency between clusters (filled with zeroes for now)
  clusters <- as.character(sort(as.numeric(unique(df_merged$cluster_id))))
  interaction_mtx <- matrix(0,
                            nrow = length(clusters),
                            ncol = length(clusters),
                            dimnames = list(clusters, clusters))
  #Filter out all of the NA neighbors (as determined by our max_dist earlier):
  neighbors <- lapply(seq_len(nrow(nn_idx)), function(i) {
    nn_idx[i, !is.na(nn_idx[i, ])]
  })
  
  #Iterate over all our neighborhoods to fill the interaction matrix:
  for(i in seq_len(nrow(df_merged))) {
    ci <- df_merged$cluster_id[i]
    for (j in neighbors[[i]]) {
      cj <- df_merged$cluster_id[j]
      interaction_mtx[ci,cj] <- interaction_mtx[ci,cj] + 1
    }
  }
  
  #Convert the raw counts to fraction of total neighbors for each cluster.
  #This is important to normalize the interaction effects across varying cluster sizes.
  row_totals <- rowSums(interaction_mtx)
  interaction_mtx_norm <- sweep(interaction_mtx, 1, row_totals, "/")
  interaction_mtx_norm[is.na(interaction_mtx_norm)] <- 0
  
  
  #To further account for differing cluster sizes, we can do a random permutation of cluster positions to establish a baseline for 
  #random cluster interactions. Then we can identify interactions between clusters that occur more or less likely than expectd by chance.
  
  #We create a 22x22x100 object to store 100 different permutations in the form of 22x22 interaction matrices.
  perm_array <- array(0, 
                      dim = c(length(clusters), length(clusters), n_perm), 
                      dimnames = list(clusters, clusters, NULL))
  set.seed(123)
  
  #Iteratively generate the permutations and store them into the perm_array:
  for (p in seq_len(n_perm)) {
    perm_cluster_id <- sample(df_merged$cluster_id)
    
    perm_mtx <- matrix(0, 
                       nrow = length(clusters), 
                       ncol = length(clusters),
                       dimnames = list(clusters, clusters))
    
    for (i in seq_len(nrow(df_merged))) {
      ci <- perm_cluster_id[i]
      
      for (j in neighbors[[i]]) {
        cj <- perm_cluster_id[j]
        perm_mtx[ci, cj] <- perm_mtx[ci, cj] + 1
      }
    }
    #Before adding each permutation matrix we have to convert it to fraction of total neighbors again:
    perm_row_totals <- rowSums(perm_mtx)
    perm_mtx_norm <- sweep(perm_mtx, 1, perm_row_totals, "/")
    perm_mtx_norm[is.na(perm_mtx_norm)] <- 0
    perm_array[,,p] <- perm_mtx_norm
  }
  
  #Next we can compare the observed interaction effect against the null model (permutations).
  perm_mean <- apply(perm_array, MARGIN = c(1,2), mean)
  perm_sd <- apply(perm_array, MARGIN = c(1,2), sd)
  
  perm_upper <- perm_mean + 2*perm_sd
  perm_lower <- perm_mean - 2*perm_sd
  
  #Perform binomial tests to determine the likelihood of attaining our observations under the null.
  binom_mat <- matrix(0, 
                      nrow = length(clusters),
                      ncol = length(clusters),
                      dimnames = list(clusters, clusters))
  
  for (i in seq_len(length(clusters))) {
    for (j in seq_len(length(clusters))) {
      binom_value <- binom.test(interaction_mtx[i,j], 
                                row_totals[i], 
                                perm_mean[i,j],
                                alternative = c("two.sided"),
                                conf.level = 0.95)
      binom_mat[i,j] <- binom_value$p.value
    }
  }
  
  #Next we have to get directionality in our p-values to determine under/over represented interactions:
  sign_mat <- sign(interaction_mtx_norm - perm_mean)
  score_mat <- -log10(binom_mat)*sign_mat
  
  #Final heatmap of over/under-representation in the observed model compared to the null model.
  pheatmap(
    score_mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    main = "Cluster Interaction Analysis",
    color = colorRampPalette(c("deepskyblue", "white", "red"))(100),
    breaks = seq(-5,5,length.out = 101)
  )
}

