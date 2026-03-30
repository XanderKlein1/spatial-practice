#qc_plots.R
#Generates basic QC plots for the pre-processed data.

qc_p <- SpatialFeaturePlot(intestine, features = c("nCount_Spatial.008um"),
                           ncol = 2)
ggsave(filename = "qc_plot.pdf",
       plot = qc_p,
       width = 7,
       height = 5,
       dpi = 300)


#Generate an elbow plot to visualize which dimensions from PCA are important.
elb_p <- ElbowPlot(intestine)
ggsave(filename = "elbow_plot.pdf",
       plot = elb_p,
       width = 7,
       height = 5,
       dpi = 300)