set.seed(88888888) # maximum luck

library(magrittr)
library(parallelDist)

start_time <- Sys.time()

OUT_DIR <- "/Users/burkhajo/Software/reticula/data/aim1/output/"#"/home/burkhart/Software/reticula/data/aim1/output/"

gtex_tissue_detail.vec <- readRDS(paste(OUT_DIR,"gtex_tissue_detail_vec.Rds",sep=""))
vst.count.mtx.train <- readRDS(paste(OUT_DIR,"vst_count_mtx_train.Rds",sep=""))
rxn_pca.nls <- readRDS(paste(OUT_DIR,"rxn_pca_nls.Rds",sep=""))

rxn_pca.df <- as.data.frame(
  sapply(as.data.frame(
    do.call(rbind, rxn_pca.nls)),
    as.numeric))
rownames(rxn_pca.df) <- names(rxn_pca.nls)

df <- scale(rxn_pca.df)
d <- parallelDist::parallelDist(df, method = "euclidean")
saveRDS(d,file=paste(OUT_DIR,"rxn_pca_dist_obj.Rds",sep=""))
hc1 <- hclust(d, method = "ward.D2" )
saveRDS(hc1,file=paste(OUT_DIR,"rxn_pca_hc_obj.Rds",sep=""))

df_t <- scale(vst.count.mtx.train)
d_t <- parallelDist::parallelDist(df_t, method = "euclidean")
saveRDS(d_t,file=paste(OUT_DIR,"transcript_count_dist_obj.Rds",sep=""))
hc2 <- hclust(d_t,method = "ward.D2")
saveRDS(hc2,file=paste(OUT_DIR,"transcript_count_hc_obj.Rds",sep=""))

# Plot the obtained dendrogram
#plot(hc1, cex = 0.6, hang = -1)
#plot(hc2, cex = 0.6, hang = -1)

end_time <- Sys.time()
print(paste(
  "Start: ",
  start_time,
  " End: ",
  end_time,
  " Difference: ",
  end_time - start_time
))