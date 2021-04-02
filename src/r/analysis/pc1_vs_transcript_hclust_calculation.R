set.seed(88888888) # maximum luck

library(magrittr)
library(dendextend)
library(parallelDist)

start_time <- Sys.time()

OUT_DIR <- "/home/burkhart/Software/reticula/data/aim1/output/"
#OUT_DIR <- "/Users/burkhajo/Software/reticula/data/aim1/output/"

#gtex_tissue_detail.vec <- readRDS(paste(OUT_DIR,"gtex_tissue_detail_vec.Rds",sep=""))
vst.count.mtx.train <- readRDS(paste(OUT_DIR,"vst_count_mtx_train.Rds",sep=""))
rxn_pca.nls <- readRDS(paste(OUT_DIR,"rxn_pca_nls.Rds",sep=""))

rxn_pca.df <- as.data.frame(
  sapply(as.data.frame(
    do.call(rbind, rxn_pca.nls)),
    as.numeric))
rownames(rxn_pca.df) <- names(rxn_pca.nls)

# hierarchical clustering on reaction 1st pcs
df <- scale(t(rxn_pca.df))
d <- parallelDist::parallelDist(df, method = "euclidean")
saveRDS(d,file=paste(OUT_DIR,"rxn_pca_dist_obj.Rds",sep=""))
hc1 <- hclust(d, method = "ward.D2" )
saveRDS(hc1,file=paste(OUT_DIR,"rxn_pca_hc_obj.Rds",sep=""))

# hierarchical clustering on transcript counts
df_t <- scale(t(vst.count.mtx.train))
d_t <- parallelDist::parallelDist(df_t, method = "euclidean")
saveRDS(d_t,file=paste(OUT_DIR,"transcript_count_dist_obj.Rds",sep=""))
hc2 <- hclust(d_t,method = "ward.D2")
saveRDS(hc2,file=paste(OUT_DIR,"transcript_count_hc_obj.Rds",sep=""))

# R session crashes here. Load hierarchical clustering objects hc1, hc2 from disk to continue

# Plot the obtained dendrogram
#plot(hc1, cex = 0.6, hang = -1)
#plot(hc2, cex = 0.6, hang = -1)

hc1 <- readRDS(paste(OUT_DIR,"rxn_pca_hc_obj.Rds",sep=""))
hc2 <- readRDS(paste(OUT_DIR,"transcript_count_hc_obj.Rds",sep=""))

dend1 <- as.dendrogram(hc1)
dend2 <- as.dendrogram(hc2)

# view the dendrograms with plot(dend1) etc.

dend_list <- dendextend::dendlist(dend1, dend2)

#stable across pearson/spearman methods
cor <- dendextend::cor_cophenetic(dend1,dend2)

print(paste("Cophenetic correlation = ",cor,".",sep=""))

THIS_WORKS = FALSE
if(THIS_WORKS){
  library(foreach)
  library(doParallel)

  cores = detectCores()
  n_cl_cores <- cores[1] - 2
  print(n_cl_cores)
  cl = makeCluster(n_cl_cores)
  registerDoParallel(cl)
  
  # and then change the below for loop into a foreach with %dopar%
}

permutation_rxn_correlations <- numeric()
n_permutations <- 10000 #10K should be enough
for(i in 1:n_permutations){
  hc_cur <- hc1
  hc_cur$labels <- sample(hc_cur$labels)
  res <- dendextend::cor_cophenetic(as.dendrogram(hc_cur), dend2)
  permutation_rxn_correlations <- c(permutation_rxn_correlations,
                                    res)
  print(paste("Permutation ",i," of ",n_permutations,"(",round((i/n_permutations)*100,digits = 4),"%). Correlation: ",res,". Max: ",max(permutation_rxn_correlations),sep=""))
}

saveRDS(permutation_rxn_correlations,file=paste(OUT_DIR,"permutation_correlations_vec.Rds",sep=""))

# from North BV, Curtis D, Sham PC. A note on the calculation of empirical P values from Monte Carlo procedures. The American Journal of Human Genetics. 2002 Aug 1;71(2):439-41.

r <- sum(cor <= permutation_rxn_correlations)
n <- length(permutation_rxn_correlations)

permutation_pvalue <- (r + 1)/(n + 1)

print(paste("permutation p-value:",formatC(permutation_pvalue, format = "e", digits = 2)))

end_time <- Sys.time()
print(paste(
  "Start: ",
  start_time,
  " End: ",
  end_time,
  " Difference: ",
  end_time - start_time
))
