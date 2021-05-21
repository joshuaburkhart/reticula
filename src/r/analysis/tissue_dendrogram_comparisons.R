OUT_DIR <- "/home/burkhart/Software/reticula/data/aim1/output/"
IN_DIR <- "/home/burkhart/Software/reticula/data/aim2/input/"

# tissue dendrogram comparisons

# misclass
#    df <- scale(t(misclass_only.df))
#    d <- parallelDist::parallelDist(df, method = "euclidean")
#    saveRDS(d,file=paste(OUT_DIR,"misclass_dist_obj.Rds",sep=""))
#    hc1 <- hclust(d, method = "ward.D2" )
#    saveRDS(hc1,file=paste(OUT_DIR,"misclass_hc_obj.Rds",sep=""))
hc1 <- readRDS(paste(OUT_DIR,"misclass_hc_obj.Rds",sep=""))
dend1 <- as.dendrogram(hc1)
plot(hc1, cex = .6)

# reaction network edge weights
rxn_edge_weights.df <- read.csv(paste(IN_DIR,"labelled_edge_weights.csv",sep=""))

ig <- rxn_edge_weights.df[,6:56]
df <- scale(t(ig))
d <- parallelDist::parallelDist(df,method="euclidean")
hc1 <- hclust(d,method="ward.D2")
dend1 <- as.dendrogram(hc1)
plot(hc1,cex=.6)

sal <- rxn_edge_weights.df[,57:107]
df <- scale(t(sal))
d <- parallelDist::parallelDist(df,method="euclidean")
hc1 <- hclust(d,method="ward.D2")
dend1 <- as.dendrogram(hc1)
plot(hc1,cex=.6)

# pathway hierarchy edge weights
ph_edge_weights.df <- read.csv(paste(IN_DIR,"pathway_hierarchy_labelled_edge_weights.csv",sep=""))
ig <- ph_edge_weights.df[,6:56]
df <- scale(t(ig))
d <- parallelDist::parallelDist(df,method="euclidean")
hc1 <- hclust(d,method="ward.D2")
dend1 <- as.dendrogram(hc1)
plot(hc1,cex=.6)

sal <- ph_edge_weights.df[,57:107]
df <- scale(t(sal))
d <- parallelDist::parallelDist(df,method="euclidean")
hc1 <- hclust(d,method="ward.D2")
dend1 <- as.dendrogram(hc1)
plot(hc1,cex=.6)

# reaction pc1 (two-tailed wilcoxon)
mean_rxn_pc1.df <- read.csv(paste(OUT_DIR,"tissuewise_mean_rxn_df.csv",sep=""))
df <- scale(t(mean_rxn_pc1.df))
d <- parallelDist::parallelDist(df,method="euclidean")
hc1 <- hclust(d,method="ward.D2")
dend1 <- as.dendrogram(hc1)
plot(hc1,cex=.6)

# update below for tissue dendrogram comparisons
hc1 <- readRDS(paste(OUT_DIR,"rxn_pca_hc_obj.Rds",sep=""))
hc2 <- readRDS(paste(OUT_DIR,"transcript_count_hc_obj.Rds",sep=""))

dend1 <- as.dendrogram(hc1)
dend2 <- as.dendrogram(hc2)

ord.h1 <- order.hclust(hc1)
ord.h2 <- order.hclust(hc2)

ord.d1 <- order.dendrogram(dend1)
ord.d2 <- order.dendrogram(dend2)

sum(ord.h1 == ord.d1)
sum(ord.h2 == ord.d2)

ord.df <- df[ord.h1,]
ord.df_t <- df_t[ord.h2,]

ord.df.downsampled <- ord.df[seq(from=1,to=nrow(ord.df),by=50),
                             seq(from=1,to=ncol(ord.df),by=100)]
ord.df_t.downsampled <- ord.df_t[seq(from=1,to=nrow(ord.df_t),by=50),
                                 seq(from=1,to=ncol(ord.df_t),by=100)]

ord.df.downsampled %>% dim()
ord.df_t.downsampled %>% dim()

library(pheatmap)
df.dcols = dist(t(ord.df.downsampled), method = "minkowski")
df_t.dcols = dist(t(ord.df_t.downsampled), method="minkowski")
callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}
pheatmap(ord.df.downsampled,
         border_color = NA,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         clustering_distance_cols = df.dcols,
         legend = FALSE,
         show_rownames = FALSE,
         show_colnames = FALSE)
pheatmap(ord.df_t.downsampled,
         border_color = NA,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         clustering_distance_cols = df_t.dcols,
         legend = FALSE,
         show_rownames = FALSE,
         show_colnames = FALSE)
#library(plot.matrix)

#plot(ord.df)
#plot(ord.df_t)

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

permutation_pvalue <- (r + 1)/(n + 1) # exact monte carlo p-value from https://arxiv.org/pdf/1603.05766.pdf (merely additional citation...see reference above)

print(paste("permutation p-value:",formatC(permutation_pvalue, format = "e", digits = 2)))