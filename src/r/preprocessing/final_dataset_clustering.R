
library(DESeq2)
library(ggplot2)
library(viridis)
library(magrittr)
library(pheatmap)
library(DescTools)
library(pdfCluster)
library(RColorBrewer)
library(SummarizedExperiment)

start_time <- Sys.time()

OUT_DIR <- "/home/burkhart/Software/reticula/data/aim1/output/"

kmeans_maj_vote <- function(expr_matrix,k){ 
  km_obj <- kmeans(expr_matrix,centers = k)

  km_cluster_calls <- km_obj$cluster
  km_clust_labels <- km_obj$cluster %>% unique()
  for(km_clust_label in km_clust_labels){
    km_clust_label_gtex_annotations <- gtex_tissue_detail.vec[which(km_obj$cluster %in% km_clust_label)]
    majority_annotation <- km_clust_label_gtex_annotations %>% table() %>% .[1] %>% names()
    km_cluster_calls[which(km_cluster_calls %in% km_clust_label)] <- majority_annotation
  }
  return(unname(km_cluster_calls))
}

#DESeq2 PCA plot
vst.counts <- readRDS(paste(OUT_DIR,"vst_counts.Rds",sep=""))

vst.count.mtx <- vst.counts %>% SummarizedExperiment::assay() %>% as.data.frame()
gtex_tissue_detail.vec <- readRDS(paste(OUT_DIR,"gtex_tissue_detail_vec.Rds",sep=""))
rxn2ensembls.nls <- readRDS(paste(OUT_DIR,"rxn2ensembls_nls.Rds",sep=""))
rxns <- rxn2ensembls.nls %>% names()

#calculate clustering coefficents for each reaction
rxn_km_cluster_calls.nls <- list()
rxn_ari.nls <- list()
rxn_gini.nls <- list()
rxn_ensembl_counts.nls <- list()
count <- 0
expr_mat <- t(vst.count.mtx[ensembl_ids,])
for(rxn_id in rxns){
 ensembl_ids <- rxn2ensembls.nls[[rxn_id]]
 rxn_km_cluster_calls.nls[[rxn_id]] <- c()
 
 for(k_value in seq(from=31,to=51,by=2)){ # start with number of tissues, odd values up to number of tissue details
   # load km cluster calls
   k_calls <- rxn_km_cluster_calls.nls[[rxn_id]]
   # store km cluster calls
   rxn_km_cluster_calls.nls[[rxn_id]] <- c(k_calls,kmeans_maj_vote(expr_mat,k_value))
 }
 
 # calculate & store adjusted rand index
 ari <- pdfCluster::adj.rand.index(unnamed_km_cluster_calls,
                               gtex_tissue_detail.vec)
 rxn_ari.nls[[rxn_id]] <- ari
 
 # calculate & store gini index
 gini <- DescTools::Gini(as.numeric(as.factor(unnamed_km_cluster_calls)))
 rxn_gini.nls[[rxn_id]] <- gini
 
 #store ensembl transcript count
 ecount <- length(ensembl_ids)
 rxn_ensembl_counts.nls[[rxn_id]] <- ecount
 
 count <- count + 1
 if(mod(count,100) == 0){
  print(paste("Last ARI:GINI:ECOUNT = ",ari,":",gini,":",ecount,". Clustered ",count," of ",length(rxns)," reactions..."))
  flush.console()
 }
}

saveRDS(rxn_km_cluster_calls.nls,paste(OUT_DIR,"rxn_km_cluster_calls.Rds",sep=""))
saveRDS(rxn_ari.nls,paste(OUT_DIR,"rxn_ari_nls.Rds",sep=""))
saveRDS(rxn_gini.nls,paste(OUT_DIR,"rxn_gini_nls.Rds",sep=""))
saveRDS(rxn_ensembl_counts.nls,paste(OUT_DIR,"rxn_ensembl_counts.nls",sep=""))

end_time <- Sys.time()
print(paste("Start: ",start_time," End: ",end_time," Difference: ",end_time - start_time))