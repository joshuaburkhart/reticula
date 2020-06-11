
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

#DESeq2 PCA plot
vst.counts <- readRDS(paste(OUT_DIR,"vst_counts.Rds",sep=""))

vst.count.mtx <- vst.counts %>% SummarizedExperiment::assay() %>% as.data.frame()
gtex_tissue_detail.vec <- readRDS(paste(OUT_DIR,"gtex_tissue_detail_vec.Rds",sep=""))
rxn2ensembls.nls <- readRDS(paste(OUT_DIR,"rxn2ensembls_nls.Rds",sep=""))
rxns <- rxn2ensembls.nls %>% names()

#calculate clustering coefficents for each reaction
rxn_ari.nls <- list()
rxn_gini.nls <- list()
rxn_ensembl_counts.nls <- list()
count <- 0
for(rxn_id in rxns){
 ensembl_ids <- rxn2ensembls.nls[[rxn_id]]
 # combined.df columns ordered by clustering
 expr_matrix <- vst.count.mtx[ensembl_ids,]
 km_obj <- kmeans(t(expr_matrix),centers = 54)
 
 km_cluster_calls <- km_obj$cluster
 km_clust_labels <- km_obj$cluster %>% unique()
 for(km_clust_label in km_clust_labels){
         km_clust_label_gtex_annotations <- gtex_tissue_detail.vec[which(km_obj$cluster %in% km_clust_label)]
         majority_annotation <- km_clust_label_gtex_annotations %>% table() %>% .[1] %>% names()
         km_cluster_calls[which(km_cluster_calls %in% km_clust_label)] <- majority_annotation
 }
 unnamed_km_cluster_calls <- unname(km_cluster_calls)
 
 # calculate & store adjusted rand index
 ari <- pdfCluster::adj.rand.index(unnamed_km_cluster_calls,
                               gtex_tissue_detail.vec)
 rxn_ari.nls[[rxn_id]] <- ari
 
 #calculate & store gini index
 gini <- DescTools::Gini(unnamed_km_cluster_calls)
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

saveRDS(rxn_ari.nls,paste(OUT_DIR,"rxn_ari_nls.Rds",sep=""))
saveRDS(rxn_gini.nls,paste(OUT_DIR,"rxn_gini_nls.Rds",sep=""))
saveRDS(rxn_ensembl_counts.nls,paste(OUT_DIR,"rxn_ensembl_counts.nls",sep=""))

end_time <- Sys.time()
print(paste("Start: ",start_time," End: ",end_time," Difference: ",end_time - start_time))