
library(DESeq2)
library(ggplot2)
library(viridis)
library(magrittr)
library(pheatmap)
library(DescTools)
library(pdfCluster)
library(RColorBrewer)
library(SummarizedExperiment)
library(caret)
library(class)

start_time <- Sys.time()

OUT_DIR <- "/home/burkhart/Software/reticula/data/aim1/output/"

#DESeq2 PCA plot
vst.counts <- readRDS(paste(OUT_DIR,"vst_counts.Rds",sep=""))

vst.count.mtx <- vst.counts %>% SummarizedExperiment::assay() %>% as.data.frame()
gtex_tissue_detail.vec <- readRDS(paste(OUT_DIR,"gtex_tissue_detail_vec.Rds",sep=""))
rxn2ensembls.nls <- readRDS(paste(OUT_DIR,"rxn2ensembls_nls.Rds",sep=""))
rxns <- rxn2ensembls.nls %>% names()

#calculate clustering coefficents for each reaction
rxn_knn_ari.nls <- list()
count <- 0

training_indices <- caret::createDataPartition(gtex_tissue_detail.vec,
                           times = 1,
                           p = 0.8,
                           list = FALSE)

vst.count.mtx.train <- vst.count.mtx[training_indices,]
vst.count.mtx.test  <- vst.count.mtx[-training_indices,]

gtex_tissue_detail.vec.train <- gtex_tissue_detail.vec[training_indices]
gtex_tissue_detail.vec.test <- gtex_tissue_detail.vec[-training_indices]

for(rxn_id in rxns){
 ensembl_ids <- rxn2ensembls.nls[[rxn_id]]
 train.expr_mat <- t(vst.count.train[ensembl_ids,])
 test.expr_mat <- t(vst.count.test[ensembl_ids,])
 
 rxn_knn_calls <- class::knn(train = train.expr_mat,
                           test = test.expr_mat,
                           cl = gtex_tissue_detail.vec.train)
   
 # calculate & store adjusted rand index
 cur_ari <- pdfCluster::adj.rand.index(rxn_knn_calls,
                                     gtex_tissue_detail.vec.test)
 rxn_knn_ari.nls[[rxn_id]] <- cur_ari
 
 #store ensembl transcript count
 ecount <- length(ensembl_ids)
 
 count <- count + 1
 if(mod(count,10) == 0){
   print(paste("Last ARI:ECOUNT = ",cur_ari,":",ecount,". Clustered ",count," of ",length(rxns)," reactions..."))
  flush.console()
 }
}

saveRDS(rxn_knn_ari.nls,paste(OUT_DIR,"rxn_knn_ari_nls.Rds",sep=""))

end_time <- Sys.time()
print(paste("Start: ",start_time," End: ",end_time," Difference: ",end_time - start_time))