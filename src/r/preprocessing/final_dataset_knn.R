
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
N_FOLDS <- 5

#DESeq2 PCA plot
vst.counts <- readRDS(paste(OUT_DIR,"vst_counts.Rds",sep=""))

vst.count.mtx <- vst.counts %>% SummarizedExperiment::assay() %>% as.data.frame()
gtex_tissue_detail.vec <- readRDS(paste(OUT_DIR,"gtex_tissue_detail_vec.Rds",sep=""))
rxn2ensembls.nls <- readRDS(paste(OUT_DIR,"rxn2ensembls_nls.Rds",sep=""))
rxns <- rxn2ensembls.nls %>% names()

#calculate clustering coefficents for each reaction
rxn_knn_misclass_rate.nls <- list()
rxn_knn_ari.nls <- list()
count <- 0

# take a look at colon samples...
gtex_tissue_detail_vec_colon_s_t_only <- gtex_tissue_detail.vec[which(gtex_tissue_detail.vec == "Colon - Sigmoid" |
                                                                      gtex_tissue_detail.vec == "Colon - Transverse")]

training_indices <- caret::createDataPartition(gtex_tissue_detail_vec_colon_s_t_only,
                           times = 1,
                           p = 0.8,
                           list = FALSE)

vst.count.mtx.train <- vst.count.mtx[,training_indices] #4/5ths of data
vst.count.mtx.test  <- vst.count.mtx[,-training_indices] #1/5th of data

gtex_tissue_detail.vec.train <- gtex_tissue_detail_vec_colon_s_t_only[training_indices]
gtex_tissue_detail.vec.test <- gtex_tissue_detail_vec_colon_s_t_only[-training_indices]

cv_fold_indices <- caret::createFolds(gtex_tissue_detail.vec.train,
                                      k = N_FOLDS)

for(rxn_id in rxns){
 ensembl_ids <- rxn2ensembls.nls[[rxn_id]]
 
 sum_misclass_rate <- 0
 sum_ari <- 0
 
 for(cv_fold in names(cv_fold_indices)){
   
   cur_cv_fold_indices <- cv_fold_indices[[cv_fold]] 
    
   vst.count.mtx.train.cv_train <- vst.count.mtx.train[,-cur_cv_fold_indices] # 4/5ths of training data
   vst.count.mtx.train.cv_test <- vst.count.mtx.train[,cur_cv_fold_indices] #1/5th of training data
   
   gtex_tissue_detail.vec.train.cv_train <- gtex_tissue_detail.vec.train[-cur_cv_fold_indices]
   gtex_tissue_detail.vec.train.cv_test <- gtex_tissue_detail.vec.train[cur_cv_fold_indices]
   
   cv_train.expr_mat <- t(vst.count.mtx.train.cv_train[ensembl_ids,])
   cv_test.expr_mat <- t(vst.count.mtx.train.cv_test[ensembl_ids,])
 
   rxn_knn_calls <- class::knn(train = cv_train.expr_mat,
                           test = cv_test.expr_mat,
                           cl = gtex_tissue_detail.vec.train.cv_train)
   
   # calculate misclassification rate for "Colon - Transverse" and "Colon - Sigmoid" (https://stat.ethz.ch/pipermail/r-help/2011-September/288885.html)
   tab <- table(rxn_knn_calls,gtex_tissue_detail.vec.train.cv_test)
   cur_misclass_rate <- 1-sum(diag(tab))/sum(tab)
   
   # calculate & store adjusted rand index
   cur_ari <- pdfCluster::adj.rand.index(rxn_knn_calls,
                                     gtex_tissue_detail.vec.train.cv_test)
   
   sum_misclass_rate <- cur_misclass_rate + sum_misclass_rate
   sum_ari <- cur_ari + sum_ari
 }
 
 mean_misclass_rate <- sum_misclass_rate / N_FOLDS
 mean_ari <- sum_ari / N_FOLDS
 
 rxn_knn_misclass_rate.nls[[rxn_id]] <- mean_misclass_rate
 rxn_knn_ari.nls[[rxn_id]] <- mean_ari
 
 #store ensembl transcript count
 ecount <- length(ensembl_ids)
 
 count <- count + 1
 if(mod(count,10) == 0){
   print(paste("Last RXN_ID = ",rxn_id,
               ": Last MISCLASS = ",mean_misclass_rate,
               ": Last ARI = ",mean_ari,
               ": Last ECOUNT = ",ecount,
               ". Clustered ",count," of ",length(rxns)," reactions..."))
  flush.console()
 }
}

saveRDS(rxn_knn_misclass_rate.nls,paste(OUT_DIR,"colon_rxn_knn_misclass_rate_nls.Rds",sep=""))
saveRDS(rxn_knn_ari.nls,paste(OUT_DIR,"colon_rxn_knn_ari_nls.Rds",sep=""))

end_time <- Sys.time()
print(paste("Start: ",start_time," End: ",end_time," Difference: ",end_time - start_time))