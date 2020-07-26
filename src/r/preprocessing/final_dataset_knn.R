set.seed(88888888) # maximum luck

library(DESeq2)
library(plotly)
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
rxn_knn_ecount.nls <- list()
rxn_pca.nls <- list()
count <- 0

# gi tract ->
#"Esophagus - Mucosa"
#"Esophagus - Gastroesophageal Junction"
#"Stomach"
#"Colon - Transverse"
#"Colon - Sigmoid"

# muscle ->
#"Esophagus - Muscularis"
#"Heart - Atrial Appendage"
#"Heart - Left Ventricle"
#"Muscle - Skeletal"

# common cancers ->
#"Breast - Mammary Tissue"
#"Lung"
#"Prostate"

# sanity check ->
#"Brain - Cerebellum"
#"Muscle - Skeletal"

# take a look at toi samples...
toi_indices <- which(gtex_tissue_detail.vec == "Colon - Transverse" |
                     gtex_tissue_detail.vec == "Colon - Sigmoid")

# filter annotations
gtex_tissue_detail_vec_tis_of_interest <- gtex_tissue_detail.vec[toi_indices]

# filter expression data
vst.count.mtx.tis_of_interest <- vst.count.mtx[,toi_indices]

training_indices <- caret::createDataPartition(gtex_tissue_detail_vec_tis_of_interest,
                           times = 1,
                           p = 0.9,
                           list = FALSE)

vst.count.mtx.train <- vst.count.mtx.tis_of_interest[,training_indices] #9/10ths of data
vst.count.mtx.test  <- vst.count.mtx.tis_of_interest[,-training_indices] #1/10th of data

gtex_tissue_detail.vec.train <- gtex_tissue_detail_vec_tis_of_interest[training_indices]
gtex_tissue_detail.vec.test <- gtex_tissue_detail_vec_tis_of_interest[-training_indices]

cv_fold_indices <- caret::createFolds(gtex_tissue_detail.vec.train,
                                      k = N_FOLDS)

for(rxn_id in rxns){
 ensembl_ids <- rxn2ensembls.nls[[rxn_id]]
 
 sum_misclass_rate <- 0
 sum_ari <- 0
 rxn_pca <- prcomp(t(vst.count.mtx.train[rxn2ensembls.nls[[rxn_id]],]),scale. = T)
 rxn_pca.nls[[rxn_id]] <- rxn_pca$x[,1] # 1st principal component of this reaction for each sample
 
 for(cv_fold in names(cv_fold_indices)){
   
   cur_cv_fold_indices <- cv_fold_indices[[cv_fold]] 
    
   vst.count.mtx.train.cv_train <- vst.count.mtx.train[,-cur_cv_fold_indices] # 4/5ths of training data
   vst.count.mtx.train.cv_test <- vst.count.mtx.train[,cur_cv_fold_indices] # 1/5th of training data
   
   gtex_tissue_detail.vec.train.cv_train <- gtex_tissue_detail.vec.train[-cur_cv_fold_indices]
   gtex_tissue_detail.vec.train.cv_test <- gtex_tissue_detail.vec.train[cur_cv_fold_indices]
   
   cv_train.expr_mat <- t(vst.count.mtx.train.cv_train[ensembl_ids,])
   cv_test.expr_mat <- t(vst.count.mtx.train.cv_test[ensembl_ids,])
 
   rxn_knn_calls <- class::knn(train = cv_train.expr_mat,
                           test = cv_test.expr_mat,
                           cl = gtex_tissue_detail.vec.train.cv_train)
   
   # calculate misclassification rate (https://stat.ethz.ch/pipermail/r-help/2011-September/288885.html)
   tab <- table(rxn_knn_calls,
                gtex_tissue_detail.vec.train.cv_test)
   cur_misclass_rate <- 1-sum(diag(tab))/sum(tab)
   
   # calculate & store adjusted rand index
   cur_ari <- pdfCluster::adj.rand.index(rxn_knn_calls,
                                         gtex_tissue_detail.vec.train.cv_test)
   
   sum_misclass_rate <- cur_misclass_rate + sum_misclass_rate
   sum_ari <- cur_ari + sum_ari
 }
 
 mean_misclass_rate <- sum_misclass_rate / N_FOLDS
 mean_ari <- sum_ari / N_FOLDS
 ecount <- length(ensembl_ids)
 
 rxn_knn_misclass_rate.nls[[rxn_id]] <- mean_misclass_rate
 rxn_knn_ari.nls[[rxn_id]] <- mean_ari
 rxn_knn_ecount.nls[[rxn_id]] <- ecount
 
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

saveRDS(rxn_knn_misclass_rate.nls,paste(OUT_DIR,"toi_rxn_knn_misclass_rate_nls.Rds",sep=""))
saveRDS(rxn_knn_ari.nls,paste(OUT_DIR,"toi_rxn_knn_ari_nls.Rds",sep=""))
saveRDS(rxn_knn_ecount.nls,paste(OUT_DIR,"toi_rxn_knn_ecount_nls.Rds",sep=""))

# compare informaction content of below files with pca plots or similar
saveRDS(rxn_pca.nls,paste(OUT_DIR,"rxn_pca_nls.Rds",sep=""))
saveRDS(vst.count.mtx.train,paste(OUT_DIR,"vst_count_mtx_train.Rds",sep=""))

d <- data.frame(RXN_ID = names(rxn2ensembls.nls),
                MISCLASS = unlist(rxn_knn_misclass_rate.nls),
                ARI = unlist(rxn_knn_ari.nls),
                ECOUNT = unlist(rxn_knn_ecount.nls))

saveRDS(d,paste(OUT_DIR,"toi_summary_df.Rds",sep=""))

# generate sample figures
min_misclass_pca <- prcomp(t(vst.count.mtx.train[rxn2ensembls.nls[["R-HSA-983147"]],]),scale. = T)
pca.d <- data.frame(PC1 = min_misclass_pca$x[,1],
                    PC2 = min_misclass_pca$x[,2],
                    PC3 = min_misclass_pca$x[,3],
                    Section = gtex_tissue_detail.vec.train)
ggplot(pca.d) +
   geom_point(aes(x=PC1,y=PC2,colour=Section)) +
   theme_bw()
plot_ly(x=pca.d$PC1, y=pca.d$PC2, z=pca.d$PC3, type="scatter3d", mode="markers", color=pca.d$Section, size = 1)

end_time <- Sys.time()
print(paste("Start: ",start_time," End: ",end_time," Difference: ",end_time - start_time))