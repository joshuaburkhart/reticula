set.seed(88888888)

library(magrittr)
library(ggplot2)

OUT_DIR <- "/home/burkhart/Software/reticula/data/aim1/output/"
toi_summary.df <- readRDS(file=paste(OUT_DIR,"toi_summary_df.Rds",sep="")) # reaction accuracy X tissue
vst_count.mtx <- as.matrix(readRDS(file=paste(OUT_DIR,"vst_count_mtx_train.Rds",sep=""))) # transcript expression value X sample
gtex_tissue_detail.vec <- readRDS(file=paste(OUT_DIR,"gtex_tissue_detail_vec_train.Rds",sep="")) # tissue labels for samples
rxn2ensembl.nls <- readRDS(file=paste(OUT_DIR,"rxn2ensembls_nls.Rds",sep="")) # reaction -> transcript list mapping

tissuewise_mean_rxn_exp.nls <- list()

mean_rxn_exp <- function(reaction_id,tissue_sample_idxs){
  rxn_ensembl_ids <- rxn2ensembl.nls[[reaction_id]]
  rxn_ensembl_idxs <- which(rownames(vst_count.mtx) %in% rxn_ensembl_ids)
  return(mean(vst_count.mtx[rxn_ensembl_idxs,tissue_sample_idxs]))
}

for(tissue_label in unique(gtex_tissue_detail.vec)){
  print(paste("Calculating results for ",tissue_label,"...",sep=""))
  tissue_sample_idxs <- which(gtex_tissue_detail.vec == tissue_label)
  for(reaction_id in names(rxn2ensembl.nls)){
    tissuewise_mean_rxn_exp.nls[[tissue_label]][[reaction_id]] <- mean_rxn_exp(reaction_id,tissue_sample_idxs)
  }
}

saveRDS(tissuewise_mean_rxn_exp.nls,file=paste(OUT_DIR,"tissuewise_mean_rxn_nls.Rds",sep = ""))

z <- data.frame(default=numeric(10726)) #need "default" column to use cbind() with non-zero length vector

for(tis_name in names(tissuewise_mean_rxn_exp.nls)){
  old_colnames <- colnames(z)
  z <- cbind(z,unlist(tissuewise_mean_rxn_exp.nls[[tis_name]]))
  colnames(z) <- c(old_colnames,tis_name)
}

z <- z[-1] #remove "default" column

write.csv(z,file=paste(OUT_DIR,"tissuewise_mean_rxn_df.csv",sep = ""))

tissuewise_dif_mean_rxn_exp.nls <- list()

# calculate difference of mean reaction expression for tissue (mean tissue of interest expression - mean other tissues expression)
dif_mean_rxn_exp <- function(reaction_id,tissue_sample_idxs){
  rxn_ensembl_ids <- rxn2ensembl.nls[[reaction_id]]
  rxn_ensembl_idxs <- which(rownames(vst_count.mtx) %in% rxn_ensembl_ids)
  mean_rxn_exp_toi <- mean(vst_count.mtx[rxn_ensembl_idxs,tissue_sample_idxs])
  mean_rxn_exp_other <- mean(vst_count.mtx[rxn_ensembl_idxs,-tissue_sample_idxs])
  return(mean_rxn_exp_toi - mean_rxn_exp_other)
}

for(tissue_label in unique(gtex_tissue_detail.vec)){
  print(paste("Calculating results for ",tissue_label,"...",sep=""))
  tissue_sample_idxs <- which(gtex_tissue_detail.vec == tissue_label)
  for(reaction_id in names(rxn2ensembl.nls)){
    tissuewise_dif_mean_rxn_exp.nls[[tissue_label]][[reaction_id]] <- dif_mean_rxn_exp(reaction_id,tissue_sample_idxs)
  }
}

saveRDS(tissuewise_dif_mean_rxn_exp.nls,file=paste(OUT_DIR,"tissuewise_dif_mean_rxn_nls.Rds",sep = ""))

pearson_correlations <- list()
for(tissue_label in unique(gtex_tissue_detail.vec)){
name_matches <- tissuewise_dif_mean_rxn_exp.nls[[tissue_label]] %>% names() != rownames(toi_summary.df)
name_matches %>% as.numeric() %>% sum() == 0
df <- data.frame(ExprDiff = tissuewise_dif_mean_rxn_exp.nls[[tissue_label]],
                 RxnAccur = 1 - toi_summary.df[,tissue_label])
pearson_correlation <- cor(df$ExprDiff,df$RxnAccur)
pearson_correlations[[tissue_label]] <- pearson_correlation
print(paste(tissue_label," reaction accuracy vs mean expression difference correlation: ",pearson_correlation,sep=""))
ggplot2::ggplot(df, aes(x=ExprDiff,y=RxnAccur)) +
  geom_point() +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle(paste(tissue_label," (pearson=",pearson_correlation,")",sep=""))
ggsave(file=paste(OUT_DIR,tissue_label,"_rxn_acc_vs_mean_expr_dif.png",sep=""))
}

par(mar=c(15,3,1,1))
pearson_correlations %>% unlist() %>% sort() %>% barplot(las=2,cex.names = 0.75)
