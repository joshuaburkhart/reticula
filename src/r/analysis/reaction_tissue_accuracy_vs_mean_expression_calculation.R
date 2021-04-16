set.seed(88888888)

library(magrittr)

OUT_DIR <- "/home/burkhart/Software/reticula/data/aim1/output/"
toi_summary.df <- readRDS(file=paste(OUT_DIR,"toi_summary_df.Rds",sep="")) # reaction accuracy X tissue
vst_count.df <- readRDS(file=paste(OUT_DIR,"vst_count_mtx_train.Rds",sep="")) # transcript expression value X sample
gtex_tissue_detail.vec <- readRDS(file=paste(OUT_DIR,"gtex_tissue_detail_vec_train.Rds",sep="")) # tissue labels for samples
rxn2ensembl.nls <- readRDS(file=paste(OUT_DIR,"rxn2ensembls_nls.Rds",sep="")) # reaction -> transcript list mapping

tissuewise_dif_mean_rxn_exp.nls <- list()

# calculate difference of mean reaction expression for tissue (mean tissue of interest expression - mean other tissues expression)
dif_mean_rxn_exp <- function(reaction_id,tissue_label){
  tissue_sample_idxs <- which(gtex_tissue_detail.vec == tissue_label)
  rxn_ensembl_ids <- rxn2ensembl.nls[[reaction_id]]
  rxn_ensembl_idxs <- which(rownames(vst_count.df) %in% c(rxn_ensembl_ids))
  mean_rxn_exp_toi <- mean(vst_count.df[rxn_ensembl_idxs,tissue_sample_idxs])
  mean_rxn_exp_other <- mean(vst_count.df[rxn_ensembl_idxs,-tissue_sample_idxs])
  return(mean_rxn_exp_toi - mean_rxn_exp_other)
}

for(tissue_label in unique(gtex_tissue_detail.vec)){
  print(paste("Calculating results for ",tissue_label,"...",sep=""))
  for(reaction_id in names(rxn2ensembl.nls)){
    tissuewise_dif_mean_rxn_exp.nls[[tissue_label]][[reaction_id]] <- dif_mean_rxn_exp(reaction_id,tissue_label)
  }
}
