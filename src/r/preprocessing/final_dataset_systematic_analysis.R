set.seed(88888888) # maximum luck

library(magrittr)
library(ggplot2)
library(ggiraph)

start_time <- Sys.time()

OUT_DIR <- "/Users/burkhajo/Software/reticula/data/aim1/output/"#"/home/burkhart/Software/reticula/data/aim1/output/"

rxn2ensembls.nls <- readRDS(paste(OUT_DIR, "rxn2ensembls_nls.Rds", sep = ""))
rxn_knn_misclass_rate.nls <- readRDS(paste(OUT_DIR, "toi_rxn_knn_misclass_rate_nls.Rds", sep = ""))
rxn_knn_ari.nls <- readRDS(paste(OUT_DIR, "toi_rxn_knn_ari_nls.Rds", sep = ""))
rxn_knn_ecount.nls <- readRDS(paste(OUT_DIR, "toi_rxn_knn_ecount_nls.Rds", sep = ""))
gtex_tissue_detail.vec <- readRDS(paste(OUT_DIR,"gtex_tissue_detail_vec.Rds",sep=""))
vst.count.mtx.train <- readRDS(paste(OUT_DIR,"vst_count_mtx_train.Rds",sep=""))
rxn_pca.nls <- readRDS(paste(OUT_DIR,"rxn_pca_nls.Rds",sep=""))

# construct summary data frame
rxn_tissue_mean_misclass.df <- as.data.frame(
                                      sapply(as.data.frame(
                                                   do.call(rbind, rxn_knn_misclass_rate.nls)),
                                             as.numeric))
rownames(rxn_tissue_mean_misclass.df) <- names(rxn_knn_misclass_rate.nls)
rxn_tissue_mean_misclass.df$RXN_ID <- names(rxn_knn_misclass_rate.nls)
rxn_tissue_mean_misclass.df$ARI <- unlist(rxn_knn_ari.nls)
rxn_tissue_mean_misclass.df$ECOUNT <- unlist(rxn_knn_ecount.nls)

# store summary data frame
saveRDS(rxn_tissue_mean_misclass.df, paste(OUT_DIR, "toi_summary_df.Rds", sep = ""))

# generate figures using summary data frame
for(tis_idx in seq(1:51)){
  tis_name <- colnames(rxn_tissue_mean_misclass.df) %>% .[tis_idx]
  sorted.df <- rxn_tissue_mean_misclass.df %>% dplyr::arrange(ECOUNT)

  plot.obj <- ggplot2::ggplot(sorted.df) + 
    ggiraph::geom_point_interactive(aes(x=ARI,
                                      y=1 - !!as.name(tis_name),
                                      colour=ECOUNT,
                                      tooltip=RXN_ID,
                                      data_id=RXN_ID)) +
    theme_bw() + 
    ggtitle(paste("ARI vs ",tis_name," 1 - misclassification rate",sep=""))

  #girafe(ggobj = plot.obj)
  
  ggsave(paste(OUT_DIR,"ARI_v_",tis_name,"_misclassification.png"),device = png())  
  dev.off()
}
  
# top n reactions for each tissue
n <- 10
for(tis_idx in seq(1:51)){
  tis_name <- colnames(rxn_tissue_mean_misclass.df) %>% .[tis_idx]
  sorted.df <- rxn_tissue_mean_misclass.df %>% dplyr::arrange(!!as.name(tis_name)) %>%
    dplyr::slice(1:n)
  d <- data.frame(RXN_ID = sorted.df$RXN_ID,
                  TIS = sorted.df[tis_idx],
                  ARI = sorted.df$ARI)
  write.csv(d,file=paste(OUT_DIR,"top_",n,"_",tis_name,"_rxns.csv",sep=""))
}

# compare wilcoxon rank sums for each reaction across proliferative & non-proliferative tissues
prolif <- c("Liver",
            "Lung",
            "Breast - Mammary Tissue",
            "Colon - Transverse",
            "Colon - Sigmoid",
            "Esophagus - Mucosa",
            "Esophagus - Gastroesophageal Junction")
non_prolif <- c("Brain - Substantia nigra",
                "Brain - Nucleus accumbens (basal ganglia)",
                "Brain - Cortex",
                "Muscle - Skeletal",
                "Heart - Left Ventricle",
                "Heart - Atrial Appendage",
                "Nerve - Tibial")
wilcox_res.nls <- list()
for(rxn_idx in seq(1:nrow(rxn_tissue_mean_misclass.df))){
  w <- wilcox.test(x=as.numeric(rxn_tissue_mean_misclass.df[rxn_idx,prolif]),
                   y=as.numeric(rxn_tissue_mean_misclass.df[rxn_idx,non_prolif]))
  wilcox_res.nls[[rxn_tissue_mean_misclass.df$RXN_ID[rxn_idx]]] <- w$p.value
}

saveRDS(wilcox_res.nls,file=paste(OUT_DIR,"wilcox_res_nls.Rds",sep=""))
