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
rxn_pca.nls <- readRDS(paste(OUT_DIR, "rxn_pca_nls.Rds", sep = ""))

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

# create pca-pca df

# clustering on rxn_pca.nls (review compare_toi_counts_to_pca.R for extracting contributions)

# hierarchical clustering, compare with raw transcript counts

end_time <- Sys.time()
print(paste(
  "Start: ",
  start_time,
  " End: ",
  end_time,
  " Difference: ",
  end_time - start_time
))