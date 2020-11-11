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
# citation: Richardson RB, Allan DS, Le Y. Greater organ involution in highly proliferative tissues associated with the early onset and acceleration of ageing in humans. Experimental gerontology. 2014 Jul 1;55:80-91.

# tissue (mean turnover)
#heart (14800)
#muscle (5510)
#adipose (2448)
#thyroid (3180)
#adrenal (455)

#testis (64)
#breast (47)
#ovary (14)
#uterus (13)
#cervix (5.7)
#vagina (3.9)

#liver (327)
#kidney (270)
#pancreas (265)
#lung (200)
#skin (64)
#bladder (49)
#esophagus (10)
#spleen (7.8)
#small intestine (4)
#colorectum (3.4)
#bone marrow (3.2)
#thymus (2.4)
#stomach (1.4)


prolif <- c("Stomach",
            "Colon - Sigmoid",
            "Colon - Transverse",
            "Small Intestine - Terminal Ileum",
            "Spleen",
            "Esophagus - Mucosa",
            "Esophagus - Gastroesophageal Junction",
            "Bladder",
            "Skin - Not Sun Exposed (Suprapubic)",
            "Lung",
            "Liver",
            "Pancreas",
            "Kidney - Cortex")
non_prolif <- c("Adipose - Visceral (Omentum)",
                "Adipose - Subcutaneous",
                "Thyroid",
                "Muscle - Skeletal",
                "Heart - Left Ventricle",
                "Heart - Atrial Appendage",
                "Adrenal Gland")

wilcox_res.nls <- list()
for(rxn_idx in seq(1:nrow(rxn_tissue_mean_misclass.df))){
  w <- wilcox.test(x=as.numeric(rxn_tissue_mean_misclass.df[rxn_idx,prolif]),
                   y=as.numeric(rxn_tissue_mean_misclass.df[rxn_idx,non_prolif]))
  wilcox_res.nls[[rxn_tissue_mean_misclass.df$RXN_ID[rxn_idx]]] <- w$p.value
}

saveRDS(wilcox_res.nls,file=paste(OUT_DIR,"wilcox_res_nls.Rds",sep=""))

wilcox_res.nls <- readRDS(file = paste(OUT_DIR,"wilcox_res_nls.Rds",sep=""))

# convert to df
wilcox_res.df <- as.data.frame(
  sapply(as.data.frame(
    do.call(rbind, wilcox_res.nls)),
    as.numeric))
rownames(wilcox_res.df) <- names(wilcox_res.nls)

saveRDS(wilcox_res.df,file=paste(OUT_DIR,"wilcox_res_df.Rds",sep=""))
wilcox_res.df$fdr <- p.adjust(wilcox_res.df$V1,method = "fdr")
colnames(wilcox_res.df) <- c("Wilcox test p-value","False Discovery Rate")
wilcox_res.df %>% write.csv(file=paste(OUT_DIR,"wilcox_res.csv",sep=""))
