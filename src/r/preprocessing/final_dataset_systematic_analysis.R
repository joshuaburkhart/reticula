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
gtex_tissue_detail.vec.train <- readRDS(paste(OUT_DIR,"gtex_tissue_detail_vec_train.Rds",sep=""))
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

#testis (64)
#breast (47)
#ovary (14)
#uterus (13)
#cervix (5.7)
#vagina (3.9)

#heart (14800)
#muscle (5510)
#adipose (2448)
#thyroid (3180)
#adrenal (455)

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


high_prolif <- c("Stomach",
            "Colon - Sigmoid",
            "Colon - Transverse",
            "Small Intestine - Terminal Ileum",
            "Spleen",
            "Esophagus - Mucosa",
            "Esophagus - Gastroesophageal Junction")
med_prolif <- c("Bladder",
                "Skin - Not Sun Exposed (Suprapubic)",
                "Lung",
                "Liver",
                "Pancreas",
                "Kidney - Cortex")
low_prolif <- c("Adipose - Visceral (Omentum)",
                "Adipose - Subcutaneous",
                "Thyroid",
                "Muscle - Skeletal",
                "Heart - Left Ventricle",
                "Heart - Atrial Appendage",
                "Adrenal Gland")

# convert initial rxn pca nls to df
rxn_pca.df <- as.data.frame(
  sapply(as.data.frame(
    do.call(rbind, rxn_pca.nls)),
    as.numeric))
rownames(rxn_pca.df) <- names(rxn_pca.nls)
colnames(rxn_pca.df) <- names(rxn_pca.nls[[1]])
rxn_pca.df$RXN_ID <- names(rxn_pca.nls)
rxn_pca.df %>% write.csv(file=paste(OUT_DIR,"rxn_pca.csv",sep=""))

high_prolif_samples <- which(gtex_tissue_detail.vec.train %in% high_prolif)
med_prolif_samples <- which(gtex_tissue_detail.vec.train %in% med_prolif)
low_prolif_samples <- which(gtex_tissue_detail.vec.train %in% low_prolif)

high_v_med_wilcox_res.nls <- list()
for(rxn_idx in seq(1:nrow(rxn_pca.df))){
  w <- wilcox.test(x=as.numeric(rxn_pca.df[rxn_idx,high_prolif_samples]),
                   y=as.numeric(rxn_pca.df[rxn_idx,med_prolif_samples]))
  high_v_med_wilcox_res.nls[[rxn_pca.df$RXN_ID[rxn_idx]]] <- w$p.value
}
med_v_low_wilcox_res.nls <- list()
for(rxn_idx in seq(1:nrow(rxn_pca.df))){
  w <- wilcox.test(x=as.numeric(rxn_pca.df[rxn_idx,med_prolif_samples]),
                   y=as.numeric(rxn_pca.df[rxn_idx,low_prolif_samples]))
  med_v_low_wilcox_res.nls[[rxn_pca.df$RXN_ID[rxn_idx]]] <- w$p.value
}

saveRDS(high_v_med_wilcox_res.nls,file=paste(OUT_DIR,"high_v_med_wilcox_res_nls.Rds",sep=""))
saveRDS(med_v_low_wilcox_res.nls,file=paste(OUT_DIR,"med_v_low_wilcox_res_nls.Rds",sep=""))

high_v_med_wilcox_res.nls <- readRDS(file = paste(OUT_DIR,"high_v_med_wilcox_res_nls.Rds",sep=""))
med_v_low_wilcox_res.nls <- readRDS(file = paste(OUT_DIR,"med_v_low_wilcox_res_nls.Rds",sep=""))

# convert to df
high_v_med_wilcox_res.df <- as.data.frame(
  sapply(as.data.frame(
    do.call(rbind, high_v_med_wilcox_res.nls)),
    as.numeric))
rownames(high_v_med_wilcox_res.df) <- names(high_v_med_wilcox_res.nls)

med_v_low_wilcox_res.df <- as.data.frame(
  sapply(as.data.frame(
    do.call(rbind, med_v_low_wilcox_res.nls)),
    as.numeric))
rownames(med_v_low_wilcox_res.df) <- names(med_v_low_wilcox_res.nls)

saveRDS(high_v_med_wilcox_res.df,file=paste(OUT_DIR,"high_v_med_wilcox_res_df.Rds",sep=""))
high_v_med_wilcox_res.df$fdr <- p.adjust(high_v_med_wilcox_res.df$V1,method = "fdr")
colnames(high_v_med_wilcox_res.df) <- c("Wilcox test p-value","False Discovery Rate")
high_v_med_wilcox_res.df %>% write.csv(file=paste(OUT_DIR,"high_v_med_wilcox_res.csv",sep=""))

saveRDS(med_v_low_wilcox_res.df,file=paste(OUT_DIR,"med_v_low_wilcox_res_df.Rds",sep=""))
med_v_low_wilcox_res.df$fdr <- p.adjust(med_v_low_wilcox_res.df$V1,method = "fdr")
colnames(med_v_low_wilcox_res.df) <- c("Wilcox test p-value","False Discovery Rate")
med_v_low_wilcox_res.df %>% write.csv(file=paste(OUT_DIR,"med_v_low_wilcox_res.csv",sep=""))

combined_wilcox_res.df <- data.frame("rxn_n1" = rownames(high_v_med_wilcox_res.df),
                                     "rxn_n2" = rownames(med_v_low_wilcox_res.df),
                                     "High_v_med_p" = high_v_med_wilcox_res.df$`Wilcox test p-value`,
                                     "Med_v_low_p" = med_v_low_wilcox_res.df$`Wilcox test p-value`)
library(metap)
combined_w_fisher <- combined_wilcox_res.df %>%
  dplyr::rowwise() %>%
  dplyr::mutate(combined_p = as.numeric((metap::sumlog(c(High_v_med_p,Med_v_low_p)) %>% .[3])))

combined_w_fisher$fdr <- p.adjust(combined_w_fisher$combined_p,method = "fdr")

combined_w_fisher %>% write.csv(file=paste(OUT_DIR,"combined_w_fisher.csv",sep=""))

# write initial count matirx df
vst.count.mtx.train <- readRDS(paste(OUT_DIR,"vst_count_mtx_train.Rds",sep=""))
vst.count.mtx.train %>% write.csv(file=paste(OUT_DIR,"vst_count_mtx_train.csv",sep=""))

end_time <- Sys.time()
print(paste(
  "Start: ",
  start_time,
  " End: ",
  end_time,
  " Difference: ",
  end_time - start_time
))
