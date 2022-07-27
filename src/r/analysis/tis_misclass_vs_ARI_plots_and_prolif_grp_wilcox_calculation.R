set.seed(88888888) # maximum luck

library(magrittr)
library(ggplot2)
library(ggiraph)

start_time <- Sys.time()

#OUT_DIR <- "/Users/burkhajo/Software/reticula/data/aim1/output/"
OUT_DIR <- "/home/burkhart/Software/reticula/data/aim1/output/"

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

misclass_only.df <- rxn_tissue_mean_misclass.df[1:51]

# generate dendrogram
df <- scale(t(misclass_only.df))
d <- parallelDist::parallelDist(df, method = "euclidean")
saveRDS(d,file=paste(OUT_DIR,"misclass_dist_obj.Rds",sep=""))
hc1 <- hclust(d, method = "ward.D2" )
saveRDS(hc1,file=paste(OUT_DIR,"misclass_hc_obj.Rds",sep=""))
hc1 <- readRDS(paste(OUT_DIR,"misclass_hc_obj.Rds",sep=""))
dend1 <- as.dendrogram(hc1)
plot(hc1, cex = 2)

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

#plots of classification accuracy - ARI
library(scales)
ari_sorted_rxn_tissue_mean_misclass.df <- rxn_tissue_mean_misclass.df %>% dplyr::arrange(ARI)
ari_sorted_rxn_tissue_mean_misclass.df[,53] %>% scales::rescale() -> x1
plot(x1)
1 - ari_sorted_rxn_tissue_mean_misclass.df %>% .[,36] %>% scales::rescale() -> y1 #36 is breast - mammary tissue
plot(y1)
y1 - x1 -> z1
plot(z1)
hist(z1)

breast_df <- data.frame(rxns = z1)

ggplot(breast_df,aes(x=rxns)) +
  geom_histogram(bins=20) +
  theme_bw()

1 - ari_sorted_rxn_tissue_mean_misclass.df %>% .[,1] %>% scales::rescale() -> y1 #1 is lung
plot(y1)
y1 - x1 -> z1
plot(z1)
hist(z1)

lung_df <- data.frame(rxns = z1)

ggplot(lung_df,aes(x=rxns)) +
  geom_histogram(bins=20) +
  theme_bw()

#ARI histogram
hist(ari_sorted_rxn_tissue_mean_misclass.df$ARI)

rxn_ari_df <- data.frame(ARI = ari_sorted_rxn_tissue_mean_misclass.df$ARI)

ggplot(rxn_ari_df,aes(x=ARI)) +
  geom_smooth(stat = "count") +
  geom_histogram(bins=20) +
  theme_bw()

#rxn log10 ecount histogram
hist(log10(ari_sorted_rxn_tissue_mean_misclass.df$ECOUNT))

rxn_ecount_df <- data.frame(ECOUNT = (ari_sorted_rxn_tissue_mean_misclass.df$ECOUNT))

# plot Transcript Count x Reaction Count
ggplot(rxn_ecount_df,aes(x=ECOUNT)) +
  geom_smooth(stat = "count") +
  geom_bar(stat="count") +
  scale_y_continuous(trans = "log10") +
  scale_x_continuous(trans = "log10") +
  theme_bw()

# plot Reaction Count x Transcript Count
tabl <- rxn_ecount_df %>% table() %>% as.data.frame()
colnames(tabl) <- c("Transcript_Count","Num_Reactions")

ggplot(tabl,aes(x=Num_Reactions,y=Transcript_Count,label=Num_Reactions)) +
  #geom_smooth(stat = "identity") +
  geom_bar(stat="identity",fill = "black") +
  #scale_y_continuous(trans = "log2") +
  geom_point(shape = 23, size = 2.25) +
  geom_text(hjust = -1.5/log10(tabl$Num_Reactions + 10), vjust = 0.5, size = 3) +
  scale_x_continuous(trans = "log10") +
  theme_bw()

#histogram of sample counts across tissue
hist(as.numeric(as.factor(gtex_tissue_detail.vec.train)))

tis_df <- data.frame(Tissue = as.factor(gtex_tissue_detail.vec.train))

ggplot(tis_df,aes(x=reorder(Tissue,desc(Tissue)))) +
  geom_bar(stat = "count") +
  scale_y_discrete(expand = expansion(mult=c(0,.05))) +
  coord_flip() +
  geom_text(stat = "count",
            aes(label=..count..),
            hjust = -0.1) +
  theme_bw() +
  theme(legend.position="none",
        panel.grid.major = element_blank())

# assume the probability as 1 - misclassification rate observed in our data
# (1 - rxn_tissue_mean_misclass.df %>% .[,1:51] %>% colMeans())^51 # "take the probability of correctly classifying each tissue type as the joint mean probability across reactions"
# Lung             Brain - Cerebellar Hemisphere                    Heart - Left Ventricle            Skin - Sun Exposed (Lower leg) 
# 0.05607242                                0.41345513                                0.11555096                                0.03777499 
# Brain - Amygdala                    Adipose - Subcutaneous                            Brain - Cortex                                    Uterus 
# 0.44037486                                0.03604756                                0.32171692                                0.43113388 
# Nerve - Tibial                         Muscle - Skeletal                                     Ovary                               Whole Blood 
# 0.07693189                                0.06433959                                0.41006013                                0.15161890 
# Brain - Nucleus accumbens (basal ganglia)                        Colon - Transverse              Adipose - Visceral (Omentum)                             Adrenal Gland 
# 0.31222191                                0.15767621                                0.11592596                                0.30372533 
# Brain - Anterior cingulate cortex (BA24)                  Brain - Substantia nigra                                   Thyroid                        Esophagus - Mucosa 
# 0.37272218                                0.49630595                                0.06443830                                0.08811945 
# Artery - Coronary                    Esophagus - Muscularis           Brain - Caudate (basal ganglia)                  Heart - Atrial Appendage 
# 0.24283629                                0.07400780                                0.27060816                                0.15548794 
# Esophagus - Gastroesophageal Junction                           Colon - Sigmoid                           Artery - Tibial                                     Liver 
# 0.16336525                                0.18589360                                0.05608241                                0.39626316 
# Prostate                                    Testis           Brain - Putamen (basal ganglia)                           Kidney - Cortex 
# 0.33412688                                0.34127418                                0.35895371                                0.71610065 
# Pancreas                                   Stomach          Small Intestine - Terminal Ileum                   Breast - Mammary Tissue 
# 0.24079190                                0.16166527                                0.38587747                                0.11981344 
# Brain - Hippocampus                        Brain - Cerebellum                                 Pituitary       Skin - Not Sun Exposed (Suprapubic) 
# 0.36472635                                0.35408958                                0.38805811                                0.08493649 
# Brain - Frontal Cortex (BA9)                            Artery - Aorta                       Cervix - Ectocervix                                    Vagina 
# 0.32781540                                0.11555485                                0.93547057                                0.38084504 
# Brain - Hypothalamus                                    Spleen        Brain - Spinal cord (cervical c-1)                      Minor Salivary Gland 
# 0.37914645                                0.42187813                                0.52019008                                0.51332601 
# Bladder                       Cervix - Endocervix                            Fallopian Tube 
# 0.89289252                                0.94669607                                0.92826076 

#(1 - rxn_tissue_mean_misclass.df %>% .[,1:51] %>% colMeans())^51 %>% mean() "then take the mean probability of correctly classifying a novel tissue sample as the mean probability across tissue types"
#[1] 0.3181029

# compare to ARI using:
# > rxn_tissue_mean_misclass.df$ARI %>% mean()
#[1] 0.2116457

# rowMeans() could give us relative reaction probabilities...

# > -log2(0.80)
# [1] 0.3219281
# > -log2(0.3181029)
# [1] 1.652435
# > a <- -log2(0.80)
# > b <- -log2(0.3181029)
# > b - a
# [1] 1.330506
# > (b - a)/40000
# [1] 3.326266e-05
# 
# > (b - a)/0.05
# [1] 26.61013 # at this edgewiese information rate, it may be possible to detect a difference using this method after 27 edges have been added... or not, also PCA$PC1 aggregation was introduced here too

# no weighting, no parameters, just an estimate of probabilities based on this training data using a knn model

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
rxn_pca.df$RXN_ID <- names(rxn_pca.nls) #final column in rxn_pca.df is reaction id
rxn_pca.df %>% write.csv(file=paste(OUT_DIR,"rxn_pca.csv",sep=""))

high_prolif_samples <- which(gtex_tissue_detail.vec.train %in% high_prolif)
med_prolif_samples <- which(gtex_tissue_detail.vec.train %in% med_prolif)
low_prolif_samples <- which(gtex_tissue_detail.vec.train %in% low_prolif)

# compare reaction principal component coordinates
high_v_med_wilcox_res.nls <- list()
med_v_low_wilcox_res.nls <- list()
rxn_pca_direction.nls <- list()
rxn_pca_difference.nls <- list()
transcript_in_rxn.nls <- list()
n_rxn_pca <- nrow(rxn_pca.df)
for(rxn_idx in seq(1:n_rxn_pca)){
  w1 <- wilcox.test(x=as.numeric(rxn_pca.df[rxn_idx,high_prolif_samples]),
                   y=as.numeric(rxn_pca.df[rxn_idx,med_prolif_samples]))
  w2 <- wilcox.test(x=as.numeric(rxn_pca.df[rxn_idx,med_prolif_samples]),
                    y=as.numeric(rxn_pca.df[rxn_idx,low_prolif_samples]))
  pca_direction <- "none"
  mean_hi <- mean(as.numeric(rxn_pca.df[rxn_idx,high_prolif_samples]))
  mean_med <- mean(as.numeric(rxn_pca.df[rxn_idx,med_prolif_samples]))
  mean_low <- mean(as.numeric(rxn_pca.df[rxn_idx,low_prolif_samples]))
  if(mean_hi >= mean_med & mean_med >= mean_low){
    pca_direction <- "positive"
  }else if(mean_low >= mean_med & mean_med >= mean_hi){
    pca_direction <- "negative"
  }
  rxn_id <- rxn_pca.df$RXN_ID[rxn_idx]
  rxn_pca_direction.nls[[rxn_id]] <- pca_direction
  rxn_pca_difference.nls[[rxn_id]] <- mean_hi - mean_low
  transcript_in_rxn.nls[[rxn_id]] <- length(rxn2ensembls.nls[[rxn_id]])
  high_v_med_wilcox_res.nls[[rxn_id]] <- w1$p.value
  med_v_low_wilcox_res.nls[[rxn_id]] <- w2$p.value
  if(mod(rxn_idx,50) == 0){
    print(paste("Processed ",rxn_idx,
                " of ",n_rxn_pca,
                " reactions (",round((rxn_idx + 1)/n_rxn_pca,digits = 3) * 100,"%)...",
                sep=""))
    flush.console()
  }
}
saveRDS(transcript_in_rxn.nls,file=paste(OUT_DIR,"transcript_in_rxn_nls.Rds",sep=""))
saveRDS(rxn_pca_direction.nls,file=paste(OUT_DIR,"rxn_pca_direction_nls.Rds",sep=""))
saveRDS(rxn_pca_difference.nls,file=paste(OUT_DIR,"rxn_pca_difference_nls.Rds",sep=""))
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
high_v_med_wilcox_res.df$fdr <- p.adjust(high_v_med_wilcox_res.df$V1,method = "bonferroni")
colnames(high_v_med_wilcox_res.df) <- c("Wilcox test p-value","False Discovery Rate")
high_v_med_wilcox_res.df %>% write.csv(file=paste(OUT_DIR,"high_v_med_wilcox_res.csv",sep=""))

saveRDS(med_v_low_wilcox_res.df,file=paste(OUT_DIR,"med_v_low_wilcox_res_df.Rds",sep=""))
med_v_low_wilcox_res.df$fdr <- p.adjust(med_v_low_wilcox_res.df$V1,method = "bonferroni")
colnames(med_v_low_wilcox_res.df) <- c("Wilcox test p-value","False Discovery Rate")
med_v_low_wilcox_res.df %>% write.csv(file=paste(OUT_DIR,"med_v_low_wilcox_res.csv",sep=""))

combined_wilcox_res.df <- data.frame("rxn_n1" = rownames(high_v_med_wilcox_res.df),
                                     "rxn_n2" = rownames(med_v_low_wilcox_res.df),
                                     "direction" = as.character(rxn_pca_direction.nls),
                                     "difference" = as.numeric(rxn_pca_difference.nls),
                                     "transcripts" = as.character(transcript_in_rxn.nls),
                                     "High_v_med_p" = high_v_med_wilcox_res.df$`Wilcox test p-value`,
                                     "Med_v_low_p" = med_v_low_wilcox_res.df$`Wilcox test p-value`)

# more conservative to select higher p-value than to combine them
#library(metap)
combined_w_fisher <- combined_wilcox_res.df %>%
  dplyr::rowwise() %>%
  #dplyr::mutate(combined_p = max(High_v_med_p,Med_v_low_p))
  dplyr::mutate(combined_p = as.numeric((metap::sumlog(c(High_v_med_p,Med_v_low_p)) %>% .[3])))

combined_w_fisher$fdr <- p.adjust(combined_w_fisher$combined_p,method = "fdr")

combined_w_fisher %>% write.csv(file=paste(OUT_DIR,"rxn_combined_w_fisher.csv",sep=""))

# compare gene transcript counts
vst.count.mtx.train$ENS_ID <- rownames(vst.count.mtx.train)

high_v_med_wilcox_res.nls <- list()
med_v_low_wilcox_res.nls <- list()
transcript_direction.nls <- list()
transcript_difference.nls <- list()
n_vst_train <- nrow(vst.count.mtx.train)
for(ens_idx in seq(1:n_vst_train)){
  w1 <- wilcox.test(x=as.numeric(vst.count.mtx.train[ens_idx,high_prolif_samples]),
                   y=as.numeric(vst.count.mtx.train[ens_idx,med_prolif_samples]))
  w2 <- wilcox.test(x=as.numeric(vst.count.mtx.train[ens_idx,med_prolif_samples]),
                   y=as.numeric(vst.count.mtx.train[ens_idx,low_prolif_samples]))
  transcript_direction <- "none"
  mean_hi <- mean(as.numeric(vst.count.mtx.train[ens_idx,high_prolif_samples]))
  mean_med <- mean(as.numeric(vst.count.mtx.train[ens_idx,med_prolif_samples]))
  mean_low <- mean(as.numeric(vst.count.mtx.train[ens_idx,low_prolif_samples]))
  if(mean_hi >= mean_med & mean_med >= mean_low){
    transcript_direction <- "positive"
  }else if(mean_low >= mean_med & mean_med >= mean_hi){
    transcript_direction <- "negative"
  }
  transcript_direction.nls[[vst.count.mtx.train$ENS_ID[ens_idx]]] <- transcript_direction
  transcript_difference.nls[[vst.count.mtx.train$ENS_ID[ens_idx]]] <- mean_hi - mean_low
  high_v_med_wilcox_res.nls[[vst.count.mtx.train$ENS_ID[ens_idx]]] <- w1$p.value
  med_v_low_wilcox_res.nls[[vst.count.mtx.train$ENS_ID[ens_idx]]] <- w2$p.value
  if(mod(ens_idx,50) == 0){
    print(paste("Processed ",ens_idx,
                " of ",n_vst_train,
                " transcripts (",round((ens_idx + 1)/n_vst_train,digits = 3) * 100,"%)...",
                sep=""))
    flush.console()
  }
}

saveRDS(transcript_direction.nls,file=paste(OUT_DIR,"transcript_direction_nls.Rds",sep=""))
saveRDS(transcript_difference.nls,file=paste(OUT_DIR,"transcript_difference_nls.Rds",sep=""))
saveRDS(high_v_med_wilcox_res.nls,file=paste(OUT_DIR,"ens_high_v_med_wilcox_res_nls.Rds",sep=""))
saveRDS(med_v_low_wilcox_res.nls,file=paste(OUT_DIR,"ens_med_v_low_wilcox_res_nls.Rds",sep=""))

high_v_med_wilcox_res.nls <- readRDS(file = paste(OUT_DIR,"ens_high_v_med_wilcox_res_nls.Rds",sep=""))
med_v_low_wilcox_res.nls <- readRDS(file = paste(OUT_DIR,"ens_med_v_low_wilcox_res_nls.Rds",sep=""))

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

saveRDS(high_v_med_wilcox_res.df,file=paste(OUT_DIR,"ens_high_v_med_wilcox_res_df.Rds",sep=""))
high_v_med_wilcox_res.df$fdr <- p.adjust(high_v_med_wilcox_res.df$V1,method = "bonferroni")
colnames(high_v_med_wilcox_res.df) <- c("Wilcox test p-value","False Discovery Rate")
high_v_med_wilcox_res.df %>% write.csv(file=paste(OUT_DIR,"ens_high_v_med_wilcox_res.csv",sep=""))

saveRDS(med_v_low_wilcox_res.df,file=paste(OUT_DIR,"ens_med_v_low_wilcox_res_df.Rds",sep=""))
med_v_low_wilcox_res.df$fdr <- p.adjust(med_v_low_wilcox_res.df$V1,method = "bonferroni")
colnames(med_v_low_wilcox_res.df) <- c("Wilcox test p-value","False Discovery Rate")
med_v_low_wilcox_res.df %>% write.csv(file=paste(OUT_DIR,"ens_med_v_low_wilcox_res.csv",sep=""))

combined_wilcox_res.df <- data.frame("ens_n1" = rownames(high_v_med_wilcox_res.df),
                                     "ens_n2" = rownames(med_v_low_wilcox_res.df),
                                     "direction" = as.character(transcript_direction.nls),
                                     "difference" = as.numeric(transcript_difference.nls),
                                     "High_v_med_p" = high_v_med_wilcox_res.df$`Wilcox test p-value`,
                                     "Med_v_low_p" = med_v_low_wilcox_res.df$`Wilcox test p-value`)

# more conservative to select higher p-value than to combine them
#library(metap)
combined_w_fisher <- combined_wilcox_res.df %>%
  dplyr::rowwise() %>%
  #dplyr::mutate(combined_p = max(High_v_med_p,Med_v_low_p))
  dplyr::mutate(combined_p = as.numeric((metap::sumlog(c(High_v_med_p,Med_v_low_p)) %>% .[3])))

combined_w_fisher$fdr <- p.adjust(combined_w_fisher$combined_p,method = "fdr")

combined_w_fisher %>% write.csv(file=paste(OUT_DIR,"ens_combined_w_fisher.csv",sep=""))

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

