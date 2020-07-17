library(DESeq2)
library(ggplot2)
library(magrittr)
library(ggfortify)
library(SummarizedExperiment)

start_time <- Sys.time()

IN_DIR <- "/home/burkhart/Software/reticula/data/aim1/input/"
OUT_DIR <- "/home/burkhart/Software/reticula/data/aim1/output/"

gtex_tissue_detail.vec <- readRDS(paste(OUT_DIR,"gtex_tissue_detail_vec.Rds",sep=""))

rxn_ari.nls <- readRDS(paste(OUT_DIR,"rxn_ari_nls.Rds",sep=""))
#rxn_gini.nls <- readRDS(paste(OUT_DIR,"rxn_gini_nls.Rds",sep=""))
rxn_ensembl_counts.nls <- readRDS(paste(OUT_DIR,"rxn_ensembl_counts.nls",sep=""))
rxn_names <- rxn_ari.nls %>% names()

rxn_ari.df <- as.data.frame(t(dplyr::bind_rows(rxn_ari.nls)))
rxn_ensembl_counts.df <- as.data.frame(t(dplyr::bind_rows(rxn_ensembl_counts.nls)))

boxplot(as.list(rxn_ari.df),main="ARI distriubtion across all reactions by k")
colSums(rxn_ari.df) %>% plot(.,main="ARI sum across all reactions by k")

rxn_ari.df$ID <- rownames(rxn_ari.df)
rxn_ari.df <- rxn_ari.df %>% dplyr::arrange(dplyr::desc(V11))

rxn_cluster_stats.df <- data.frame(ID = rxn_names,
                                   ARI = unlist(rxn_ari.nls),
                                   #GINI = unlist(rxn_gini.nls),
                                   ECOUNT = unlist(rxn_ensembl_counts.nls))

# plot(x = rxn_cluster_stats.df$GINI,
#      y = rxn_cluster_stats.df$ARI,
#      pch = 16,
#      main = "Reaction clustering Gini vs ARI",
#      xlab = "Gini Index",
#      ylab = "Adjusted Rand Index (ARI)")

plot(x = rxn_cluster_stats.df$ECOUNT,
     y = rxn_cluster_stats.df$ARI,
     pch = 4,
     main = "Reaction clustering # transcripts vs ARI",
     xlab = "Reaction Transcript Count",
     ylab = "Adjusted Rand Index")

# plot(x = rxn_cluster_stats.df$ECOUNT,
#      y = rxn_cluster_stats.df$GINI,
#      pch = 16,
#      main = "Reaction clustering # transcripts vs Gini",
#      xlab = "Reaction Transcript Count",
#      ylab = "Gini Index")

rxn_ari.nls.sorted <- rxn_ari.nls[order(unlist(rxn_ari.nls),decreasing=TRUE)]

vst.counts.assay <- readRDS(paste(OUT_DIR,"vst_counts.Rds",sep="")) %>% SummarizedExperiment::assay()

ari_rxn1_pca <- prcomp(t(vst.counts.assay[rxn2ensembls.nls[[rxn_ari.nls.sorted %>% names() %>% .[1]]],]),scale.=T)

saveRDS(ari_rxn1_pca,paste(OUT_DIR,"top_ari_rxn_pca.Rds",sep=""))

plot(x = ari_rxn1_pca$x[,1],
     y = ari_rxn1_pca$x[,2],
     pch=16,
     main = paste(rxn_ari.nls.sorted %>% names() %>% .[1],": 1st ARI",sep=""),
     xlab = "PC 1",
     ylab = "PC 2",
     col = as.numeric(as.factor(gtex_tissue_detail.vec)))

ari_rxn2_pca <- prcomp(t(vst.counts.assay[rxn2ensembls.nls[[rxn_ari.nls.sorted %>% names() %>% .[2]]],]),scale.=T)

plot(x = ari_rxn2_pca$x[,1],
     y = ari_rxn2_pca$x[,2],
     pch=16,
     main = paste(rxn_ari.nls.sorted %>% names() %>% .[2],": 2nd ARI",sep=""),
     xlab = "PC 1",
     ylab = "PC 2",
     col = as.numeric(as.factor(gtex_tissue_detail.vec)))

ari_rxn3_pca <- prcomp(t(vst.counts.assay[rxn2ensembls.nls[[rxn_ari.nls.sorted %>% names() %>% .[3]]],]),scale.=T)

plot(x = ari_rxn3_pca$x[,1],
     y = ari_rxn3_pca$x[,2],
     pch=16,
     main = paste(rxn_ari.nls.sorted %>% names() %>% .[3],": 3rd ARI",sep=""),
     xlab = "PC 1",
     ylab = "PC 2",
     col = as.numeric(as.factor(gtex_tissue_detail.vec)))

ari_rxn5k_pca <- prcomp(t(vst.counts.assay[rxn2ensembls.nls[[rxn_ari.nls.sorted %>% names() %>% .[5000]]],]),scale.=T)

plot(x = ari_rxn5k_pca$x[,1],
     y = ari_rxn5k_pca$x[,2],
     pch=16,
     main = paste(rxn_ari.nls.sorted %>% names() %>% .[5000],": 5000th ARI",sep=""),
     xlab = "PC 1",
     ylab = "PC 2",
     col = as.numeric(as.factor(gtex_tissue_detail.vec)))

end_time <- Sys.time()
print(paste("Start: ",start_time," End: ",end_time," Difference: ",end_time - start_time))


