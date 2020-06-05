library(DESeq2)
library(ggplot2)
library(magrittr)
library(ggfortify)
library(SummarizedExperiment)

start_time <- Sys.time()

IN_DIR <- "/home/burkhart/Software/reticula/data/aim1/input/"
OUT_DIR <- "/home/burkhart/Software/reticula/data/aim1/output/"

rxn_ari.nls <- readRDS(paste(OUT_DIR,"rxn_ari_nls.Rds",sep=""))
rxn_ari.nls.sorted <- rxn_ari.nls[order(unlist(rxn_ari.nls),decreasing=TRUE)]

top_ari_rxn <- rxn_ari.nls.sorted %>% names() %>% .[1]
top_ari_rxn_ensembls <- rxn2ensembls.nls[[top_ari_rxn]]

vst.counts.assay <- readRDS(paste(OUT_DIR,"vst_counts.Rds",sep="")) %>% SummarizedExperiment::assay()

top_ari_rxn_pca <- prcomp(t(vst.counts.assay[top_ari_rxn_ensembls,]),scale.=T)

saveRDS(top_ari_rxn_pca,paste(OUT_DIR,"top_ari_rxn_pca.Rds",sep=""))

plot(top_ari_rxn_pca)

gtex_tissue_detail.vec <- readRDS(paste(OUT_DIR,"gtex_tissue_detail_vec.Rds",sep=""))

plot(top_ari_rxn_pca$x[,1],
     top_ari_rxn_pca$x[,2],
     pch=16,
     col = as.numeric(as.factor(gtex_tissue_detail.vec)))

end_time <- Sys.time()
print(paste("Start: ",start_time," End: ",end_time," Difference: ",end_time - start_time))