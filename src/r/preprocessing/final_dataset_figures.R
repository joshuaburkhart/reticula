library(NMF)
library(DESeq2)
library(fossil)
library(ggplot2)
library(viridis)
library(magrittr)
library(pheatmap)
library(RColorBrewer)
library(SummarizedExperiment)

OUT_DIR <- "/home/burkhart/Software/reticula/data/aim1/output/"

#DESeq2 PCA plot
vst.counts <- readRDS(paste(OUT_DIR,"vst_counts.Rds",sep=""))

DESeq2::plotPCA(vst.counts,intgroup="Tissue")

vst.count.mtx <- vst.counts %>% SummarizedExperiment::assay() %>% as.data.frame()
gtex_tissue_detail.vec <- readRDS(paste(OUT_DIR,"gtex_tissue_detail_vec.Rds",sep=""))
rxn2ensembls.nls <- readRDS(paste(OUT_DIR,"rxn2ensembls_nls.Rds",sep=""))
rxns <- rxn2ensembls.nls %>% names()

#calculate clustering coefficents for each reaction
for(rxn_id in rxns){
 ensembl_ids <- rxn2ensembls.nls[[rxn_id]]
 # combined.df columns ordered by clustering
 expr_matrix <- vst.count.mtx[ensembl_ids,]
 km_obj <- kmeans(t(expr_matrix),centers = 54)
 ari <- fossil::adj.rand.index(as.factor(km_obj$cluster),
                               as.factor(gtex_tissue_detail.vec))
 purity <- NMF::purity(as.factor(km_obj$cluster),
                       as.factor(gtex_tissue_detail.vec))
}

#"ENSG00000010671" "ENSG00000146535"
p <- ggplot(zt,aes(BTK,GNA12))
p <- p + geom_point(aes(color = gtex_tissue_detail.vec))
p <- p + theme(legend.position = "none")
ggsave("~/test123.png")