set.seed(88888888) # maximum luck

library(magrittr)
library(ggplot2)
library(ggiraph)
library(plotly)
library(factoextra)

start_time <- Sys.time()

OUT_DIR <- "/Users/burkhajo/Software/reticula/data/aim1/output/"

gtex_tissue_detail.vec <- readRDS(paste(OUT_DIR,"gtex_tissue_detail_vec.Rds",sep=""))
rxn_pca.nls <- readRDS(paste(OUT_DIR,"rxn_pca_nls.Rds",sep=""))
vst.count.mtx.train <- readRDS(paste(OUT_DIR,"vst_count_mtx_train.Rds",sep=""))

# take a look at toi samples...
toi_indices <- which(gtex_tissue_detail.vec == "Colon - Transverse" |
                       gtex_tissue_detail.vec == "Colon - Sigmoid")

# filter annotations
gtex_tissue_detail_vec_tis_of_interest <- gtex_tissue_detail.vec[toi_indices]

training_indices <- caret::createDataPartition(gtex_tissue_detail_vec_tis_of_interest,
                                               times = 1,
                                               p = 0.9,
                                               list = FALSE)

gtex_tissue_detail.vec.train <- gtex_tissue_detail_vec_tis_of_interest[training_indices]

# create pca and df to compare this training data with vst.count.mtx.train df and plot in two and three dimensions
pca.df <- data.frame(matrix(unlist(rxn_pca.nls), nrow=length(rxn_pca.nls), byrow=T))
rownames(pca.df) <- names(rxn_pca.nls)
colnames(pca.df) <- names(rxn_pca.nls[[1]])
  
pca_pca.obj <- prcomp(t(pca.df),scale. = T)
pca_pca.df <- data.frame(PC1 = pca_pca.obj$x[,1],
                    PC2 = pca_pca.obj$x[,2],
                    PC3 = pca_pca.obj$x[,3],
                    Section = gtex_tissue_detail.vec.train)

vst.count.mtx.train_pca.obj <- prcomp(t(vst.count.mtx.train),scale. = T)
vst.count.mtx.train_pca.df <- data.frame(PC1 = vst.count.mtx.train_pca.obj$x[,1],
                                         PC2 = vst.count.mtx.train_pca.obj$x[,2],
                                         PC3 = vst.count.mtx.train_pca.obj$x[,3],
                                         Section = gtex_tissue_detail.vec.train)

# pca plotting

ggplot2::ggplot(pca_pca.df) +
  geom_point(aes(x=PC1,y=PC2,colour=Section)) +
  theme_bw() +
  ggtitle("Reaction PC1")

ggplot2::ggplot(vst.count.mtx.train_pca.df) +
  geom_point(aes(x=PC1,y=PC2,colour=Section)) +
  theme_bw() +
  ggtitle("Transcript Count")

plotly::plot_ly(x=pca_pca.df$PC1,
        y=pca_pca.df$PC2,
        z=pca_pca.df$PC3,
        type="scatter3d",
        mode="markers",
        color=pca_pca.df$Section, 
        size = 1)

plotly::plot_ly(x=vst.count.mtx.train_pca.df$PC1,
        y=vst.count.mtx.train_pca.df$PC2,
        z=vst.count.mtx.train_pca.df$PC3,
        type="scatter3d", 
        mode="markers",
        color=vst.count.mtx.train_pca.df$Section, 
        size = 1)

# pca analysis (from http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/#access-to-the-pca-results)

pca.var <- factoextra::get_pca(pca_pca.obj, "var") # Results for variable categories
count.var <- factoextra::get_pca(vst.count.mtx.train_pca.obj,"var")

pca.var.contributions <- pca.var$contrib
count.var.contributions <- count.var$contrib

pca.var.sorted_contribution <- pca.var.contributions[order(-pca.var.contributions[,"Dim.1"]),]
count.var.sorted_contribution <- count.var.contributions[order(-count.var.contributions[,"Dim.1"]),]

# print top reactions & transcripts
pca.var.sorted_contribution[1:10,1]
count.var.sorted_contribution[1:10,1]

# R-HSA-2316434 (PI3K phosphorylates PIP2 to PIP3) emerges as a top 10 separating reaction

z <- readRDS(paste(OUT_DIR,"rxn2ensembls_nls.Rds",sep=""))
# show transcripts annotated to this reaction
z[["R-HSA-2316434"]]

# show transcripts in this reaction not highly ranked by transcript pca
match(z[["R-HSA-2316434"]],rownames(count.var.sorted_contribution)) %>% sort()

factoextra::fviz_eig(pca_pca.obj)
factoextra::fviz_eig(vst.count.mtx.train_pca.obj)

# top 10 phosphorylation reaction
b <- prcomp(t(vst.count.mtx.train[z[["R-HSA-2316434"]],]),scale.=T)
bd <- data.frame(pc1 = b$x[,1],pc2 = b$x[,2],Section = gtex_tissue_detail.vec.train)
ggplot2::ggplot(bd) + geom_point(aes(x=pc1,y=pc2,colour = Section)) + theme_bw() + ggtitle("phosphorylation rxn")

# 10000th reaction
b <- prcomp(t(vst.count.mtx.train[z[["R-HSA-8939335"]],]),scale.=T)
bd <- data.frame(pc1 = b$x[,1],pc2 = b$x[,2],Section = gtex_tissue_detail.vec.train)
ggplot2::ggplot(bd) + geom_point(aes(x=pc1,y=pc2,colour = Section)) + theme_bw() + ggtitle("10k rxn")

# compare pca results to knn results
knn_res <- readRDS(paste(OUT_DIR,"toi_summary_df.Rds",sep=""))
pca_v_knn <- data.frame(RXN_ID = knn_res$RXN_ID,
                        MISCLASS = knn_res$MISCLASS,
                        ARI = knn_res$ARI,
                        ECOUNT = knn_res$ECOUNT,
                        LOADING = pca.var.contributions[,"Dim.1"])

pca_v_knn <- pca_v_knn %>% dplyr::arrange(MISCLASS)

plot.obj <- ggplot2::ggplot(pca_v_knn) + 
  ggiraph::geom_point_interactive(aes(x=ARI,
                                      y=LOADING,
                                      colour = MISCLASS,
                                      tooltip=RXN_ID,
                                      data_id = RXN_ID)) +
  theme_bw() + 
  ggtitle("ARI v LOADING")

girafe(ggobj = plot.obj)

b <- prcomp(t(vst.count.mtx.train[z[["R-HSA-983147"]],]),scale.=T)
bd <- data.frame(pc1 = b$x[,1],pc2 = b$x[,2],Section = gtex_tissue_detail.vec.train)
ggplot2::ggplot(bd) + geom_point(aes(x=pc1,y=pc2,colour = Section)) + theme_bw() + ggtitle("R-HSA-983147")

b <- prcomp(t(vst.count.mtx.train[z[["R-HSA-6809663"]],]),scale.=T)
bd <- data.frame(pc1 = b$x[,1],pc2 = b$x[,2],Section = gtex_tissue_detail.vec.train)
ggplot2::ggplot(bd) + geom_point(aes(x=pc1,y=pc2,colour = Section)) + theme_bw() + ggtitle("R-HSA-6809663")

b <- prcomp(t(vst.count.mtx.train[z[["R-HSA-6787642"]],]),scale.=T)
bd <- data.frame(pc1 = b$x[,1],pc2 = b$x[,2],Section = gtex_tissue_detail.vec.train)
ggplot2::ggplot(bd) + geom_point(aes(x=pc1,y=pc2,colour = Section)) + theme_bw() + ggtitle("R-HSA-6787642")

end_time <- Sys.time()
print(paste("Start: ",start_time," End: ",end_time," Difference: ",end_time - start_time))