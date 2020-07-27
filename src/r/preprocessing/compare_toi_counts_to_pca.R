set.seed(88888888) # maximum luck

library(magrittr)
library(ggplot2)
library(plotly)

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
rownames(df) <- names(rxn_pca.nls)
colnames(df) <- names(rxn_pca.nls[[1]])
  
pca_pca.obj <- prcomp(t(pca.df %>% dplyr::distinct()),scale. = T)
pca_pca.df <- data.frame(PC1 = pca_pca.obj$x[,1],
                    PC2 = pca_pca.obj$x[,2],
                    PC3 = pca_pca.obj$x[,3],
                    Section = gtex_tissue_detail.vec.train)

vst.count.mtx.train_pca.obj <- prcomp(t(vst.count.mtx.train %>% dplyr::distinct()),scale. = T)
vst.count.mtx.train_pca.df <- data.frame(PC1 = vst.count.mtx.train_pca.obj$x[,1],
                                         PC2 = vst.count.mtx.train_pca.obj$x[,2],
                                         PC3 = vst.count.mtx.train_pca.obj$x[,3],
                                         Section = gtex_tissue_detail.vec.train)

ggplot2::ggplot(pca_pca.df) +
  geom_point(aes(x=PC1,y=PC2,colour=Section)) +
  theme_bw()

ggplot2::ggplot(vst.count.mtx.train_pca.df) +
  geom_point(aes(x=PC1,y=PC2,colour=Section)) +
  theme_bw()

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


end_time <- Sys.time()
print(paste("Start: ",start_time," End: ",end_time," Difference: ",end_time - start_time))