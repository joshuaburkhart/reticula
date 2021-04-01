set.seed(88888888) # maximum luck

library(magrittr)
library(ggplot2)
library(ggiraph)
library(plotly)
library(plyr)
library(reshape2)
library(factoextra)

start_time <- Sys.time()

OUT_DIR <- "/home/burkhart/Software/reticula/data/aim1/output/"

gtex_tissue_detail.vec <- readRDS(paste(OUT_DIR,"gtex_tissue_detail_vec.Rds",sep=""))
rxn_pca.nls <- readRDS(paste(OUT_DIR,"rxn_pca_nls.Rds",sep=""))
vst.count.mtx.train <- readRDS(paste(OUT_DIR,"vst_count_mtx_train.Rds",sep=""))

# returns sixth prcomp object
# z <- pca_results.nls[[names(pca_results.nls)[6]]]

# returns proportion of variance for each PC
# summary(z) %>% .$importance %>% .[2,]
# PC1     PC2 
# 0.83185 0.16815 

# > z
# [1] 1 2
# > p
# [1] 1 2 3
# > f
# [1] 3
# > q <- plyr::rbind.fill(as.data.frame(t(z)),as.data.frame(t(p)),as.data.frame(t(f)))
# > q
# V1 V2 V3
# 1  1  2 NA
# 2  1  2  3
# 3  3 NA NA

pca_data_fns <- c("full_rxn_pca_results_nls0-1000.Rds",
                  "full_rxn_pca_results_nls1000-2000.Rds",
"full_rxn_pca_results_nls2000-3000.Rds",
"full_rxn_pca_results_nls3000-4000.Rds",
"full_rxn_pca_results_nls4000-5000.Rds",
"full_rxn_pca_results_nls5000-6000.Rds",
"full_rxn_pca_results_nls6000-7000.Rds",
"full_rxn_pca_results_nls7000-8000.Rds",
"full_rxn_pca_results_nls8000-9000.Rds",
"full_rxn_pca_results_nls9000-10000.Rds",
"full_rxn_pca_results_nls.Rds")

pca_importance.df <- data.frame()

for(pca_data_fn in pca_data_fns){
  pca_results.nls <- readRDS(file=paste(OUT_DIR,pca_data_fn,sep=""))
  for(rxn_name in names(pca_results.nls)){
    z <- pca_results.nls[[rxn_name]]
    pc_imp_vec <- summary(z) %>% .$importance %>% .[2,]
    pca_importance.df <- plyr::rbind.fill(pca_importance.df,
                                          as.data.frame(t(pc_imp_vec)))
  }
}
pca_importance.df %>% .[,1:10] %>% boxplot()

# take a look at toi samples...
toi_indices <- seq(1,length(gtex_tissue_detail.vec))
#which(
#   gtex_tissue_detail.vec == "Colon - Transverse" |
#      gtex_tissue_detail.vec == "Colon - Sigmoid"
#)

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
  ggtitle("Colon Reactions' 1st PC's")

ggplot2::ggplot(vst.count.mtx.train_pca.df) +
  geom_point(aes(x=PC1,y=PC2,colour=Section)) +
  theme_bw() +
  ggtitle("Colon Transcript Counts")

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
z[["R-HSA-983156"]]

# show transcripts in this reaction not highly ranked by transcript pca
match(z[["R-HSA-983156"]],rownames(count.var.sorted_contribution)) %>% sort()

factoextra::fviz_eig(pca_pca.obj)
factoextra::fviz_eig(vst.count.mtx.train_pca.obj)

# compare pca results to knn results
knn_res <- readRDS(paste(OUT_DIR,"toi_summary_df.Rds",sep=""))
pca_v_knn <- data.frame(RXN_ID = knn_res$RXN_ID,
                        MISCLASS = knn_res$MISCLASS,
                        zMISCLASS = (knn_res$MISCLASS - mean(knn_res$MISCLASS))/sd(knn_res$MISCLASS),
                        ARI = knn_res$ARI,
                        zARI = (knn_res$ARI - mean(knn_res$ARI))/sd(knn_res$ARI),
                        ECOUNT = knn_res$ECOUNT,
                        LOADING = pca.var.contributions[,"Dim.1"],
                        zLOADING = (pca.var.contributions[,"Dim.1"] - mean(pca.var.contributions[,"Dim.1"]))/sd(pca.var.contributions[,"Dim.1"]))

pca_v_knn <- pca_v_knn %>% dplyr::arrange(zMISCLASS) %>% dplyr::filter(ECOUNT >=2)

plot.obj <- ggplot2::ggplot(pca_v_knn) + 
  ggiraph::geom_point_interactive(aes(x=zLOADING,
                                      y=zARI,
                                      colour = zMISCLASS,
                                      tooltip=RXN_ID,
                                      data_id = RXN_ID)) +
  theme_bw() + 
  ggtitle("Separability v Spread (rxns >= 2 transcripts)")

girafe(ggobj = plot.obj)

# quadrant 1
b <- prcomp(t(vst.count.mtx.train[z[["R-HSA-983147"]],]),scale.=T)
bd <- data.frame(pc1 = b$x[,1],pc2 = b$x[,2],Section = gtex_tissue_detail.vec.train)
ggplot2::ggplot(bd) + geom_point(aes(x=pc1,y=pc2,colour = Section)) + theme_bw() + ggtitle("R-HSA-983147: hi spread, hi separability")

# quadrant 2
b <- prcomp(t(vst.count.mtx.train[z[["R-HSA-6809663"]],]),scale.=T)
bd <- data.frame(pc1 = b$x[,1],pc2 = b$x[,2],Section = gtex_tissue_detail.vec.train)
ggplot2::ggplot(bd) + geom_point(aes(x=pc1,y=pc2,colour = Section)) + theme_bw() + ggtitle("R-HSA-6809663: lo spread, hi separability")

# quadrant 4
b <- prcomp(t(vst.count.mtx.train[z[["R-HSA-8937844"]],]),scale.=T)
bd <- data.frame(pc1 = b$x[,1],pc2 = b$x[,2],Section = gtex_tissue_detail.vec.train)
ggplot2::ggplot(bd) + geom_point(aes(x=pc1,y=pc2,colour = Section)) + theme_bw() + ggtitle("R-HSA-8937844: hi spread, lo separability")

# quadrant 3
b <- prcomp(t(vst.count.mtx.train[z[["R-HSA-378978"]],]),scale.=T)
bd <- data.frame(pc1 = b$x[,1],pc2 = b$x[,2],Section = gtex_tissue_detail.vec.train)
ggplot2::ggplot(bd) + geom_point(aes(x=pc1,y=pc2,colour = Section)) + theme_bw() + ggtitle("R-HSA-378978: lo spread, lo separability")

# review boxplots
# quadrant 1
b <- as.data.frame(t(vst.count.mtx.train[z[["R-HSA-983147"]],]))
b$Section <- unlist(gtex_tissue_detail.vec.train)
dat.m <- melt(b,id.vars='Section', measure.vars=colnames(b[,-ncol(b)]))
ggplot(dat.m) + geom_boxplot(aes(x=variable, y=value, fill=Section)) +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("R-HSA-983147 Transcript Distributions")

# quadrant 2
b <- as.data.frame(t(vst.count.mtx.train[z[["R-HSA-6809663"]],]))
b$Section <- unlist(gtex_tissue_detail.vec.train)
dat.m <- melt(b,id.vars='Section', measure.vars=colnames(b[,-ncol(b)]))
ggplot(dat.m) + geom_boxplot(aes(x=variable, y=value, fill=Section)) +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("R-HSA-6809663 Transcript Distributions")

# quadrant 4
b <- as.data.frame(t(vst.count.mtx.train[z[["R-HSA-8937844"]],]))
b$Section <- unlist(gtex_tissue_detail.vec.train)
dat.m <- melt(b,id.vars='Section', measure.vars=colnames(b[,-ncol(b)]))
ggplot(dat.m) + geom_boxplot(aes(x=variable, y=value, fill=Section)) +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("R-HSA-8937844 Transcript Distributions")

# quadrant 3
b <- as.data.frame(t(vst.count.mtx.train[z[["R-HSA-378978"]],]))
b$Section <- unlist(gtex_tissue_detail.vec.train)
dat.m <- melt(b,id.vars='Section', measure.vars=colnames(b[,-ncol(b)]))
ggplot(dat.m) + geom_boxplot(aes(x=variable, y=value, fill=Section)) +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("R-HSA-378978 Transcript Distributions")


end_time <- Sys.time()
print(paste("Start: ",start_time," End: ",end_time," Difference: ",end_time - start_time))