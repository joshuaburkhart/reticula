# Review plot_ly function used in https://github.com/joshuaburkhart/reticula/blob/master/src/r/validation/srp035988/pca_and_knn_calculation.R 
library(pheatmap)
library(magrittr)
library(dplyr)
library(plotly)
library(ggplot2)

SRC_RXN <- "R-HSA-8956140"
HUB_RXN <- "R-HSA-8956184"
OUT_DIR <- "/home/jgburk/PycharmProjects/reticula/data/SRP035988/output/"
Relaxed <- c("ENSG00000108671","ENSG00000143106","ENSG00000092010","ENSG00000126067","ENSG00000100567","ENSG00000173692","ENSG00000163636","ENSG00000142507")

rxn2ensembls.nls <- readRDS(paste(OUT_DIR,"rxn2ensembls_nls.Rds",sep=""))
rxn_pca.df <- readRDS(paste(OUT_DIR,"rxn_pca_nls.Rds",sep="")) %>% as.data.frame()
tissue.vec <- readRDS(paste(OUT_DIR,"srp035988_tissue_vec.Rds",sep=""))
labelled_edge_weights.df <- read.table(paste(OUT_DIR,"labelled_edge_weights.csv",sep=""),
                                       stringsAsFactors = FALSE,sep = ",",header = TRUE)

# extract unique reactions from edges of hub reaction
hub_edges.df <- labelled_edge_weights.df %>%
  dplyr::filter(Preceeding_Reaction == HUB_RXN | Following_Reaction == HUB_RXN) %>%
  dplyr::select(Preceeding_Reaction, Following_Reaction)

hub_edge_names.vec <- c(hub_edges.df$Preceeding_Reaction,hub_edges.df$Following_Reaction) %>% unique() %>% make.names()

# 3d plot of those reactions pc1 values colored by phenotype
hub_edge_rxn_pca.df <- rxn_pca.df[,hub_edge_names.vec]

hub_edge_rxn_pca <-
  prcomp(hub_edge_rxn_pca.df,
         scale. = T)

pca.d <- data.frame(
  PC1 = hub_edge_rxn_pca$x[, 1],
  PC2 = hub_edge_rxn_pca$x[, 2],
  PC3 = hub_edge_rxn_pca$x[, 3],
  Section = tissue.vec
)

ggplot(pca.d) +
  geom_point(aes(x = PC1, y = PC2, colour = Section)) +
  theme_bw()

psoriasis_color = "#990000"
normal_color = "#BA957A"

plot_ly(
  x = pca.d$PC1,
  y = pca.d$PC2,
  z = pca.d$PC3,
  type = "scatter3d",
  mode = "markers",
  color = pca.d$Section,
  colors = c(psoriasis_color,normal_color),
  text = rownames(hub_edge_rxn_pca.df),
  size = 1
)

# investigate transcripts with significant wilcoxon p-values, update with significant positive/negative association directions
ens_ids <- rxn2ensembls.nls[[HUB_RXN]]
src_ens_ids <- rxn2ensembls.nls[[SRC_RXN]]

vst.count.mtx.train <- readRDS(paste(OUT_DIR,"vst_count_mtx_train.Rds",sep=""))

sig_expr.df <- data.frame(tissue_group = factor(tissue.vec))

for(ens_id in ens_ids){
  samples_by_count.df <- data.frame(
    expression_value = vst.count.mtx.train[ens_id, ] %>% t() %>% .[,1],
    tissue_group = factor(tissue.vec))
  
  sig_expr.df[[ens_id]] <- samples_by_count.df$expression_value
  
  wilcox_result = pairwise.wilcox.test(samples_by_count.df$expression_value,
                             samples_by_count.df$tissue_group,
                             p.adjust.method = "fdr")
  
  svg(filename = paste(OUT_DIR,ens_id,"_expr_v_grps.png",sep=""),width = 11, height = 10)
  plt <- ggplot(samples_by_count.df, aes(x=expression_value,
                                  y=tissue_group,
                                  color= tissue_group,
                                  fill = tissue_group)) +
    geom_violin() +
    coord_flip() +
    geom_boxplot(width  =0.15,aes(fill="x",color="y")) +
    labs(title=paste(ens_id," expression wilcox p-value = ",signif(wilcox_result$p.value[1],3),sep=""),
         x="Normalized Expression Value",
         y = "Tissue Group") +
    theme_bw() +
    theme(legend.position = "none") +
    scale_color_manual(values = c(psoriasis_color,normal_color, "black")) +
    scale_fill_manual(values = c(psoriasis_color, normal_color, "white"))
  print(plt)
  dev.off()
  kruskal.test(expression_value ~ tissue_group, data = samples_by_count.df)
}

pearson_cor <- cor(sig_expr.df[,Relaxed])
svg(filename =paste(OUT_DIR,"pheatmap.svg",sep=""),width = 11,height = 10)
pheatmap(pearson_cor, display_numbers = round(pearson_cor,3), fontsize_number = 14, number_color = "black")
dev.off()

library(tsne)
library(plotly)
library(viridis)

# from https://plotly.com/r/t-sne-and-umap-projections/
tsne_obj <- tsne(t(vst.count.mtx.train), initial_dims = 2) # takes ~8 minutes
tsne_df <- data.frame(tsne_obj)  
labelled_tsne_df <- cbind(tsne_df,as.factor(tissue.vec))

svg(filename = paste(OUT_DIR,"psoriasis_tsne.svg",sep=""),width = 11,height = 10)
plt <- ggplot2::ggplot(labelled_tsne_df, aes(x = X1, y = X2)) +
  ggplot2::geom_point(aes(color = as.factor(tissue.vec)),size=4) +
  scale_color_manual(values = c(psoriasis_color,normal_color)) +
  theme_minimal() +
  theme(legend.position = "none", title = element_text("Psoriasis TSNE"))
print(plt)
dev.off()

for(sig_ens in Relaxed){
  svg(filename = paste(OUT_DIR,sig_ens,"_expression_tsne.svg",sep=""),width = 11,height = 10)
  plt <- ggplot2::ggplot(labelled_tsne_df, aes(x = X1, y = X2)) +
    ggplot2::geom_point(aes(color = sig_expr.df[,sig_ens]),size=4) +
    scale_color_viridis(option = "inferno", direction = -1) +
    theme_minimal() +
    theme(legend.position = "none", title = element_text(sig_ens))
  print(plt)
  dev.off()
}
