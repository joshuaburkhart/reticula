library(magrittr)

OUT_DIR <- "/home/burkhart/Software/reticula/data/aim1/output/"
IN_DIR <- "/home/burkhart/Software/reticula/data/aim2/input/"

X <- readRDS(paste(OUT_DIR, "rxn_pca_nls.Rds", sep = ""))
Y <- readRDS(paste(OUT_DIR,"gtex_tissue_detail_vec_train.Rds",sep=""))

X <- as.data.frame(X)

Y <- as.data.frame(Y)

write.table(X,
          file=paste(IN_DIR,"resnet_node_features.txt",sep=""),
          row.names = FALSE,
          col.names = FALSE)
write.table(Y,
          file=paste(IN_DIR,"resnet_graph_targets.txt",sep=""),
          row.names = FALSE,
          col.names = FALSE)
