library(magrittr)

OUT_DIR <- "/home/burkhart/Software/reticula/data/aim1/output/"
IN_DIR <- "/home/burkhart/Software/reticula/data/aim2/input/"

X <- readRDS(paste(OUT_DIR, "rxn_pca_nls.Rds", sep = ""))
Y <- readRDS(paste(OUT_DIR,"gtex_tissue_detail_vec_train.Rds",sep=""))
E <- read.table(paste(IN_DIR,"ReactionNetwork_Rel_71_122820.txt",sep=""))

rxn2nodeLabel.nls <- list()
nodeLabel2rxn.nls <- list()
for(i in 1:length(X)){
  rxn2nodeLabel.nls[[names(X)[i]]] <- i
  nodeLabel2rxn.nls[[i]] <- names(X)[i]
}

E <- E %>%
  dplyr::filter(V1 %in% names(rxn2nodeLabel.nls)) %>%
  dplyr::filter(V3 %in% names(rxn2nodeLabel.nls)) %>%
  dplyr::select(V1,V3)

write.table(E,
            file=paste(IN_DIR,"edgeLabels.csv",sep=""),
            row.names = FALSE,
            col.names = FALSE)

node1 <- numeric()
node2 <- numeric()
for(i in 1:nrow(E)){
  node1 <- c(node1,rxn2nodeLabel.nls[[as.character(E$V1[i])]])
  node2 <- c(node2,rxn2nodeLabel.nls[[as.character(E$V3[i])]])
}

z <- unlist(rxn2nodeLabel.nls)
y <- unlist(nodeLabel2rxn.nls)

write.table(z,
            file=paste(IN_DIR,"rxn2nodeLabel_nls.csv",sep=""),
            row.names = TRUE,
            col.names = FALSE)
write.table(y,
            file=paste(IN_DIR,"nodeLabel2rxn_nls.csv",sep=""),
            row.names = TRUE,
            col.names = FALSE)

E <- data.frame(node1 = node1,
                node2 = node2)

X <- as.data.frame(X)

Y <- as.data.frame(Y)

write.table(E,
          file=paste(IN_DIR,"edges.txt",sep=""),
          row.names = FALSE,
          col.names = FALSE)
write.table(X,
          file=paste(IN_DIR,"node_features.txt",sep=""),
          row.names = FALSE,
          col.names = FALSE)
write.table(Y,
          file=paste(IN_DIR,"graph_targets.txt",sep=""),
          row.names = FALSE,
          col.names = FALSE)
