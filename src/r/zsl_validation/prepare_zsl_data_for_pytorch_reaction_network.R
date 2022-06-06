library(magrittr)

TCGA_OUT_DIR <- "/home/jgburk/PycharmProjects/reticula/data/tcga/output/"
TCGA_IN_DIR <- "/home/jgburk/PycharmProjects/reticula/data/tcga/input/"

X <- readRDS(paste(TCGA_OUT_DIR, "zsl_tcga_rxn_pca_nls.Rds", sep = ""))
Y <- readRDS(paste(TCGA_OUT_DIR,"zsl_tcga_tissue_vec_train.Rds",sep=""))

tissue_alphabetical_order <- order(Y)

X <- as.data.frame(X) %>% .[tissue_alphabetical_order,]

Y <- as.data.frame(Y[tissue_alphabetical_order])

write.table(X,
          file=paste(TCGA_IN_DIR,"zsl_tcga_node_features.txt",sep=""),
          row.names = FALSE,
          col.names = FALSE)
write.table(Y,
          file=paste(TCGA_IN_DIR,"zsl_tcga_graph_targets.txt",sep=""),
          row.names = FALSE,
          col.names = FALSE)

call_check.df <- data.frame(Sample_ID = rownames(X),
                            Tissue = Y$`Y[tissue_alphabetical_order]`,
                            Dataset = rep(c("Tuning Set","Test Set"),length(X[,1])/2))
write.table(call_check.df,
            file=paste(TCGA_OUT_DIR,"zsl_tcga_call_check.csv",sep=""),
            row.names = FALSE,
            col.names = TRUE,
            sep = ",")

GTEX_OUT_DIR <- "/home/jgburk/PycharmProjects/reticula/data/gtex/output/"
GTEX_IN_DIR <- "/home/jgburk/PycharmProjects/reticula/data/gtex/input/"

X <- readRDS(paste(GTEX_OUT_DIR, "zsl_gtex_rxn_pca_nls.Rds", sep = ""))
Y <- readRDS(paste(GTEX_OUT_DIR,"zsl_gtex_tissue_vec_train.Rds",sep=""))
E <- read.table(paste(GTEX_IN_DIR,"ReactionNetwork_Rel_71_122820.txt",sep=""))

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
            file=paste(GTEX_IN_DIR,"edgeLabels.csv",sep=""),
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

assertthat::are_equal(length(z),length(y))

write.table(z,
            file=paste(GTEX_IN_DIR,"rxn2nodeLabel_nls.csv",sep=""),
            row.names = TRUE,
            col.names = FALSE)
write.table(y,
            file=paste(GTEX_IN_DIR,"nodeLabel2rxn_nls.csv",sep=""),
            row.names = TRUE,
            col.names = FALSE)

E <- data.frame(node1 = node1,
                node2 = node2)

tissue_alphabetical_order <- order(Y)

X <- as.data.frame(X) %>% .[tissue_alphabetical_order,]

Y <- as.data.frame(Y[tissue_alphabetical_order])

write.table(E,
            file=paste(GTEX_IN_DIR,"edges.txt",sep=""),
            row.names = FALSE,
            col.names = FALSE)
write.table(X,
            file=paste(GTEX_IN_DIR,"zsl_gtex_node_features.txt",sep=""),
            row.names = FALSE,
            col.names = FALSE)
write.table(Y,
            file=paste(GTEX_IN_DIR,"zsl_gtex_graph_targets.txt",sep=""),
            row.names = FALSE,
            col.names = FALSE)
