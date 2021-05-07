library(magrittr)

OUT_DIR <- "/home/burkhart/Software/reticula/data/aim1/output/"
IN_DIR <- "/home/burkhart/Software/reticula/data/aim2/input/"

X <- readRDS(paste(OUT_DIR, "rxn_pca_nls.Rds", sep = ""))
Y <- readRDS(paste(OUT_DIR,"gtex_tissue_detail_vec_train.Rds",sep=""))
E <- read.table(paste(IN_DIR,"pathway_reaction_id_edges.txt",sep=""),sep = "\t")

# ensure no reactions have any children
assertthat::are_equal(E$V2 %in% names(X) %>% as.numeric() %>% -1 %>% abs() %>% sum(),nrow(E))

pathway_node_ids <- E$V2 %>% unique() %>% as.character()

rxnAndPthwy2nodeLabel.nls <- list()
nodeLabel2rxnAndPthwy.nls <- list()
n_rxns <- length(X)
for(i in 1:n_rxns){
  rxnAndPthwy2nodeLabel.nls[[names(X)[i]]] <- i
  nodeLabel2rxnAndPthwy.nls[[i]] <- names(X)[i]
}
n_pthwys <- length(pathway_node_ids)
for(i in 1:n_pthwys){
  rxnAndPthwy2nodeLabel.nls[[pathway_node_ids[i]]] <- i + n_rxns
  nodeLabel2rxnAndPthwy.nls[[i + n_rxns]] <- pathway_node_ids[i]
  X[[pathway_node_ids[i]]] <- 1 #set all pathway weights to 1
}

E <- E %>%
  dplyr::filter(V1 %in% names(rxnAndPthwy2nodeLabel.nls)) %>%
  dplyr::filter(V2 %in% names(rxnAndPthwy2nodeLabel.nls)) %>%
  dplyr::select(V1,V2)

write.table(E,
            file=paste(IN_DIR,"pathway_hierarchy_edgeLabels.csv",sep=""),
            row.names = FALSE,
            col.names = FALSE)

node1 <- numeric()
node2 <- numeric()
for(i in 1:nrow(E)){
  node1 <- c(node1,rxnAndPthwy2nodeLabel.nls[[as.character(E$V1[i])]])
  node2 <- c(node2,rxnAndPthwy2nodeLabel.nls[[as.character(E$V2[i])]])
}

z <- unlist(rxnAndPthwy2nodeLabel.nls)
y <- unlist(nodeLabel2rxnAndPthwy.nls)

assertthat::are_equal(length(z),length(y))

write.table(z,
            file=paste(IN_DIR,"pathway_hierarchy_rxn2nodeLabel_nls.csv",sep=""),
            row.names = TRUE,
            col.names = FALSE)
write.table(y,
            file=paste(IN_DIR,"pathway_hierarchy_nodeLabel2rxn_nls.csv",sep=""),
            row.names = TRUE,
            col.names = FALSE)

E <- data.frame(node1 = node1,
                node2 = node2)

X <- as.data.frame(X)

Y <- as.data.frame(Y)

write.table(E,
          file=paste(IN_DIR,"pathway_hierarchy_edges.txt",sep=""),
          row.names = FALSE,
          col.names = FALSE)
write.table(X,
          file=paste(IN_DIR,"pathway_hierarchy_node_features.txt",sep=""),
          row.names = FALSE,
          col.names = FALSE)
write.table(Y,
          file=paste(IN_DIR,"pathway_hierarchy_graph_targets.txt",sep=""),
          row.names = FALSE,
          col.names = FALSE)

# ensure node edge indices do not exceed input bounds
assertthat::assert_that(max(E$node1) <= ncol(X) & max(E$node2) <= ncol(X))

print("done.")
