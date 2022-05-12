library(magrittr)

OUT_DIR <- "/home/jgburk/PycharmProjects/reticula/data/tcga/output/"
IN_DIR <- "/home/jgburk/PycharmProjects/reticula/data/tcga/input/"

X <- readRDS(paste(OUT_DIR, "zsl_rxn_pca_nls.Rds", sep = ""))
Y <- readRDS(paste(OUT_DIR,"zsl_tcga_tissue_vec_train.Rds",sep=""))

tissue_alphabetical_order <- order(Y)

X <- as.data.frame(X) %>% .[tissue_alphabetical_order,]

Y <- as.data.frame(Y[tissue_alphabetical_order])

write.table(X,
          file=paste(IN_DIR,"zsl_node_features.txt",sep=""),
          row.names = FALSE,
          col.names = FALSE)
write.table(Y,
          file=paste(IN_DIR,"zsl_graph_targets.txt",sep=""),
          row.names = FALSE,
          col.names = FALSE)

call_check.df <- data.frame(Sample_ID = rownames(X),
                            Tissue = Y$`Y[tissue_alphabetical_order]`,
                            Dataset = rep(c("Tuning Set","Test Set"),length(X[,1])/2))
write.table(call_check.df,
            file=paste(OUT_DIR,"zsl_call_check.csv",sep=""),
            row.names = FALSE,
            col.names = TRUE,
            sep = ",")
