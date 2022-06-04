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

tissue_alphabetical_order <- order(Y)

X <- as.data.frame(X) %>% .[tissue_alphabetical_order,]

Y <- as.data.frame(Y[tissue_alphabetical_order])

write.table(X,
            file=paste(GTEX_IN_DIR,"zsl_gtex_node_features.txt",sep=""),
            row.names = FALSE,
            col.names = FALSE)
write.table(Y,
            file=paste(GTEX_IN_DIR,"zsl_gtex_graph_targets.txt",sep=""),
            row.names = FALSE,
            col.names = FALSE)
