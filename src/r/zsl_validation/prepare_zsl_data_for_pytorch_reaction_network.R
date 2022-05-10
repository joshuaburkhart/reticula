library(magrittr)

OUT_DIR <- "/home/jgburk/PycharmProjects/reticula/data/tcga/output/"
IN_DIR <- "/home/jgburk/PycharmProjects/reticula/data/tcga/input/"
GTEX_DIR <- "/home/jgburk/PycharmProjects/reticula/data/gtex/"

X <- readRDS(paste(OUT_DIR, "rxn_pca_nls.Rds", sep = ""))
Y <- readRDS(paste(OUT_DIR,"tcga_tissue_vec_train.Rds",sep=""))
gtex_rxn2nodeLabel.df <- read.table(paste(GTEX_DIR,
                                           "Copy\ of\ rxn2nodeLabel_nls.csv",sep=""),
                                     sep=" ")

tissue_alphabetical_order <- order(Y)

X <- as.data.frame(X) %>% .[tissue_alphabetical_order,]

Y <- as.data.frame(Y[tissue_alphabetical_order])

m <- setdiff(make.names(gtex_rxn2nodeLabel.df$V1),colnames(X)) 
X[m] <- 0
X <- X[make.names(gtex_rxn2nodeLabel.df$V1)]

write.table(X,
          file=paste(IN_DIR,"node_features.txt",sep=""),
          row.names = FALSE,
          col.names = FALSE)
write.table(Y,
          file=paste(IN_DIR,"graph_targets.txt",sep=""),
          row.names = FALSE,
          col.names = FALSE)

call_check.df <- data.frame(Sample_ID = rownames(X),
                            Tissue = Y$`Y[tissue_alphabetical_order]`,
                            Dataset = rep(c("Tuning Set","Test Set"),length(X[,1])/2))
write.table(call_check.df,
            file=paste(OUT_DIR,"call_check.csv",sep=""),
            row.names = FALSE,
            col.names = TRUE,
            sep = ",")
