library(magrittr)

OUT_DIR <- "/home/jgburk/PycharmProjects/reticula/data/SRP042228/output/"
IN_DIR <- "/home/jgburk/PycharmProjects/reticula/data/SRP042228/input/"
GTEX_DIR <- "/home/jgburk/PycharmProjects/reticula/data/gtex/"

X <- readRDS(paste(OUT_DIR, "rxn_pca_nls.Rds", sep = ""))
Y <- readRDS(paste(OUT_DIR,"srp042228_tissue_vec_train.Rds",sep=""))
gtex_rxn2nodeLabel.df <- read.table(paste(GTEX_DIR,
                                           "Copy\ of\ rxn2nodeLabel_nls.csv",sep=""),
                                     sep=" ")

X <- as.data.frame(X)

Y <- as.data.frame(Y)

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
