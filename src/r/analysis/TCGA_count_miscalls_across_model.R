library(dplyr)
library(magrittr)

DATA_DIR <- "/home/jgburk/PycharmProjects/reticula/data/tcga/output/"

incorrect_calls_df <- read.table(paste(DATA_DIR,"res_gnn_incorrect.csv",sep=""),sep = ",",header = TRUE) 
colnames(incorrect_calls_df) <- c("Tissue", "Both_Misclass", "Only_Resnet_Misclass", "Only_GNN_Misclass")
n <- ncol(incorrect_calls_df)
for(i in 2:n){
  incorrect_calls_df[,i] <- as.numeric(incorrect_calls_df[,i])
}

incorrect_totals.df <- incorrect_calls_df %>% dplyr::group_by(Tissue) %>%
  dplyr::summarise(Both_Misclass_Total = sum(Both_Misclass),
                   Only_Resnet_Misclass_Total = sum(Only_Resnet_Misclass),
                   Only_GNN_Misclass_Total = sum(Only_GNN_Misclass))

print(incorrect_totals.df[,2:4] %>% colSums())
