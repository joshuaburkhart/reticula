library(DESeq2)
library(magrittr)
library(SummarizedExperiment)

start_time <- Sys.time()

IN_DIR <- "/home/jgburk/PycharmProjects/reticula/data/SRP050223/input/"
OUT_DIR <- "/home/jgburk/PycharmProjects/reticula/data/SRP050223/output/"

SRP050223_DATA_FIL <- "rse_gene(2).Rdata"

ensembl2rxns.df <- read.table(paste(IN_DIR,"Ensembl2ReactomeReactions.txt",sep=""),
                              sep="\t",
                              comment.char = "")

load(paste(IN_DIR,SRP050223_DATA_FIL,sep=""))

saveRDS(rse_gene$characteristics %>% as.data.frame() %>% .$value,file=paste(OUT_DIR,"srp050223_tissue_vec.Rds",sep=""))
srp050223.tissue.vec <- readRDS(paste(OUT_DIR,"srp050223_tissue_vec.Rds",sep=""))

srp050223.df <- rse_gene %>% SummarizedExperiment::assay() %>% as.data.frame()
ensembl_wo_ids <- gsub("\\.[0-9]+","",rownames(srp050223.df))
rownames(srp050223.df) <- ensembl_wo_ids
reactome_ensembl_ids <- intersect(ensembl2rxns.df$V1,ensembl_wo_ids)
saveRDS(reactome_ensembl_ids,file=paste(OUT_DIR,"reactome_ensembl_ids.Rds",sep=""))
srp050223.df <- srp050223.df[reactome_ensembl_ids,]
saveRDS(srp050223.df,file=paste(OUT_DIR,"srp050223_df.Rds",sep=""))

end_time <- Sys.time()
print(paste("Start: ",start_time," End: ",end_time," Difference: ",end_time - start_time))
