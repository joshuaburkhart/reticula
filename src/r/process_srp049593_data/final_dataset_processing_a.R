library(DESeq2)
library(magrittr)
library(SummarizedExperiment)

start_time <- Sys.time()

IN_DIR <- "/home/jgburk/PycharmProjects/reticula/data/SRP049593/input/"
OUT_DIR <- "/home/jgburk/PycharmProjects/reticula/data/SRP049593/output/"

SRP049593_DATA_FIL <- "rse_gene(1).Rdata"

ensembl2rxns.df <- read.table(paste(IN_DIR,"Ensembl2ReactomeReactions.txt",sep=""),
                              sep="\t",
                              comment.char = "")

load(paste(IN_DIR,SRP049593_DATA_FIL,sep=""))

g <- character()
for(i in 1:length(rse_gene$characteristics)){
  g <- c(g,unlist(strsplit(rse_gene$characteristics %>% .[[i]] %>% .[1],": "))[2]) #1 = "genotype"
}

assertthat::are_equal(length(g),length(rse_gene$characteristics))

saveRDS(g,file=paste(OUT_DIR,"srp049593_tissue_vec.Rds",sep=""))
srp049593.tissue.vec <- readRDS(paste(OUT_DIR,"srp049593_tissue_vec.Rds",sep=""))

srp049593.df <- rse_gene %>% SummarizedExperiment::assay() %>% as.data.frame()
ensembl_wo_ids <- gsub("\\.[0-9]+","",rownames(srp049593.df))
rownames(srp049593.df) <- ensembl_wo_ids
reactome_ensembl_ids <- intersect(ensembl2rxns.df$V1,ensembl_wo_ids)
saveRDS(reactome_ensembl_ids,file=paste(OUT_DIR,"reactome_ensembl_ids.Rds",sep=""))
srp049593.df <- srp049593.df[reactome_ensembl_ids,]
saveRDS(srp049593.df,file=paste(OUT_DIR,"srp049593_df.Rds",sep=""))

end_time <- Sys.time()
print(paste("Start: ",start_time," End: ",end_time," Difference: ",end_time - start_time))
