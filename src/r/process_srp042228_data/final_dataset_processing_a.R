library(DESeq2)
library(magrittr)
library(SummarizedExperiment)

start_time <- Sys.time()

IN_DIR <- "/home/jgburk/PycharmProjects/reticula/data/SRP042228/input/"
OUT_DIR <- "/home/jgburk/PycharmProjects/reticula/data/SRP042228/output/"

SRP042228_DATA_FIL <- "rse_gene.Rdata"

ensembl2rxns.df <- read.table(paste(IN_DIR,"Ensembl2ReactomeReactions.txt",sep=""),
                              sep="\t",
                              comment.char = "")

load(paste(IN_DIR,SRP042228_DATA_FIL,sep=""))

g <- character()
for(i in 1:length(rse_gene$characteristics)){
  g <- c(g,unlist(strsplit(rse_gene$characteristics %>% .[[i]] %>% .[6],": "))[2]) #6 = "I2 type"
}

healthy_and_colorectal_sample_idxs <- sort(c(which(g == "iCD"),
                                             which(g == "cCD"),
                                             which(g == "Not IBD"),
                                             which(g == "not IBD"),
                                             which(g == "UC")))

#healthy_and_colorectal_sample_idxs <- setdiff(healthy_and_colorectal_sample_idxs,which(grepl("Sample_1S",rse_gene$title))) # Keep only colorectal stages II-IV

g <- g[healthy_and_colorectal_sample_idxs]

stopifnot(assertthat::are_equal(length(g),314))

saveRDS(g,file=paste(OUT_DIR,"srp042228_tissue_vec.Rds",sep=""))
srp042228.tissue.vec <- readRDS(paste(OUT_DIR,"srp042228_tissue_vec.Rds",sep=""))

srp042228_healthy_and_colorectal.df <- rse_gene %>% SummarizedExperiment::assay() %>% as.data.frame() %>% .[,healthy_and_colorectal_sample_idxs]
ensembl_wo_ids <- gsub("\\.[0-9]+","",rownames(srp042228_healthy_and_colorectal.df))
rownames(srp042228_healthy_and_colorectal.df) <- ensembl_wo_ids
reactome_ensembl_ids <- intersect(ensembl2rxns.df$V1,ensembl_wo_ids)
saveRDS(reactome_ensembl_ids,file=paste(OUT_DIR,"reactome_ensembl_ids.Rds",sep=""))
srp042228_healthy_and_colorectal.df <- srp042228_healthy_and_colorectal.df[reactome_ensembl_ids,]
saveRDS(srp042228_healthy_and_colorectal.df,file=paste(OUT_DIR,"srp042228_df.Rds",sep=""))

end_time <- Sys.time()
print(paste("Start: ",start_time," End: ",end_time," Difference: ",end_time - start_time))
