library(DESeq2)
library(magrittr)
library(SummarizedExperiment)

start_time <- Sys.time()

IN_DIR <- "/home/jgburk/PycharmProjects/reticula/data/SRP061240/input/"
OUT_DIR <- "/home/jgburk/PycharmProjects/reticula/data/SRP061240/output/"

SRP061240_DATA_FIL <- "rse_gene.Rdata"

ensembl2rxns.df <- read.table(paste(IN_DIR,"Ensembl2ReactomeReactions.txt",sep=""),
                              sep="\t",
                              comment.char = "")

load(paste(IN_DIR,SRP061240_DATA_FIL,sep=""))

"_1S" %in% rse_gene$title

g <- character()
for(i in 1:length(rse_gene$characteristics)){
  g <- c(g,unlist(strsplit(rse_gene$characteristics %>% .[[i]] %>% .[4],": "))[2]) #3 = "disease type"
}

healthy_and_colorectal_sample_idxs <- sort(c(which(g == "Healthy Control"),which(g == "Prostate Cancer")))#"Colorectal Cancer")))

#healthy_and_colorectal_sample_idxs <- setdiff(healthy_and_colorectal_sample_idxs,which(grepl("Sample_1S",rse_gene$title))) # Keep only colorectal stages II-IV

g <- g[healthy_and_colorectal_sample_idxs]

assertthat::are_equal(length(g),172) # 100 Healthy + 72 Prostate Samples

saveRDS(g,file=paste(OUT_DIR,"srp061240_tissue_vec.Rds",sep=""))
srp061240.tissue.vec <- readRDS(paste(OUT_DIR,"srp061240_tissue_vec.Rds",sep=""))

srp061240_healthy_and_colorectal.df <- rse_gene %>% SummarizedExperiment::assay() %>% as.data.frame() %>% .[,healthy_and_colorectal_sample_idxs]
ensembl_wo_ids <- gsub("\\.[0-9]+","",rownames(srp061240_healthy_and_colorectal.df))
rownames(srp061240_healthy_and_colorectal.df) <- ensembl_wo_ids
reactome_ensembl_ids <- intersect(ensembl2rxns.df$V1,ensembl_wo_ids)
saveRDS(reactome_ensembl_ids,file=paste(OUT_DIR,"reactome_ensembl_ids.Rds",sep=""))
srp061240_healthy_and_colorectal.df <- srp061240_healthy_and_colorectal.df[reactome_ensembl_ids,]
saveRDS(srp061240_healthy_and_colorectal.df,file=paste(OUT_DIR,"srp061240_df.Rds",sep=""))

end_time <- Sys.time()
print(paste("Start: ",start_time," End: ",end_time," Difference: ",end_time - start_time))
