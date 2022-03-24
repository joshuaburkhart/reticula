library(DESeq2)
library(magrittr)
library(SummarizedExperiment)

start_time <- Sys.time()

IN_DIR <- "/home/jgburk/PycharmProjects/reticula/data/tcga/input/"
OUT_DIR <- "/home/jgburk/PycharmProjects/reticula/data/tcga/output/"

TCGA_DATA_FIL <- "rse_gene(4).Rdata"

ensembl2rxns.df <- read.table(paste(IN_DIR,"Ensembl2ReactomeReactions.txt",sep=""),
                              sep="\t",
                              comment.char = "")

load(paste(IN_DIR,TCGA_DATA_FIL,sep=""))
rse_gene <- rse_gene[,(rse_gene$gdc_cases.samples.sample_type == "Solid Tissue Normal")]

tcga.cols <- rse_gene %>% colData()
saveRDS(tcga.cols$gdc_cases.project.primary_site,file=paste(OUT_DIR,"tcga_tissue_vec.Rds",sep=""))
tcga.tissue.detail.vec <- readRDS(paste(OUT_DIR,"tcga_tissue_vec.Rds",sep=""))

tcga.df <- rse_gene %>% SummarizedExperiment::assay() %>% as.data.frame()
ensembl_wo_ids <- gsub("\\.[0-9]+","",rownames(tcga.df))
rownames(tcga.df) <- ensembl_wo_ids
reactome_ensembl_ids <- intersect(ensembl2rxns.df$V1,ensembl_wo_ids)
saveRDS(reactome_ensembl_ids,file=paste(OUT_DIR,"reactome_ensembl_ids.Rds",sep=""))
tcga.df <- tcga.df[reactome_ensembl_ids,]
saveRDS(tcga.df,file=paste(OUT_DIR,"tcga_df.Rds",sep=""))

end_time <- Sys.time()
print(paste("Start: ",start_time," End: ",end_time," Difference: ",end_time - start_time))
