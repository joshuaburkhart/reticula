library(DESeq2)
library(magrittr)
library(SummarizedExperiment)

start_time <- Sys.time()

IN_DIR <- "/home/jgburk/PycharmProjects/reticula/data/gtex/input/"
OUT_DIR <- "/home/jgburk/PycharmProjects/reticula/data/gtex/output/"

#GTEx_DATA_DIR <- "/home/users/burkhajo/WuLab/WuLabLustreDir/reticula/input/recount2/recount2edGTEx/"
#GTEx_DATA_DIR <- "/Users/burkhajo/Software/reticula/data/aim1/input/recount2/recount2edGTEx/"
#GTEx_DATA_DIR <- "/home/burkhart/Software/reticula/data/aim1/input/recount2/recountedGTEx/"
GTEx_DATA_DIR <- IN_DIR
GTEx_DATA_FIL <- "rse_gene.Rdata"

ensembl2rxns.df <- read.table(paste(IN_DIR,"Ensembl2ReactomeReactions.txt",sep=""),
                              sep="\t",
                              comment.char = "")

load(paste(GTEx_DATA_DIR,GTEx_DATA_FIL,sep=""))
rse_gene <- rse_gene[,!(rse_gene$smtsd == "Cells - Transformed fibroblasts" | # 306 of these samples
                          rse_gene$smtsd == "Cells - Leukemia cell line (CML)" | # 102 of these samples
                          rse_gene$smtsd == "Cells - EBV-transformed lymphocytes" | # 139 of these samples
                          rse_gene$smtsd == "" )] # 5 of these samples (details below) 

# > rse_gene[,rse_gene$smts == ""] %>% .$smtsd
# [1] "Esophagus - Mucosa"             "Skin - Sun Exposed (Lower leg)" "Stomach"                        "Skin - Sun Exposed (Lower leg)" "Esophagus - Mucosa"            
# > rse_gene[,rse_gene$smts == ""] %>% .$sample
# [1] "SRS638130" "SRS626174" "SRS626199" "SRS648268" "SRS637678"

gtex.cols <- rse_gene %>% colData()
saveRDS(gtex.cols$smtsd,file=paste(OUT_DIR,"gtex_tissue_detail_vec.Rds",sep=""))
gtex.tissue.detail.vec <- readRDS(paste(OUT_DIR,"gtex_tissue_detail_vec.Rds",sep=""))

gtex.df <- rse_gene %>% SummarizedExperiment::assay() %>% as.data.frame()
ensembl_wo_ids <- gsub("\\.[0-9]+","",rownames(gtex.df))
rownames(gtex.df) <- ensembl_wo_ids
reactome_ensembl_ids <- intersect(ensembl2rxns.df$V1,ensembl_wo_ids)
saveRDS(reactome_ensembl_ids,file=paste(OUT_DIR,"reactome_ensembl_ids.Rds",sep=""))
gtex.df <- gtex.df[reactome_ensembl_ids,]
saveRDS(gtex.df,file=paste(OUT_DIR,"gtex_df.Rds",sep=""))

end_time <- Sys.time()
print(paste("Start: ",start_time," End: ",end_time," Difference: ",end_time - start_time))