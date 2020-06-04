library(DESeq2)
library(magrittr)
library(SummarizedExperiment)

start_time <- Sys.time()

IN_DIR <- "/home/burkhart/Software/reticula/data/aim1/input/"
OUT_DIR <- "/home/burkhart/Software/reticula/data/aim1/output/"

#GTEx_DATA_DIR <- "/home/users/burkhajo/WuLab/WuLabLustreDir/reticula/input/recount2/recount2edGTEx/"
#GTEx_DATA_DIR <- "/Users/burkhajo/Software/reticula/data/aim1/input/recount2/recount2edGTEx/"
GTEx_DATA_DIR <- "/home/burkhart/Software/reticula/data/aim1/input/recount2/recountedGTEx/"
GTEx_DATA_FIL <- "rse_gene.Rdata"

gtex.df <- readRDS(paste(OUT_DIR,"gtex_df.Rds",sep=""))
gtex.tissue.detail.vec <- readRDS(paste(OUT_DIR,"gtex_tissue_detail_vec.Rds",sep=""))

#minimum shrinkage, leaving max() == integer max
scale.factor <- (.Machine$integer.max - 1) / max(gtex.df)
gtex.df <- round(gtex.df * scale.factor)
gtex.df <- gtex.df + 1

dds <- DESeq2::DESeqDataSetFromMatrix(countData = as.matrix(gtex.df),
                                      colData = data.frame(Sample=colnames(gtex.df),
                                                           Tissue=gtex.tissue.detail.vec),
                                      design = ~ Tissue)
saveRDS(dds,
        paste(OUT_DIR,"dds.Rds",sep=""))

end_time <- Sys.time()
print(paste("Start: ",start_time," End: ",end_time," Difference: ",end_time - start_time))
