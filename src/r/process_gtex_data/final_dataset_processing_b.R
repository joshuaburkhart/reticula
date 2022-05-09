library(DESeq2)
library(magrittr)
library(SummarizedExperiment)

start_time <- Sys.time()

IN_DIR <- "/home/jgburk/PycharmProjects/reticula/data/gtex/input/"
OUT_DIR <- "/home/jgburk/PycharmProjects/reticula/data/gtex/output/"

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
