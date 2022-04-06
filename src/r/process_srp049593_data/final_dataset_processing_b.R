library(DESeq2)
library(magrittr)
library(SummarizedExperiment)

start_time <- Sys.time()

IN_DIR <- "/home/jgburk/PycharmProjects/reticula/data/SRP049593/input/"
OUT_DIR <- "/home/jgburk/PycharmProjects/reticula/data/SRP049593/output/"

srp049593.df <- readRDS(paste(OUT_DIR,"srp049593_df.Rds",sep=""))
srp049593.tissue.vec <- readRDS(paste(OUT_DIR,"srp049593_tissue_vec.Rds",sep=""))

#minimum shrinkage, leaving max() == integer max
scale.factor <- (.Machine$integer.max - 1) / max(srp049593.df)
srp049593.df <- round(srp049593.df * scale.factor)
srp049593.df <- srp049593.df + 1

dds <- DESeq2::DESeqDataSetFromMatrix(countData = as.matrix(srp049593.df),
                                      colData = data.frame(Sample=colnames(srp049593.df),
                                                           Tissue=srp049593.tissue.vec),
                                      design = ~ Tissue)
saveRDS(dds,
        paste(OUT_DIR,"dds.Rds",sep=""))

end_time <- Sys.time()
print(paste("Start: ",start_time," End: ",end_time," Difference: ",end_time - start_time))
