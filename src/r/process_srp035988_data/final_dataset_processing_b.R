library(DESeq2)
library(magrittr)
library(SummarizedExperiment)

start_time <- Sys.time()

IN_DIR <- "/home/jgburk/PycharmProjects/reticula/data/SRP035988/input/"
OUT_DIR <- "/home/jgburk/PycharmProjects/reticula/data/SRP035988/output/"

srp035988.df <- readRDS(paste(OUT_DIR,"srp035988_df.Rds",sep=""))
srp035988.tissue.vec <- readRDS(paste(OUT_DIR,"srp035988_tissue_vec.Rds",sep=""))

#minimum shrinkage, leaving max() == integer max
scale.factor <- (.Machine$integer.max - 1) / max(srp035988.df)
srp035988.df <- round(srp035988.df * scale.factor)
srp035988.df <- srp035988.df + 1

dds <- DESeq2::DESeqDataSetFromMatrix(countData = as.matrix(srp035988.df),
                                      colData = data.frame(Sample=colnames(srp035988.df),
                                                           Tissue=srp035988.tissue.vec),
                                      design = ~ Tissue)
saveRDS(dds,
        paste(OUT_DIR,"dds.Rds",sep=""))

end_time <- Sys.time()
print(paste("Start: ",start_time," End: ",end_time," Difference: ",end_time - start_time))
