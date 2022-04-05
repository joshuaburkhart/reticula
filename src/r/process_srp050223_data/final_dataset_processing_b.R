library(DESeq2)
library(magrittr)
library(SummarizedExperiment)

start_time <- Sys.time()

IN_DIR <- "/home/jgburk/PycharmProjects/reticula/data/SRP050223/input/"
OUT_DIR <- "/home/jgburk/PycharmProjects/reticula/data/SRP050223/output/"

srp050223.df <- readRDS(paste(OUT_DIR,"srp050223_df.Rds",sep=""))
srp050223.tissue.vec <- readRDS(paste(OUT_DIR,"srp050223_tissue_vec.Rds",sep=""))

#minimum shrinkage, leaving max() == integer max
scale.factor <- (.Machine$integer.max - 1) / max(srp050223.df)
srp050223.df <- round(srp050223.df * scale.factor)
srp050223.df <- srp050223.df + 1

dds <- DESeq2::DESeqDataSetFromMatrix(countData = as.matrix(srp050223.df),
                                      colData = data.frame(Sample=colnames(srp050223.df),
                                                           Tissue=srp050223.tissue.vec),
                                      design = ~ Tissue)
saveRDS(dds,
        paste(OUT_DIR,"dds.Rds",sep=""))

end_time <- Sys.time()
print(paste("Start: ",start_time," End: ",end_time," Difference: ",end_time - start_time))
