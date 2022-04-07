library(DESeq2)
library(magrittr)
library(SummarizedExperiment)

start_time <- Sys.time()

IN_DIR <- "/home/jgburk/PycharmProjects/reticula/data/SRP042228/input/"
OUT_DIR <- "/home/jgburk/PycharmProjects/reticula/data/SRP042228/output/"

srp042228.df <- readRDS(paste(OUT_DIR,"srp042228_df.Rds",sep=""))
srp042228.tissue.vec <- readRDS(paste(OUT_DIR,"srp042228_tissue_vec.Rds",sep=""))

#minimum shrinkage, leaving max() == integer max
scale.factor <- (.Machine$integer.max - 1) / max(srp042228.df)
srp042228.df <- round(srp042228.df * scale.factor)
srp042228.df <- srp042228.df + 1

dds <- DESeq2::DESeqDataSetFromMatrix(countData = as.matrix(srp042228.df),
                                      colData = data.frame(Sample=colnames(srp042228.df),
                                                           Tissue=srp042228.tissue.vec),
                                      design = ~ Tissue)
saveRDS(dds,
        paste(OUT_DIR,"dds.Rds",sep=""))

end_time <- Sys.time()
print(paste("Start: ",start_time," End: ",end_time," Difference: ",end_time - start_time))
