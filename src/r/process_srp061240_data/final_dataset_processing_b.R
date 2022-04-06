library(DESeq2)
library(magrittr)
library(SummarizedExperiment)

start_time <- Sys.time()

IN_DIR <- "/home/jgburk/PycharmProjects/reticula/data/SRP061240/input/"
OUT_DIR <- "/home/jgburk/PycharmProjects/reticula/data/SRP061240/output/"

srp061240.df <- readRDS(paste(OUT_DIR,"srp061240_df.Rds",sep=""))
srp061240.tissue.vec <- readRDS(paste(OUT_DIR,"srp061240_tissue_vec.Rds",sep=""))

#minimum shrinkage, leaving max() == integer max
scale.factor <- (.Machine$integer.max - 1) / max(srp061240.df)
srp061240.df <- round(srp061240.df * scale.factor)
srp061240.df <- srp061240.df + 1

dds <- DESeq2::DESeqDataSetFromMatrix(countData = as.matrix(srp061240.df),
                                      colData = data.frame(Sample=colnames(srp061240.df),
                                                           Tissue=srp061240.tissue.vec),
                                      design = ~ Tissue)
saveRDS(dds,
        paste(OUT_DIR,"dds.Rds",sep=""))

end_time <- Sys.time()
print(paste("Start: ",start_time," End: ",end_time," Difference: ",end_time - start_time))
