library(DESeq2)
library(magrittr)
library(SummarizedExperiment)

start_time <- Sys.time()

IN_DIR <- "/home/jgburk/PycharmProjects/reticula/data/tcga/input/"
OUT_DIR <- "/home/jgburk/PycharmProjects/reticula/data/tcga/output/"

tcga.df <- readRDS(paste(OUT_DIR,"tcga_df.Rds",sep=""))
tcga.tissue.vec <- readRDS(paste(OUT_DIR,"tcga_tissue_vec.Rds",sep=""))

#minimum shrinkage, leaving max() == integer max
scale.factor <- (.Machine$integer.max - 1) / max(tcga.df)
tcga.df <- round(tcga.df * scale.factor)
tcga.df <- tcga.df + 1

dds <- DESeq2::DESeqDataSetFromMatrix(countData = as.matrix(tcga.df),
                                      colData = data.frame(Sample=colnames(tcga.df),
                                                           Tissue=tcga.tissue.vec),
                                      design = ~ Tissue)
saveRDS(dds,
        paste(OUT_DIR,"dds.Rds",sep=""))

end_time <- Sys.time()
print(paste("Start: ",start_time," End: ",end_time," Difference: ",end_time - start_time))
