library(DESeq2)
library(magrittr)
library(SummarizedExperiment)

start_time <- Sys.time()

TCGA_IN_DIR <- "/home/jgburk/PycharmProjects/reticula/data/tcga/input/"
TCGA_OUT_DIR <- "/home/jgburk/PycharmProjects/reticula/data/tcga/output/"
GTEX_IN_DIR <- "/home/jgburk/PycharmProjects/reticula/data/gtex/input/"
GTEX_OUT_DIR <- "/home/jgburk/PycharmProjects/reticula/data/gtex/output/"

gtex.df <- readRDS(paste(GTEX_OUT_DIR,"gtex_df.Rds",sep=""))
gtex.tissue.detail.vec <- readRDS(paste(GTEX_OUT_DIR,"gtex_tissue_detail_vec.Rds",sep=""))

gtex.scale.factor <- (.Machine$integer.max - 1) / max(gtex.df)
gtex.df <- round(gtex.df * gtex.scale.factor)
gtex.df <- gtex.df + 1.0

dds <- DESeq2::DESeqDataSetFromMatrix(countData = as.matrix(gtex.df),
                                      colData = data.frame(Sample=colnames(gtex.df),
                                                           Tissue=gtex.tissue.detail.vec),
                                      design = ~ Tissue)
saveRDS(dds,
        paste(GTEX_OUT_DIR,"sf_dds.Rds",sep=""))

tcga.df <- readRDS(paste(TCGA_OUT_DIR,"tcga_df.Rds",sep=""))
tcga.tissue.vec <- readRDS(paste(TCGA_OUT_DIR,"tcga_tissue_vec.Rds",sep=""))

tcga.df <- round(tcga.df * gtex.scale.factor)
tcga.df <- tcga.df + 1.0

dds <- DESeq2::DESeqDataSetFromMatrix(countData = as.matrix(tcga.df),
                                      colData = data.frame(Sample=colnames(tcga.df),
                                                           Tissue=tcga.tissue.vec),
                                      design = ~ Tissue)
saveRDS(dds,
        paste(TCGA_OUT_DIR,"sf_dds.Rds",sep=""))

end_time <- Sys.time()
print(paste("Start: ",start_time," End: ",end_time," Difference: ",end_time - start_time))
