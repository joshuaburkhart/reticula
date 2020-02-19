library(DESeq2)

combined.df <- readRDS("~/combined_df.Rds")
tissue.vec <- readRDS("~/tissue_vec.Rds")
datasource.vec <- readRDS("~/datasource_vec.Rds")

#1.5 is a minimum shrinkage, leaving max() == integer max
scale.factor <- (.Machine$integer.max - .Machine$double.eps) / max(combined.df)
combined.scaled.df <- round(combined.df * scale.factor)
combined.nozero.df <- combined.scaled.df + .Machine$double.eps

dds <- DESeq2::DESeqDataSetFromMatrix(countData = as.matrix(combined.nozero.df),
                                      colData = data.frame(Sample=colnames(combined.df),
                                                           Tissue=tissue.vec,
                                                           Datasource=datasource.vec),
                                      design = ~ Tissue + Datasource)
saveRDS(dds,
        "~/dds.Rds")

vst.counts <- DESeq2::vst(dds,
                          blind = FALSE,
                          fitType = "local")

saveRDS(vst.counts,
        "~/vst_counts.Rds")
