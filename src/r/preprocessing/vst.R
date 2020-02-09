library(DESeq2)

combined.df <- readRDS("~/combined_df.Rds")
tissue.vec <- readRDS("~/tissue_vec.Rds")
datasource.vec <- readRDS("~/datasource_vec.Rds")

combined.scaled.df <- floor(combined.df * ((.Machine$integer.max - 1)/max(combined.df)))
combined.nozero.df <- combined.scaled.df + 1

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
