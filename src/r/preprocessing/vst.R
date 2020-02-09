library(DESeq2)

combined.df <- readRDS("~/combined_df.Rds")
tissue.vec <- readRDS("~/tissue_vec.Rds")
datasource.vec <- readRDS("~/datasource_vec.Rds")

#1.5 is a minimum shrinkage, leaving max() == integer max
combined.scaled.df <- ceiling(combined.df * ((.Machine$integer.max - 1.5)/max(combined.df)))
combined.nozero.df <- combined.scaled.df + 1

#number of 1's in nozero should match number of 0's in scaled
sum(combined.scaled.df == 0) == sum(combined.nozero.df == 1)

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
