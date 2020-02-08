library(DESeq2)

combined.df <- readRDS("~/combined_df.Rds")
tissue.vec <- readRDS("~/tissue_vec.Rds")
datasource.vec <- readRDS("~/datasource_vec.Rds")

dds <- DESeq2::DESeqDataSetFromMatrix(countData = as.matrix(combined.df),
                                      colData = colnames(combined.df),
                                      design = ~ datasource.vec + tissue.vec)
saveRDS(dds,
        "~/dds.Rds")
rlog_counts <- DESeq2::rlog(dds)
saveRDS(rlog_counts,
        "~/rlog_counts.Rds")