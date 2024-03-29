library(DESeq2)
library(magrittr)
library(EnhancedVolcano)

ALPHA <- 0.05
OUT_DIR <- "/home/jgburk/PycharmProjects/reticula/data/SRP050223/output/"

dds <- readRDS(paste(OUT_DIR,"dds.Rds",sep=""))

dds$Tissue <- relevel(dds$Tissue,
                      ref= "tissue: normal thymus")
dds_de <- DESeq(dds,
                betaPrior = FALSE)
res <- results(dds_de,
               contrast = c('Tissue',
                            'tissue: T cell acute lymphoblastic leukemia',
                            'tissue: normal thymus'))
res_lfcSh <- lfcShrink(dds_de,
                       contrast = c('Tissue',
                                    'tissue: T cell acute lymphoblastic leukemia',
                                    'tissue: normal thymus'),
                       res=res,
                       type='normal')
EnhancedVolcano(res,
                lab=rownames(res),
                x='log2FoldChange',
                y='pvalue')

res.df <- res_lfcSh[order(res_lfcSh$padj),] %>% .[.$padj < ALPHA,] %>% as.data.frame()

res.df$EnsemblID <- rownames(res.df)

res.df <- res.df[,c("EnsemblID","log2FoldChange","pvalue","padj")]

write.table(res.df, file=paste(OUT_DIR,"res_df.csv",sep=""),quote = FALSE)
