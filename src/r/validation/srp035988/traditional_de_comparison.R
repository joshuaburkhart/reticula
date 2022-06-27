library(DESeq2)
library(magrittr)
library(EnhancedVolcano)

ALPHA <- 0.05
OUT_DIR <- "/home/jgburk/PycharmProjects/reticula/data/SRP035988/output/"

dds <- readRDS(paste(OUT_DIR,"dds.Rds",sep=""))

dds$Tissue <- relevel(dds$Tissue,
                      ref= "normal skin")
dds_de <- DESeq(dds,
                betaPrior = FALSE)
res <- results(dds_de,
               contrast = c('Tissue',
                            'lesional psoriatic skin',
                            'normal skin'))
res_lfcSh <- lfcShrink(dds_de,
                       contrast = c('Tissue',
                                    'lesional psoriatic skin',
                                    'normal skin'),
                       res=res,
                       type='normal')
svg(filename = paste(OUT_DIR,"EnhancedVolcano.svg",sep=""), width = 11, height = 10)
EnhancedVolcano(res,
                lab=rownames(res),
                x='log2FoldChange',
                y='pvalue')
dev.off()

res.df <- res_lfcSh[order(res_lfcSh$padj),] %>% as.data.frame()

res.df$Sig_W_Default_Thresholds = ((res.df$log2FoldChange > 2) & (res.df$pvalue < 10e-6))
res.df$Sig_W_Strict_Thresholds = ((res.df$log2FoldChange > 2) & (res.df$padj < 10e-32))
res.df$Sig_W_Relaxed_Thresholds = ((res.df$log2FoldChange > 0.5) & (res.df$padj < ALPHA))

res.df$EnsemblID <- rownames(res.df)

res.df <- res.df[,c("EnsemblID","log2FoldChange","pvalue","padj","Sig_W_Default_Thresholds","Sig_W_Strict_Thresholds","Sig_W_Relaxed_Thresholds")]

write.table(res.df, file=paste(OUT_DIR,"res_df.csv",sep=""),row.names = FALSE, quote = FALSE)
