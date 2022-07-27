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

# Should match manuscript Table S2

gene_name_labs <- rownames(res)
gene_name_labs[which(gene_name_labs == "ENSG00000114251")] <- "WNT5A"
gene_name_labs[which(gene_name_labs == "ENSG00000050730")] <- "TNIP3"
gene_name_labs[which(gene_name_labs == "ENSG00000113749")] <- "HRH2"
gene_name_labs[which(gene_name_labs == "ENSG00000172350")] <- "ABCG4"
gene_name_labs[which(gene_name_labs == "ENSG00000174808")] <- "BTC"
gene_name_labs[which(gene_name_labs == "ENSG00000182481")] <- "KPNA2"
gene_name_labs[which(gene_name_labs == "ENSG00000134827")] <- "TCN1"
gene_name_labs[which(gene_name_labs == "ENSG00000105173")] <- "CCNE1"
gene_name_labs[which(gene_name_labs == "ENSG00000136695")] <- "IL36RN"
gene_name_labs[which(gene_name_labs == "ENSG00000134755")] <- "DSC2"
gene_name_labs[which(gene_name_labs == "ENSG00000007171")] <- "NOS2"
gene_name_labs[which(gene_name_labs == "ENSG00000092295")] <- "TGM1"
gene_name_labs[which(gene_name_labs == "ENSG00000143365")] <- "RORC"
gene_name_labs[which(gene_name_labs == "ENSG00000100342")] <- "APOL1"
gene_name_labs[which(gene_name_labs == "ENSG00000177257")] <- "DEFB4B"
gene_name_labs[which(gene_name_labs == "ENSG00000049249")] <- "TNFRSF9"
gene_name_labs[which(gene_name_labs == "ENSG00000057657")] <- "PRDM1"
gene_name_labs[which(gene_name_labs == "ENSG00000185962")] <- "LCE3A"
gene_name_labs[which(gene_name_labs == "ENSG00000007402")] <- "CACNA2D2"

svg(filename = paste(OUT_DIR,"EnhancedVolcano.svg",sep=""), width = 11, height = 10)
EnhancedVolcano(res,
                lab=gene_name_labs,
                selectLab=c("WNT5A",
                            "TNIP3",
                            "HRH2",
                            "ABCG4",
                            "BTC",
                            "KPNA2",
                            "TCN1",
                            "CCNE1",
                            "IL36RN",
                            "DSC2",
                            "NOS2",
                            "TGM1",
                            "RORC",
                            "APOL1",
                            "DEFB4B",
                            "TNFRSF9",
                            "PRDM1",
                            "LCE3A",
                            "CACNA2D2"), 
                x='log2FoldChange',
                y='pvalue')
dev.off()

# Should match manuscript Table S3

gene_name_labs <- rownames(res)
gene_name_labs[which(gene_name_labs == "ENSG00000108671")] <- "PSMD11"
gene_name_labs[which(gene_name_labs == "ENSG00000143106")] <- "PSMA5"
gene_name_labs[which(gene_name_labs == "ENSG00000092010")] <- "PSME1"
gene_name_labs[which(gene_name_labs == "ENSG00000126067")] <- "PSMB2"
gene_name_labs[which(gene_name_labs == "ENSG00000100567")] <- "PSMA3"
gene_name_labs[which(gene_name_labs == "ENSG00000173692")] <- "PSMD1"
gene_name_labs[which(gene_name_labs == "ENSG00000163636")] <- "PSMD6"
gene_name_labs[which(gene_name_labs == "ENSG00000142507")] <- "PSMB6"

svg(filename = paste(OUT_DIR,"EnhancedVolcano_CentralReaction.svg",sep=""), width = 11, height = 10)
EnhancedVolcano(res,
                lab=gene_name_labs,
                selectLab=c("PSMD11",
                            "PSMA5",
                            "PSME1",
                            "PSMB2",
                            "PSMA3",
                            "PSMD1",
                            "PSMD6",
                            "PSMB6"), 
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
