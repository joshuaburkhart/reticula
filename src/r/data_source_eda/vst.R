library(DESeq2)

combined.df <- readRDS("~/combined_df.Rds")
tissue.vec <- readRDS("~/tissue_vec.Rds")
datasource.vec <- readRDS("~/datasource_vec.Rds")
study.vec <- readRDS("~/study_vec.Rds")

#minimum shrinkage, leaving max() == integer max
scale.factor <- (.Machine$integer.max - 1) / max(combined.df)
combined.scaled.df <- round(combined.df * scale.factor)
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

vst.counts <- readRDS("~/vst_counts.Rds")

vst.mtx.counts <- as.matrix(assay(vst.counts))

rbe.design <- model.matrix(~tissue.vec)

vst.rbe <- limma::removeBatchEffect(vst.mtx.counts,
                                    batch=datasource.vec,
                                    design=rbe.design)
saveRDS(vst.rbe,
        "~/vst_rbe.Rds")

vst.rbe.counts <- vst.counts
assay(vst.rbe.counts) <- vst.rbe

DESeq2::plotPCA(vst.counts,intgroup="Tissue")
DESeq2::plotPCA(vst.rbe.counts,intgroup="Tissue")

# umap takes a long time and isn't necessary here
#umap.com <- umap::umap(t(vst.mtx.counts))
#plot(umap.com$layout,col=as.numeric(as.factor(datasource.vec)),main="raw")

#umap.rbe <- umap::umap(t(vst.rbe))
#plot(umap.rbe$layout,col=as.numeric(as.factor(datasource.vec)),main="raw")

filter <- tissue.vec == "ESGS" |
  tissue.vec == "Esophagus - Gastroesophageal Junction" |
  tissue.vec == "Esophagus - Mucosa" |
  tissue.vec == "Esophagus - Muscularis"
f.tissue <- tissue.vec[filter]
f.datasource <- datasource.vec[filter]
f.study <- study.vec[filter]

dds.f <- DESeq2::DESeqDataSetFromMatrix(as.matrix(combined.nozero.df[,filter]),
                                        colData = data.frame(Sample=colnames(combined.df[,filter]),
                                                             Datasource=f.datasource,
                                                             Study=f.study,
                                                             Tissue=f.tissue),
                                        design = ~Datasource)

vst.counts.f <- DESeq2::vst(dds.f,
                            blind=FALSE,
                            fitType="local")

DESeq2::plotPCA(vst.counts.f,intgroup="Datasource")
DESeq2::plotPCA(vst.counts.f,intgroup="Tissue")
DESeq2::plotPCA(vst.counts.f,intgroup="Study")

vst.mtx.counts.f <- as.matrix(assay(vst.counts.f))
vst.rbe.f <- limma::removeBatchEffect(vst.mtx.counts.f,
                                      batch = f.datasource)
assay(vst.counts.f) <- vst.rbe.f

DESeq2::plotPCA(vst.counts.f,intgroup="Datasource")
DESeq2::plotPCA(vst.counts.f,intgroup="Tissue")
DESeq2::plotPCA(vst.counts.f,intgroup="Study")
