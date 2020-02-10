library(snm)
library(umap)
library(limma)
library(phateR)
library(Biobase)
library(magrittr)
library(SummarizedExperiment)

tissue.vec <- readRDS("~/tissue_vec.Rds")
datasource.vec <- readRDS("~/datasource_vec.Rds")

# combined
combined.df <- readRDS("~/combined_df.Rds")

umap.com <- umap::umap(t(combined.df))
plot(umap.com$layout,col=as.numeric(as.factor(tissue.vec)),main="raw")
plot(umap.com$layout,col=as.numeric(as.factor(datasource.vec)),main="raw")

# vst
vst.counts <- readRDS("~/vst_counts.Rds")

umap.vst <- umap::umap(t(assay(vst.counts)))
plot(umap.vst$layout,col=as.numeric(as.factor(tissue.vec)),main="vst")
plot(umap.vst$layout,col=as.numeric(as.factor(datasource.vec)),main="vst")

# rbe vst
vst.rbe <- readRDS("~/vst_rbe.Rds")

umap.rbe <- umap::umap(t(vst.rbe))
plot(umap.rbe$layout,col=as.numeric(as.factor(tissue.vec)),main="rbe vst")
plot(umap.rbe$layout,col=as.numeric(as.factor(datasource.vec)),main="rbe vst")

# sqrt
sqrt.combined.mtx <- as.matrix(sqrt(combined.df))

umap.sqrt <- umap::umap(t(sqrt.combined.mtx))
plot(umap.sqrt$layout,col=as.numeric(as.factor(tissue.vec)),main="sqrt")
plot(umap.sqrt$layout,col=as.numeric(as.factor(datasource.vec)),main="sqrt")

# snm sqrt
snmR.sqrt.cad <- readRDS("~/snmR_sqrt_cad.Rds")

umap.snm <- umap::umap(t(snmR.sqrt.cad$norm.dat))
plot(umap.snm$layout,col=as.numeric(as.factor(tissue.vec)),main = "snm sqrt")
plot(umap.snm$layout,col=as.numeric(as.factor(datasource.vec)),main = "snm sqrt")
date()
Sys.Date()
snmR.sqrt.cad <- snm(sqrt.combined.mtx,
                bio.var=data.frame(tissue=tissue.vec),
                adj.var=data.frame(datasource=datasource.vec),
                rm.adj=TRUE,
                nbins=5,
                lmer.max.iter=10,
                spline.dim=3,
                num.iter=3)

saveRDS(snmR.sqrt.cad,
        "~/snmR_sqrt_cad.Rds")

vst.mtx.counts <- as.matrix(assay(vst.counts))

rbe.design <- model.matrix(~tissue.vec)

vst.rbe <- limma::removeBatchEffect(vst.mtx.counts,
                                    batch=datasource.vec,
                                    design=rbe.design)
saveRDS(vst.rbe,
        "~/vst_rbe.Rds")

vst.rbe <- readRDS("~/vst_rbe.Rds")

vst.rbe.t <- t(vst.rbe)

# rm(vst.rbe)

# vst + rbe
phate.data <- phateR::phate(vst.rbe.t,
                            n.jobs=-2)

saveRDS(phate.data,
        "~/phate.data")

# sqrt + snm
snmR.sqrt.cad <- readRDS("~/snmR_sqrt_cad.Rds")
colnames(snmR.sqrt.cad$norm.dat) <- colnames(snmR.sqrt.cad$raw.dat)

phate.snm.data <- phateR::phate(t(snmR.sqrt.cad$norm.dat),
                                n.jobs=-2)

saveRDS(phate.snm.data,
        "~/phate.snm.data.Rds")

# sqrt
saveRDS(sqrt.combined.mtx,
        "~/sqrt_combined_mtx.Rds")

phate.sqrt.data <- phateR::phate(t(sqrt.combined.mtx),
                                 n.jobs = -2)

saveRDS(phate.sqrt.data,
        "~/phate_sqrt_data.Rds")

palette(rainbow(10))
plot(phate.sqrt.data,
     col=as.numeric(as.factor(tissue.vec))^2)
plot(phate.sqrt.data,
     col=as.numeric(as.factor(datasource.vec))^2)

# vst 
vst.mtx.counts <- as.matrix(assay(vst.counts))

phate.vst.data <- phateR::phate(t(vst.mtx.counts),
                                 n.jobs = -2)

saveRDS(phate.vst.data,
        "~/phate_vst_data.Rds")

palette(rainbow(10))
plot(phate.vst.data,
     col=as.numeric(as.factor(tissue.vec))^2)
plot(phate.vst.data,
     col=as.numeric(as.factor(datasource.vec))^2)

# combined.df
com.mtx.counts <- as.matrix(combined.df)

phate.com.data <- phateR::phate(t(com.mtx.counts),
                                 n.jobs = -2)

saveRDS(phate.com.data,
        "~/phate_com_data.Rds")

palette(rainbow(10))
plot(phate.com.data,
     col=as.numeric(as.factor(tissue.vec))^2)
plot(phate.com.data,
     col=as.numeric(as.factor(datasource.vec))^2)

snmR.cad <- snm(vst.mtx.counts,
               bio.var=data.frame(tissue=tissue.vec),
               adj.var=data.frame(datasource=datasource.vec),
               rm.adj=TRUE,
               verbose=TRUE,
               diagnose=TRUE,
               nbins=5,
               lmer.max.iter=10,
               spline.dim=3,
               num.iter=3)

saveRDS(snmR.cad,
        "~/snmR_cad.Rds")

snm_eset <- snmR.cad$norm.dat %>% Biobase::ExpressionSet()

saveRDS(snm_eset,
        "~/snm_eset.Rds")
