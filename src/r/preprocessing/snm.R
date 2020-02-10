library(snm)
library(limma)
library(phateR)
library(Biobase)
library(magrittr)
library(SummarizedExperiment)

vst.counts <- readRDS("~/vst_counts.Rds")
tissue.vec <- readRDS("~/tissue_vec.Rds")
datasource.vec <- readRDS("~/datasource_vec.Rds")

combined.df <- readRDS("~/combined_df.Rds")
sqrt.combined.mtx <- as.matrix(sqrt(combined.df))

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

phate.data <- phateR::phate(vst.rbe.t,
                            n.jobs=-2)

saveRDS(phate.data,
        "~/phate.data")

palette(rainbow(10))
plot(phate.data,
     col=as.numeric(as.factor(tissue.vec)))
plot(phate.data,
     col=as.numeric(as.factor(datasource.vec)))

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
