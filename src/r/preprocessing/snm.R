library(snm)
library(Biobase)
library(magrittr)
library(SummarizedExperiment)

vst.counts <- readRDS("~/vst_counts.Rds")
tissue.vec <- readRDS("~/tissue_vec.Rds")
datasource.vec <- readRDS("~/datasource_vec.Rds")

snmR.cad <- snm(as.matrix(assay(vst.counts)),
               bio.var=data.frame(tissue=tissue.vec),
               adj.var=data.frame(datasource=datasource.vec),
               rm.adj=TRUE,
               verbose=TRUE,
               num.iter=5)

saveRDS(snmR.cad,
        "~/snmR_cad.Rds")

snm_eset <- snmR.cad$norm.dat %>% Biobase::ExpressionSet()

saveRDS(snm_eset,
        "~/snm_eset.Rds")
