library(snm)
library(Biobase)
library(magrittr)

combined.df <- readRDS("~/combined_df.Rds")
tissue.vec <- readRDS("~/tissue_vec.Rds")
datasource.vec <- readRDS("~/datasource_vec.Rds")

snmR.cad = snm(as.matrix(combined.df),
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
