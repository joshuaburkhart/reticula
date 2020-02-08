library(snm)

combine.df <- readRds("~/combine_df.Rds")
tissue.vec <- readRds("~/tissue_vec.Rds")
datasource.vec <- readRds("~/datasource_vec.Rds")

snmR.cad = snm(combine.df,
               tissue.vec,
               datasource.vec,
               rm.adj=TRUE,
               num.iter=5)

snm_eset <- snmR.cad$norm.dat %>% ExpressionSet()

saveRds(snm_eset,
        "~/snm_eset.Rds")