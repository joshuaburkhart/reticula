library(magrittr)
OUT_DIR <- "/home/jgburk/PycharmProjects/reticula/data/SRP035988/output/"
rxn2ensembls.nls <- readRDS(paste(OUT_DIR,"rxn2ensembls_nls.Rds",sep=""))
rxn_ensembl_counts.df <- rxn2ensembls.nls %>% lapply(.,length) %>% as.data.frame() %>% t() %>% as.data.frame()
rxn_ensembl_counts.df$V1 %>% hist(breaks=max(rxn_ensembl_counts.df$V1),xlim=range(0:25),labels=TRUE)

tot_represented_rxn_count <- rxn_ensembl_counts.df$V1 %>% length() #10516
gte4_represented_rxn_count <- rxn_ensembl_counts.df[rxn_ensembl_counts.df$V1 >= 4,] %>% length() #3449
gte8_represented_rxn_count <- rxn_ensembl_counts.df[rxn_ensembl_counts.df$V1 >= 8,] %>% length() #1612

gte4_represented_rxn_count/tot_represented_rxn_count * 100 #32.79764
gte8_represented_rxn_count/tot_represented_rxn_count * 100 #15.32902

rxn_ensembl_counts.nls <- rxn2ensembls.nls %>% lapply(.,length)
saveRDS(rxn_ensembl_counts.nls,
        paste(OUT_DIR,"rxn_ensembl_counts_nls.Rds",sep=""))
