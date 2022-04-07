library(DESeq2)
library(magrittr)
library(SummarizedExperiment)

start_time <- Sys.time()

IN_DIR <- "/home/jgburk/PycharmProjects/reticula/data/SRP035988/input/"
OUT_DIR <- "/home/jgburk/PycharmProjects/reticula/data/SRP035988/output/"

dds <- readRDS(paste(OUT_DIR,"dds.Rds",sep=""))

vst.counts <- DESeq2::vst(dds,
                          blind = FALSE,
                          fitType = "local")

saveRDS(vst.counts,
        paste(OUT_DIR,"vst_counts.Rds",sep=""))

ensembl2rxns.df <- read.table(paste(IN_DIR,"Ensembl2ReactomeReactions.txt",sep=""),
                              sep="\t",
                              comment.char = "")

reactome_ensembl_ids <- readRDS(paste(OUT_DIR,"reactome_ensembl_ids.Rds",sep=""))

rxn2ensembls.nls <- list()
rxns_w_srp035988_ensembls.df <- ensembl2rxns.df %>% dplyr::filter(V1 %in% reactome_ensembl_ids)
rxns_w_srp035988_ensembls.df$V1 <- as.character(rxns_w_srp035988_ensembls.df$V1)
rxns_w_srp035988_ensembls.df$V2 <- as.character(rxns_w_srp035988_ensembls.df$V2)
for(i in 1:nrow(rxns_w_srp035988_ensembls.df)){
  ens_id <- rxns_w_srp035988_ensembls.df$V1[i]
  rxn_id <- rxns_w_srp035988_ensembls.df$V2[i]
  ensembl_list_for_rxn_id <- rxn2ensembls.nls[[rxn_id]]
  if(is.null(ensembl_list_for_rxn_id)){
    ensembl_list_for_rxn_id <- c(ens_id)
  }
  rxn2ensembls.nls[[rxn_id]] <- c(ensembl_list_for_rxn_id,ens_id) %>% unique()
}
saveRDS(rxn2ensembls.nls,
        paste(OUT_DIR,"rxn2ensembls_nls.Rds",sep=""))

end_time <- Sys.time()
print(paste("Start: ",start_time," End: ",end_time," Difference: ",end_time - start_time))
