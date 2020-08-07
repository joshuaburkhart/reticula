set.seed(88888888) # maximum luck

library(magrittr)

start_time <- Sys.time()

OUT_DIR <- "/Users/burkhajo/Software/reticula/data/aim1/output/"#"/home/burkhart/Software/reticula/data/aim1/output/"

rxn2ensembls.nls <- readRDS(paste(OUT_DIR, "rxn2ensembls_nls.Rds", sep = ""))
rxn_knn_misclass_rate.nls <- readRDS(paste(OUT_DIR, "toi_rxn_knn_misclass_rate_nls.Rds", sep = ""))
rxn_knn_ari.nls <- readRDS(paste(OUT_DIR, "toi_rxn_knn_ari_nls.Rds", sep = ""))
rxn_knn_ecount.nls <- readRDS(paste(OUT_DIR, "toi_rxn_knn_ecount_nls.Rds", sep = ""))
rxn_pca.nls <- readRDS(paste(OUT_DIR, "rxn_pca_nls.Rds", sep = ""))

rxn_tissue_mean_misclass.df <- as.data.frame(do.call(rbind, rxn_knn_misclass_rate.nls))

# construct d
d <- data.frame(
  RXN_ID = names(rxn2ensembls.nls),
  MISCLASS = unlist(rxn_knn_misclass_rate.nls),
  ARI = unlist(rxn_knn_ari.nls),
  ECOUNT = unlist(rxn_knn_ecount.nls)
)

# store d
saveRDS(d, paste(OUT_DIR, "toi_summary_df.Rds", sep = ""))

# make figures using d



# clustering on rxn_pca.nls (review compare_toi_counts_to_pca.R for extracting contributions)

end_time <- Sys.time()
print(paste(
  "Start: ",
  start_time,
  " End: ",
  end_time,
  " Difference: ",
  end_time - start_time
))