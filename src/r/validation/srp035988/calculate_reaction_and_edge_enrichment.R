# Review phyper function used in https://github.com/joshuaburkhart/reticula/blob/master/src/r/analysis/hypergeometric_enrichment_analysis.R
# use phyper to calculate (positive) enrichment with fisher exact p-value (explained here: https://stackoverflow.com/questions/53051977/p-value-from-fisher-test-does-not-match-phyper)
# #phyper(q - 1, m, n, k, lower.tail = FALSE)
# q = number of white balls drawn (number of transcripts/reactions shared between selection & pathway) 
# m = number of white balls in urn (number of transcripts/reactions in pathway)
# n = number of black balls in urn (number of transcripts/reactions not in pathway)
# k = number of balls drawn (number of transcripts/reactions in selection)

library(magrittr)
library(dplyr)

ETA <- 1e-5
N_RXNS <- 10516
N_EDGES <- 40032
OUT_DIR <- "/home/jgburk/PycharmProjects/reticula/data/SRP035988/output/"

# load rxn2ensembls.nls
rxn2ensembls.nls <- readRDS(paste(OUT_DIR,"rxn2ensembls_nls.Rds",sep=""))

# load res.df
res.df <- read.table(file=paste(OUT_DIR,"res_df.csv",sep=""), header = TRUE, stringsAsFactors = FALSE)

# represented reaction enrichment

get_reaction_enrichment_p_value <- function(x,y,z){
  # x = reaction name
  # y = significant transcripts
  # z = all transcripts
  q = length(intersect(y,rxn2ensembls.nls[[x]]))
  m = length(rxn2ensembls.nls[[x]])
  n = length(z) - length(rxn2ensembls.nls[[x]])
  k = length(y)
  return(c(paste(x),stats::phyper(q - 1,m,n,k,lower.tail = FALSE)))
}

# Default Thresholds
default_thresholds_reaction_enrichment_p_values <- mapply(get_reaction_enrichment_p_value,
                             names(rxn2ensembls.nls),
                             res.df %>%
                               dplyr::filter(Sig_W_Default_Thresholds == TRUE) %>%
                               dplyr::select(EnsemblID),
                             res.df %>%
                               dplyr::select(EnsemblID))

# Strict Thresholds
strict_thresholds_reaction_enrichment_p_values <- mapply(get_reaction_enrichment_p_value,
                                                          names(rxn2ensembls.nls),
                                                          res.df %>%
                                                            dplyr::filter(Sig_W_Strict_Thresholds == TRUE) %>%
                                                            dplyr::select(EnsemblID),
                                                          res.df %>%
                                                            dplyr::select(EnsemblID))

# Relaxed Thresholds
relaxed_thresholds_reaction_enrichment_p_values <- mapply(get_reaction_enrichment_p_value,
                                                          names(rxn2ensembls.nls),
                                                          res.df %>%
                                                            dplyr::filter(Sig_W_Relaxed_Thresholds == TRUE) %>%
                                                            dplyr::select(EnsemblID),
                                                          res.df %>%
                                                            dplyr::select(EnsemblID))

# check "R-HSA-8956184" = 0.01576931
spot_test <- relaxed_thresholds_reaction_enrichment_p_values %>% t() %>% as.data.frame() %>% .[.$V1=="R-HSA-8956184",]
stopifnot(assertthat::are_equal(abs(as.numeric(spot_test$V2) - 0.01576931) < ETA,TRUE))

stopifnot(assertthat::are_equal(ncol(strict_thresholds_reaction_enrichment_p_values),N_RXNS))

stopifnot(assertthat::are_equal(dim(strict_thresholds_reaction_enrichment_p_values),
                                dim(default_thresholds_reaction_enrichment_p_values)))
stopifnot(assertthat::are_equal(dim(strict_thresholds_reaction_enrichment_p_values),
                                dim(relaxed_thresholds_reaction_enrichment_p_values)))

reaction_enrichment.df <- data.frame(row.names=strict_thresholds_reaction_enrichment_p_values %>% t() %>% as.data.frame() %>% .$V1,
                                     strict_thresholds_p_value=strict_thresholds_reaction_enrichment_p_values %>% t() %>% as.data.frame() %>% .$V2,
                                     default_thresholds_p_value=default_thresholds_reaction_enrichment_p_values %>% t() %>% as.data.frame() %>% .$V2,
                                     relaxed_thresholds_p_value=relaxed_thresholds_reaction_enrichment_p_values %>% t() %>% as.data.frame() %>% .$V2)

reaction_enrichment.df$strict_thresholds_BH_adj_p_value = p.adjust(reaction_enrichment.df$strict_thresholds_p_value, method = "BH")
reaction_enrichment.df$default_thresholds_BH_adj_p_value = p.adjust(reaction_enrichment.df$default_thresholds_p_value, method = "BH")
reaction_enrichment.df$relaxed_thresholds_BH_adj_p_value = p.adjust(reaction_enrichment.df$relaxed_thresholds_p_value, method = "BH")

write.table(reaction_enrichment.df,
            file=paste(OUT_DIR,"reaction_enrichment_df.csv",sep=""))

# edge enrichment

# load labelled_edge_weights.df
labelled_edge_weights.df <- read.table(paste(OUT_DIR,"labelled_edge_weights.csv",sep=""),
                                       stringsAsFactors = FALSE,sep = ",",header = TRUE)

get_edge_enrichment_p_value <- function(w,x,y,z){
  # w = preceding reaction name
  # x = following reaction name
  # y = significant transcripts
  # z = all transcripts
  edge_rxns = unique(c(rxn2ensembls.nls[[w]],
                       rxn2ensembls.nls[[x]]))
  q = length(intersect(y,edge_rxns))
  m = length(edge_rxns)
  n = length(z) - length(edge_rxns)
  k = length(y)
  return(c(paste(w,x,sep=" -> "),stats::phyper(q - 1,m,n,k,lower.tail = FALSE)))
}

# Default Thresholds
default_thresholds_edge_enrichment_p_values <- mapply(get_edge_enrichment_p_value,
                             labelled_edge_weights.df$Preceeding_Reaction,
                             labelled_edge_weights.df$Following_Reaction,
                             res.df %>%
                               dplyr::filter(Sig_W_Default_Thresholds == TRUE) %>%
                               dplyr::select(EnsemblID),
                             res.df %>%
                               dplyr::select(EnsemblID))

# Strict Thresholds
strict_thresholds_edge_enrichment_p_values <- mapply(get_edge_enrichment_p_value,
                                                      labelled_edge_weights.df$Preceeding_Reaction,
                                                      labelled_edge_weights.df$Following_Reaction,
                                                      res.df %>%
                                                        dplyr::filter(Sig_W_Strict_Thresholds == TRUE) %>%
                                                        dplyr::select(EnsemblID),
                                                      res.df %>%
                                                        dplyr::select(EnsemblID))

# Relaxed Thresholds
relaxed_thresholds_edge_enrichment_p_values <- mapply(get_edge_enrichment_p_value,
                                                      labelled_edge_weights.df$Preceeding_Reaction,
                                                      labelled_edge_weights.df$Following_Reaction,
                                                      res.df %>%
                                                        dplyr::filter(Sig_W_Relaxed_Thresholds == TRUE) %>%
                                                        dplyr::select(EnsemblID),
                                                      res.df %>%
                                                        dplyr::select(EnsemblID))

stopifnot(assertthat::are_equal(ncol(strict_thresholds_edge_enrichment_p_values),N_EDGES))

stopifnot(assertthat::are_equal(dim(strict_thresholds_edge_enrichment_p_values),
                                dim(default_thresholds_edge_enrichment_p_values)))
stopifnot(assertthat::are_equal(dim(strict_thresholds_edge_enrichment_p_values),
                                dim(relaxed_thresholds_edge_enrichment_p_values)))

edge_enrichment.df <- data.frame(row.names=strict_thresholds_edge_enrichment_p_values %>% t() %>% as.data.frame() %>% .$V1,
                                     strict_thresholds_p_value=strict_thresholds_edge_enrichment_p_values %>% t() %>% as.data.frame() %>% .$V2,
                                     default_thresholds_p_value=default_thresholds_edge_enrichment_p_values %>% t() %>% as.data.frame() %>% .$V2,
                                     relaxed_thresholds_p_value=relaxed_thresholds_edge_enrichment_p_values %>% t() %>% as.data.frame() %>% .$V2)

edge_enrichment.df$strict_thresholds_BH_adj_p_value = p.adjust(edge_enrichment.df$strict_thresholds_p_value, method = "BH")
edge_enrichment.df$default_thresholds_BH_adj_p_value = p.adjust(edge_enrichment.df$default_thresholds_p_value, method = "BH")
edge_enrichment.df$relaxed_thresholds_BH_adj_p_value = p.adjust(edge_enrichment.df$relaxed_thresholds_p_value, method = "BH")

write.table(edge_enrichment.df,
            file=paste(OUT_DIR,"edge_enrichment_df.csv",sep=""))
