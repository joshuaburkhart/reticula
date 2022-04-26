library(magrittr)
OUT_DIR <- "/home/jgburk/PycharmProjects/reticula/data/SRP035988/output/"
labelled_edge_weights.df <- read.table(paste(OUT_DIR,"labelled_edge_weights.csv",sep=""),
                                       stringsAsFactors = FALSE,sep = ",",header = TRUE)
rxn_ensembl_counts.nls <- readRDS(paste(OUT_DIR,"rxn_ensembl_counts_nls.Rds",sep=""))
get_transcript_count <- function(x){
  return(c(x,rxn_ensembl_counts.nls[[x]]))
}
p_rxn_transcript_count <- labelled_edge_weights.df$Preceeding_Reaction %>% lapply(.,get_transcript_count)
p_c.df <- p_rxn_transcript_count %>% as.data.frame() %>% t() %>% as.data.frame()
labelled_edge_weights.df$Preceding_T_Check <- p_c.df$V1
labelled_edge_weights.df$Preceding_T_Count <- p_c.df$V2

f_rxn_transcript_count <- labelled_edge_weights.df$Following_Reaction %>% lapply(.,get_transcript_count)
f_c.df <- f_rxn_transcript_count %>% as.data.frame() %>% t() %>% as.data.frame()
labelled_edge_weights.df$Following_T_Check <- f_c.df$V1
labelled_edge_weights.df$Following_T_Count <- f_c.df$V2

write.table(labelled_edge_weights.df,
            paste(OUT_DIR,"labelled_and_counted_edge_weights.csv",sep=""),
            row.names = FALSE,sep = ",")
