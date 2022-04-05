library(magrittr)

IN_DIR <- "/home/jgburk/PycharmProjects/reticula/data/SRP050223/output/"
GTEX_DIR <- "/home/jgburk/PycharmProjects/reticula/data/gtex/"
ALPHA <- 0.05

tissue2idx.df <- data.frame(read.table(paste(IN_DIR,"inverted_targets.txt",sep=""),
                                       stringsAsFactors = FALSE),
                            read.table(paste(IN_DIR,"transformed_targets.txt",sep=""),
                                       stringsAsFactors = FALSE)) %>%
  unique()
colnames(tissue2idx.df) <- c("Tissue","Index")

tissue_name <- function(tissue_index){
  return(tissue2idx.df$Tissue[tissue2idx.df$Index == tissue_index])
}

node2idx.df <- data.frame(read.table(paste(GTEX_DIR,"Copy\ of\ nodeLabel2rxn_nls.csv",sep="")),
                          read.table(paste(GTEX_DIR,"Copy\ of\ rxn2nodeLabel_nls.csv",sep="")))
colnames(node2idx.df) <- c("Index","Reaction","Reaction_DoubleCheck","Index_DoubleCheck")

reaction_name <- function(reaction_index){
  r_name <- node2idx.df$Reaction[node2idx.df$Index == reaction_index]
  r_check <- node2idx.df$Reaction_DoubleCheck[node2idx.df$Index_DoubleCheck == reaction_index]
  if(r_name == r_check){
    return(r_name)
  }
  return(FALSE)
}

edge2idx.df <- data.frame(read.table(paste(GTEX_DIR,"Copy\ of\ edges.txt",sep="")),
                          read.table(paste(GTEX_DIR,"Copy\ of\ edgeLabels.csv",sep="")))
colnames(edge2idx.df) <- c("Preceeding_Index","Following_Index","Preceeding_Reaction","Following_Reaction")

ig_edge_idx.df <- data.frame(read.table(paste(IN_DIR, "ig_0.txt",sep="")))
ig_colnames <- c(paste("IG_",tissue_name(0),sep=""))
ig_edge_idx.df <- cbind(ig_edge_idx.df,read.table(paste(IN_DIR,"ig_",1,".txt",sep="")))
ig_colnames <- c(ig_colnames,paste("IG_",tissue_name(1),sep=""))
colnames(ig_edge_idx.df) <- ig_colnames

tissue2idx.df %>% dim()
node2idx.df %>% dim()
edge2idx.df %>% dim()
ig_edge_idx.df %>% dim()

labelled_edge_weights.df <- cbind(edge2idx.df,ig_edge_idx.df)

# checking reaction pairs match indices
for(i in 1:nrow(labelled_edge_weights.df)){
  p_idx <- labelled_edge_weights.df$Preceeding_Index[i]
  f_idx <- labelled_edge_weights.df$Following_Index[i]
  p_rxn <- labelled_edge_weights.df$Preceeding_Reaction[i]
  f_rxn <- labelled_edge_weights.df$Following_Reaction[i]
  if(p_rxn != reaction_name(p_idx) | f_rxn != reaction_name(f_idx)){
    print(paste("Index mismatch on line ",i,
                ": p_idx = ",p_idx,
                ". f_idx = ",f_idx,
                ". p_rxn = ",p_rxn,
                ". f_rxn = ",f_rxn,
                ". reaction_name(p_idx) = ", reaction_name(p_idx),
                ". reaction_name(f_idx) = ", reaction_name(f_idx)))
    break
  }else{
   print(paste("line ",i," seems fine...")) 
  }
}

write.csv(labelled_edge_weights.df,file = paste(IN_DIR,"labelled_edge_weights.csv",sep=""))

#calculate wilcox for each edge across tissues, relative to others

#ig weights
tissuewise_ig_edge_wilcox_res.nls <- list()
for(edge_idx in 1:nrow(ig_edge_idx.df)){
  for(tissue_idx in 1:2){
    wilcox_res <- wilcox.test(x=as.numeric(ig_edge_idx.df[edge_idx,tissue_idx]),
                              y=as.numeric(ig_edge_idx.df[edge_idx,-tissue_idx]),
                              alternative = "greater")
    tissuewise_ig_edge_wilcox_res.nls[[tissue_name(tissue_idx - 1)]][[edge_idx]] <- wilcox_res$p.value
  }
  print(paste("Calculated wilcox p-values for ",edge_idx," of ",nrow(ig_edge_idx.df)," ig edges...",sep=""))
}

saveRDS(tissuewise_ig_edge_wilcox_res.nls,file=paste(IN_DIR,"tissuewise_ig_edge_wilcox_res_nls.Rds",sep=""))

# note lowest possible wilcoxon p-value for an edge is given by "wilcox.test(x = c(51), y= seq(1:50),alternative = "greater")"