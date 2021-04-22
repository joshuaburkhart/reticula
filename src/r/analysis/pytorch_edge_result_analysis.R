IN_DIR <- "/home/burkhart/Software/reticula/data/aim2/input/"

tissue2idx.df <- data.frame(read.table(paste(IN_DIR,"inverted_targets.txt",sep=""),
                                       stringsAsFactors = FALSE),
                            read.table(paste(IN_DIR,"transformed_targets.txt",sep=""),
                                       stringsAsFactors = FALSE)) %>%
  unique()
colnames(tissue2idx.df) <- c("Tissue","Index")

tissue_name <- function(tissue_index){
  return(tissue2idx.df$Tissue[tissue2idx.df$Index == tissue_index])
}

node2idx.df <- data.frame(read.table(paste(IN_DIR,"nodeLabel2rxn_nls.csv",sep="")),
                          read.table(paste(IN_DIR,"rxn2nodeLabel_nls.csv",sep="")))
colnames(node2idx.df) <- c("Index","Reaction","Reaction_DoubleCheck","Index_DoubleCheck")

reaction_name <- function(reaction_index){
  r_name <- node2idx.df$Reaction[node2idx.df$Index == reaction_index]
  r_check <- node2idx.df$Reaction_DoubleCheck[node2idx.df$Index_DoubleCheck == reaction_index]
  if(r_name == r_check){
    return(r_name)
  }
  return(FALSE)
}

edge2idx.df <- data.frame(read.table(paste(IN_DIR,"edges.txt",sep="")),
                          read.table(paste(IN_DIR,"edgeLabels.csv",sep="")))
colnames(edge2idx.df) <- c("Preceeding_Index","Following_Index","Preceeding_Reaction","Following_Reaction")

ig_edge_idx.df <- data.frame(read.table(paste(IN_DIR, "ig_0.txt",sep="")))
ig_colnames <- c(paste("IG_",tissue_name(0),sep=""))
for(i in 1:50) {
  ig_edge_idx.df <- cbind(ig_edge_idx.df,read.table(paste(IN_DIR,"ig_",i,".txt",sep="")))
  ig_colnames <- c(ig_colnames,paste("IG_",tissue_name(i),sep=""))
}
colnames(ig_edge_idx.df) <- ig_colnames

saliency_edge_idx.df <- data.frame(read.table(paste(IN_DIR, "saliency_0.txt",sep="")))
saliency_colnames <- c(paste("Saliency_",tissue_name(0),sep=""))
for(i in 1:50) {
  saliency_edge_idx.df <- cbind(saliency_edge_idx.df,read.table(paste(IN_DIR,"saliency_",i,".txt",sep="")))
  saliency_colnames <- c(saliency_colnames,paste("Saliency_",tissue_name(i),sep=""))
}
colnames(saliency_edge_idx.df) <- saliency_colnames

tissue2idx.df %>% dim()
node2idx.df %>% dim()
edge2idx.df %>% dim()
ig_edge_idx.df %>% dim()
saliency_edge_idx.df %>% dim()

labelled_edge_weights.df <- cbind(edge2idx.df,ig_edge_idx.df,saliency_edge_idx.df)

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
