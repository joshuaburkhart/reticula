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

write.csv(labelled_edge_weights.df,file = paste(IN_DIR,"labelled_edge_weights.csv"))

#calculate wilcox for each edge across tissues, relative to others

#ig weights
tissuewise_ig_edge_wilcox_res.nls <- list()
for(edge_idx in 1:nrow(ig_edge_idx.df)){
  for(tissue_idx in 1:51){
    wilcox_res <- wilcox.test(x=as.numeric(ig_edge_idx.df[edge_idx,tissue_idx]),
                              y=as.numeric(ig_edge_idx.df[edge_idx,-tissue_idx]),
                              alternative = "greater")
    tissuewise_ig_edge_wilcox_res.nls[[tissue_name(tissue_idx - 1)]][[edge_idx]] <- wilcox_res$p.value
  }
  print(paste("Calculatd wilcox p-values for ",edge_idx," of ",nrow(ig_edge_idx.df)," ig edges...",sep=""))
}

#saliency weights
tissuewise_saliency_edge_wilcox_res.nls <- list()
for(edge_idx in 1:nrow(saliency_edge_idx.df)){
  for(tissue_idx in 1:51){
    wilcox_res <- wilcox.test(x=as.numeric(saliency_edge_idx.df[edge_idx,tissue_idx]),
                              y=as.numeric(saliency_edge_idx.df[edge_idx,-tissue_idx]),
                              alternative = "greater")
    tissuewise_saliency_edge_wilcox_res.nls[[tissue_name(tissue_idx - 1)]][[edge_idx]] <- wilcox_res$p.value
  }
  print(paste("Calculatd wilcox p-values for ",edge_idx," of ",nrow(ig_edge_idx.df)," saliency edges...",sep=""))
}

#generate pathway mapping for edges
REACTION_TO_PTHWY_FN <- "/home/burkhart/Software/reticula/data/aim1/input/ReactionToPathway_Rel_71_122820.csv" # on box.com

reaction_2_pthwy.df <- read.table(file=REACTION_TO_PTHWY_FN,sep="\t",header = TRUE,
                                  colClasses = c("character","character"))

pthwy_2_reaction.nls <- list()
reaction_2_pthwy.nls <- list()
for(i in 1:nrow(reaction_2_pthwy.df)){
  rxn <- reaction_2_pthwy.df$ReactionlikeEvent[i]
  pthwy <- reaction_2_pthwy.df$Pathway[i]
  if(is.null(pthwy_2_reaction.nls[[pthwy]])){
    pthwy_2_reaction.nls[[pthwy]] <- rxn
  }else{
    pthwy_2_reaction.nls[[pthwy]] <- c(pthwy_2_reaction.nls[[pthwy]],rxn)
  }
  if(is.null(reaction_2_pthwy.nls[[rxn]])){
    reaction_2_pthwy.nls[[rxn]] <- pthwy
  }else{
    reaction_2_pthwy.nls[[rxn]] <- c(reaction_2_pthwy.nls[[rxn]],pthwy)
  }
}

pthwy_2_edge.nls <- list()
edge_2_pthwy.nls <- list()
for(i in 1:nrow(labelled_edge_weights.df)){
  p_rxn <- labelled_edge_weights.df$Preceeding_Reaction[i]
  f_rxn <- labelled_edge_weights.df$Following_Reaction[i]
  pthwy_candidates <- reaction_2_pthwy.nls[[p_rxn]]
  for(pthwy_candidate in pthwy_candidates){
    if(f_rxn %in% pthwy_2_reaction.nls[[pthwy_candidate]]){
      edge <- paste(p_rxn,"->",f_rxn)
      if(is.null(pthwy_2_edge.nls[[pthwy_candidate]])){
        pthwy_2_edge.nls[[pthwy_candidate]] <- edge
      }else{
        pthwy_2_edge.nls[[pthwy_candidate]] <- c(pthwy_2_edge.nls[[pthwy_candidate]],edge)
      }
      if(is.null(edge_2_pthwy.nls[[edge]])){
        edge_2_pthwy.nls[[edge]] <- pthwy_candidate
      }else{
        edge_2_pthwy.nls[[edge]] <- c(edge_2_pthwy.nls[[edge]],pthwy_candidate)
      }
    }
  }
  print(paste("Processed ",i," of ",nrow(labelled_edge_weights.df)," edges...",sep=""))
}

#calculate edge pathway enrichment
for(pthwy in shared_pathways){
  #print(pthwy) #debugging
  reactions_in_pathway <- reaction_pathway_list[[pthwy]]
  q = length(intersect(significant_reactions.df$rxn_n1,reactions_in_pathway))
  m = length(reactions_in_pathway)
  n = length(unique(reaction_2_pthwy.df.shared$ReactionlikeEvent)) - length(reactions_in_pathway)
  k = nrow(significant_reactions.df)
  reaction_pathway_enrichment[[pthwy]] <- stats::phyper(q-1,m,n,k,lower.tail = FALSE)
}
