library(dplyr)
library(magrittr)

IN_DIR <- "/home/jgburk/PycharmProjects/reticula/data/SRP035988/output/"
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

#generate pathway mapping for edges
REACTION_TO_PTHWY_FN <- "/home/jgburk/PycharmProjects/reticula/data/SRP035988/input/ReactionToPathway_Rel_71_122820.txt" # on box.com

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
edge_name_list <- list()
for(i in 1:nrow(labelled_edge_weights.df)){
  p_rxn <- labelled_edge_weights.df$Preceeding_Reaction[i]
  f_rxn <- labelled_edge_weights.df$Following_Reaction[i]
  edge_name <- paste(p_rxn,"->",f_rxn)
  edge_name_list[[i]] <- edge_name
  pthwy_candidates <- reaction_2_pthwy.nls[[p_rxn]]
  for(pthwy_candidate in pthwy_candidates){
    if(f_rxn %in% pthwy_2_reaction.nls[[pthwy_candidate]]){
      if(is.null(pthwy_2_edge.nls[[pthwy_candidate]])){
        pthwy_2_edge.nls[[pthwy_candidate]] <- edge_name
      }else{
        pthwy_2_edge.nls[[pthwy_candidate]] <- c(pthwy_2_edge.nls[[pthwy_candidate]],edge_name)
      }
      if(is.null(edge_2_pthwy.nls[[edge_name]])){
        edge_2_pthwy.nls[[edge_name]] <- pthwy_candidate
      }else{
        edge_2_pthwy.nls[[edge_name]] <- c(edge_2_pthwy.nls[[edge_name]],pthwy_candidate)
      }
    }
  }
  print(paste("Processed ",i," of ",nrow(labelled_edge_weights.df)," edges...",sep=""))
}

saveRDS(pthwy_2_edge.nls,file=paste(IN_DIR,"pthwy_2_edge_nls.Rds",sep=""))
saveRDS(edge_2_pthwy.nls,file=paste(IN_DIR,"edge_2_pthwy_nls.Rds",sep=""))
saveRDS(edge_name_list,file=paste(IN_DIR,"edge_name_list.Rds",sep=""))

pthwy_2_edge.nls <- readRDS(file=paste(IN_DIR,"pthwy_2_edge_nls.Rds",sep=""))
edge_2_pthwy.nls <- readRDS(file=paste(IN_DIR,"edge_2_pthwy_nls.Rds",sep=""))
edge_name_list <- readRDS(file=paste(IN_DIR,"edge_name_list.Rds",sep=""))

#calculate edge pathway enrichment
edge_pathway_enrichment.nls <- list()
for(pthwy in names(pthwy_2_edge.nls)){
  for(i in 1:2){
    significant_edges <- intersect(names(edge_2_pthwy.nls),unlist(edge_name_list[which(tissuewise_ig_edge_wilcox_res.nls[[i]] < ALPHA)]))
    edges_in_pathway <- edge_2_pthwy.nls[[pthwy]]
    q = length(intersect(significant_edges,edges_in_pathway))
    print(paste("significant edges in pathway: ",q,sep=""))
    m = length(edges_in_pathway)
    print(paste("edges in pathway: ",m,sep=""))
    n = length(names(edge_2_pthwy.nls)) - length(edges_in_pathway)
    print(paste("edges not in pathway: ",n,sep=""))
    k = length(significant_edges)
    print(paste("significant edges: ",k,sep=""))
    p <- stats::phyper(q-1,m,n,k,lower.tail = FALSE)
    print(paste("p-value: ",p,sep=""))
    edge_pathway_enrichment.nls[[tissue_name(i - 1)]][[pthwy]] <- p
  }
}

saveRDS(edge_pathway_enrichment.nls,file=paste(IN_DIR,"edge_pathway_enrichment_nls.Rds",sep=""))

rxn2ensembls.nls <- readRDS(paste(IN_DIR,"rxn2ensembls_nls.Rds",sep=""))

TOP_N = 20
top_weighted_rxns_lesional_psoriatic_skin <- labelled_edge_weights.df %>%
  dplyr::slice_max(labelled_edge_weights.df$`IG_lesional psoriatic skin`,n=TOP_N) %>%
  dplyr::select(Preceeding_Reaction)

top_weighted_rxn_lesional_psoriatic_skin_ensembls <- sort(unique(unlist(rxn2ensembls.nls[unlist(top_weighted_rxns_lesional_psoriatic_skin)])))
write.table(top_weighted_rxn_lesional_psoriatic_skin_ensembls,
            file=paste(IN_DIR,"top_weighted_rxn_psoriatic_skin_ensembls.txt",sep=""),
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE,
            sep="\n")

top_weighted_rxns_normal_skin <- labelled_edge_weights.df %>%
  dplyr::slice_max(labelled_edge_weights.df$`IG_normal skin`,n=TOP_N) %>%
  dplyr::select(Preceeding_Reaction)

top_weighted_rxn_normal_skin_ensembls <- sort(unique(unlist(rxn2ensembls.nls[unlist(top_weighted_rxns_normal_skin)])))
write.table(top_weighted_rxn_normal_skin_ensembls,
            file=paste(IN_DIR,"top_weighted_rxn_normal_skin_ensembls.txt",sep=""),
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE,
            sep="\n")
