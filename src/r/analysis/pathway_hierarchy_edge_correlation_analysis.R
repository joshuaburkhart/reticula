library(dplyr)

# breast, lung

IN_DIR <- "/home/burkhart/Software/reticula/data/aim2/input/"
OUT_DIR <- "/home/burkhart/Software/reticula/data/aim1/output/"

labelled_edge_weights.df <- read.table(file=paste(IN_DIR,"pathway_hierarchy_labelled_edge_weights.csv",sep=""),header = TRUE,sep = ",")
misclass_rates.df <- readRDS(file=paste(OUT_DIR,"toi_summary_df.Rds",sep=""))

TIS_NAME <- "Breast - Mammary Tissue"
#TIS_NAME <- "Lung"

ig_TIS_NAME <- "IG_Breast...Mammary.Tissue"
#ig_TIS_NAME <- "IG_Lung"
sal_TIS_NAME <- "Saliency_Breast...Mammary.Tissue"
#sal_TIS_NAME <- "Saliency_Lung"

tis_edges.df <- labelled_edge_weights.df %>% dplyr::select(Preceeding_Reaction,
                                                              Following_Reaction,
                                                              ig_TIS_NAME,
                                                              sal_TIS_NAME)
top_ig <- list()
top_saliency <- list()
for(i in 1:nrow(tis_edges.df)){
  p_rxn <- as.character(tis_edges.df[i,"Preceeding_Reaction"])
  f_rxn <- as.character(tis_edges.df[i,"Following_Reaction"])
  edge_ig <- tis_edges.df[i,ig_TIS_NAME]
  edge_saliency <- tis_edges.df[i,sal_TIS_NAME]
  if(is.null(top_ig[[p_rxn]])){
    top_ig[[p_rxn]] <- edge_ig
  }else if(edge_ig > top_ig[[p_rxn]]){
    top_ig[[p_rxn]] <- edge_ig
  }
  if(is.null(top_ig[[f_rxn]])){
    top_ig[[f_rxn]] <- edge_ig
  }else if(edge_ig > top_ig[[f_rxn]]){
    top_ig[[f_rxn]] <- edge_ig
  }
  if(is.null(top_saliency[[p_rxn]])){
    top_saliency[[p_rxn]] <- edge_saliency
  }else if(edge_saliency > top_saliency[[p_rxn]]){
    top_saliency[[p_rxn]] <- edge_saliency
  }
  if(is.null(top_saliency[[f_rxn]])){
    top_saliency[[f_rxn]] <- edge_saliency
  }else if(edge_saliency > top_saliency[[f_rxn]]){
    top_saliency[[f_rxn]] <- edge_saliency
  }
  print(paste("Processed row ",i," of ",nrow(tis_edges.df),"...",sep=""))
}

tis_preceeding_reaction_acc.vec <- numeric()
tis_following_reaction_acc.vec <- numeric()
preceeding_ari.vec <- numeric()
following_ari.vec <- numeric()
preceeding_top_ig.vec <- numeric()
following_top_ig.vec <- numeric()
preceeding_top_saliency.vec <- numeric()
following_top_saliency.vec <- numeric()
for(i in 1:nrow(tis_edges.df)){
  p_rxn <- as.character(tis_edges.df[i,"Preceeding_Reaction"])
  f_rxn <- as.character(tis_edges.df[i,"Following_Reaction"])
  tis_preceeding_reaction_acc.vec <- c(tis_preceeding_reaction_acc.vec,1 - misclass_rates.df[p_rxn,TIS_NAME])
  tis_following_reaction_acc.vec <- c(tis_following_reaction_acc.vec,1 - misclass_rates.df[f_rxn,TIS_NAME])
  preceeding_ari.vec <- c(preceeding_ari.vec,misclass_rates.df[p_rxn,"ARI"])
  following_ari.vec <- c(following_ari.vec,misclass_rates.df[f_rxn,"ARI"])
  preceeding_top_ig.vec <- c(preceeding_top_ig.vec,top_ig[[p_rxn]])
  following_top_ig.vec <- c(following_top_ig.vec,top_ig[[f_rxn]])
  preceeding_top_saliency.vec <- c(preceeding_top_saliency.vec,top_saliency[[p_rxn]])
  following_top_saliency.vec <- c(following_top_saliency.vec,top_saliency[[f_rxn]])
  print(paste("Processed row ",i," of ",nrow(tis_edges.df),"...",sep=""))
}

length(tis_preceeding_reaction_acc.vec)
length(tis_following_reaction_acc.vec)
length(preceeding_ari.vec)
length(following_ari.vec)
length(preceeding_top_ig.vec)
length(following_top_ig.vec)
length(preceeding_top_saliency.vec)
length(following_top_saliency.vec)

tis_edges.df$preceeding_acc <- tis_preceeding_reaction_acc.vec
tis_edges.df$following_acc <- tis_following_reaction_acc.vec
tis_edges.df$preceeding_ari <- preceeding_ari.vec
tis_edges.df$following_ari <- following_ari.vec
tis_edges.df$acc_diff <- 1 - abs(tis_preceeding_reaction_acc.vec - tis_following_reaction_acc.vec)
tis_edges.df$ari_diff <- 1 - abs(preceeding_ari.vec - following_ari.vec)
tis_edges.df$preceeding_top_ig <- preceeding_top_ig.vec
tis_edges.df$following_top_ig <- following_top_ig.vec
tis_edges.df$preceeding_top_saliency <- preceeding_top_saliency.vec
tis_edges.df$following_top_saliency <- following_top_saliency.vec

x_vector <- tis_edges.df$preceeding_top_ig
y_vector <- tis_edges.df$preceeding_acc

plot(x=x_vector,y=y_vector)
abline(lm(y_vector~x_vector),col="red")
lines(lowess(x_vector,y_vector),col="blue")
cor.test(x_vector,y_vector)

x_vector <- tis_edges.df$preceeding_top_saliency
y_vector <- tis_edges.df$preceeding_acc

plot(x=x_vector,y=y_vector)
abline(lm(y_vector~x_vector),col="red")
lines(lowess(x_vector,y_vector),col="blue")
cor.test(x_vector,y_vector)

x_vector <- tis_edges.df$preceeding_top_ig
y_vector <- tis_edges.df$preceeding_ari

plot(x=x_vector,y=y_vector)
abline(lm(y_vector~x_vector),col="red")
lines(lowess(x_vector,y_vector),col="blue")
cor.test(x_vector,y_vector)

x_vector <- tis_edges.df$preceeding_top_saliency
y_vector <- tis_edges.df$preceeding_ari

plot(x=x_vector,y=y_vector)
abline(lm(y_vector~x_vector),col="red")
lines(lowess(x_vector,y_vector),col="blue")
cor.test(x_vector,y_vector)

###
### lung results below
###

# preceeding_top_ig vs preceeding_acc


# preceeding_top_saliency vs preceeding_acc


# preceeding_top_ig vs preceeding_ari


# preceeding_top_saliency vs preceeding_ari


###
### breast results below
###

# preceeding_top_ig vs preceeding_acc


# preceeding_top_saliency vs preceeding_acc


# preceeding_top_ig vs preceeding_ari


# preceeding_top_saliency vs preceeding_ari

