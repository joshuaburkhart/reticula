library(dplyr)

# breast, lung

IN_DIR <- "/home/burkhart/Software/reticula/data/aim2/input/"
OUT_DIR <- "/home/burkhart/Software/reticula/data/aim1/output/"

labelled_edge_weights.df <- read.table(file=paste(IN_DIR,"labelled_edge_weights.csv",sep=""),header = TRUE,sep = ",")
misclass_rates.df <- readRDS(file=paste(OUT_DIR,"toi_summary_df.Rds",sep=""))


lung_edges.df <- labelled_edge_weights.df %>% dplyr::select(Preceeding_Reaction,
                                                              Following_Reaction,
                                                              IG_Lung,
                                                              Saliency_Lung)
top_ig <- list()
top_saliency <- list()
for(i in 1:nrow(lung_edges.df)){
  p_rxn <- as.character(lung_edges.df[i,"Preceeding_Reaction"])
  f_rxn <- as.character(lung_edges.df[i,"Following_Reaction"])
  edge_ig <- lung_edges.df[i,"IG_Lung"]
  edge_saliency <- lung_edges.df[i,"Saliency_Lung"]
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
  print(paste("Processed row ",i," of ",nrow(lung_edges.df),"...",sep=""))
}

lung_preceeding_reaction_acc.vec <- numeric()
lung_following_reaction_acc.vec <- numeric()
preceeding_ari.vec <- numeric()
following_ari.vec <- numeric()
preceeding_top_ig.vec <- numeric()
following_top_ig.vec <- numeric()
preceeding_top_saliency.vec <- numeric()
following_top_saliency.vec <- numeric()
for(i in 1:nrow(lung_edges.df)){
  p_rxn <- as.character(lung_edges.df[i,"Preceeding_Reaction"])
  f_rxn <- as.character(lung_edges.df[i,"Following_Reaction"])
  lung_preceeding_reaction_acc.vec <- c(lung_preceeding_reaction_acc.vec,1 - misclass_rates.df[p_rxn,"Lung"])
  lung_following_reaction_acc.vec <- c(lung_following_reaction_acc.vec,1 - misclass_rates.df[f_rxn,"Lung"])
  preceeding_ari.vec <- c(preceeding_ari.vec,misclass_rates.df[p_rxn,"ARI"])
  following_ari.vec <- c(following_ari.vec,misclass_rates.df[f_rxn,"ARI"])
  preceeding_top_ig.vec <- c(preceeding_top_ig.vec,top_ig[[p_rxn]])
  following_top_ig.vec <- c(following_top_ig.vec,top_ig[[f_rxn]])
  preceeding_top_saliency.vec <- c(preceeding_top_saliency.vec,top_saliency[[p_rxn]])
  following_top_saliency.vec <- c(following_top_saliency.vec,top_saliency[[f_rxn]])
  print(paste("Processed row ",i," of ",nrow(lung_edges.df),"...",sep=""))
}

length(lung_preceeding_reaction_acc.vec)
length(lung_following_reaction_acc.vec)
length(preceeding_ari.vec)
length(following_ari.vec)
length(preceeding_top_ig.vec)
length(following_top_ig.vec)
length(preceeding_top_saliency.vec)
length(following_top_saliency.vec)

lung_edges.df$preceeding_acc <- lung_preceeding_reaction_acc.vec
lung_edges.df$following_acc <- lung_following_reaction_acc.vec
lung_edges.df$preceeding_ari <- preceeding_ari.vec
lung_edges.df$following_ari <- following_ari.vec
lung_edges.df$acc_diff <- 1 - abs(lung_preceeding_reaction_acc.vec - lung_following_reaction_acc.vec)
lung_edges.df$ari_diff <- 1 - abs(preceeding_ari.vec - following_ari.vec)
lung_edges.df$preceeding_top_ig <- preceeding_top_ig.vec
lung_edges.df$following_top_ig <- following_top_ig.vec
lung_edges.df$preceeding_top_saliency <- preceeding_top_saliency.vec
lung_edges.df$following_top_saliency <- following_top_saliency.vec

x_vector <- lung_edges.df$preceeding_top_ig
y_vector <- lung_edges.df$preceeding_acc

plot(x=x_vector,y=y_vector)
abline(lm(y_vector~x_vector),col="red")
lines(lowess(x_vector,y_vector),col="blue")
cor.test(x_vector,y_vector)

x_vector <- lung_edges.df$preceeding_top_saliency
y_vector <- lung_edges.df$preceeding_acc

plot(x=x_vector,y=y_vector)
abline(lm(y_vector~x_vector),col="red")
lines(lowess(x_vector,y_vector),col="blue")
cor.test(x_vector,y_vector)

x_vector <- lung_edges.df$preceeding_top_ig
y_vector <- lung_edges.df$preceeding_ari

plot(x=x_vector,y=y_vector)
abline(lm(y_vector~x_vector),col="red")
lines(lowess(x_vector,y_vector),col="blue")
cor.test(x_vector,y_vector)

x_vector <- lung_edges.df$preceeding_top_saliency
y_vector <- lung_edges.df$preceeding_ari

plot(x=x_vector,y=y_vector)
abline(lm(y_vector~x_vector),col="red")
lines(lowess(x_vector,y_vector),col="blue")
cor.test(x_vector,y_vector)

###
### lung results below
###

# preceeding_top_ig vs preceeding_acc
# Pearson's product-moment correlation
# 
# data:  x_vector and y_vector
# t = 85.855, df = 40030, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.3860361 0.4025817
# sample estimates:
#       cor 
# 0.3943409

# preceeding_top_saliency vs preceeding_acc
# Pearson's product-moment correlation
# 
# data:  x_vector and y_vector
# t = 90.479, df = 40030, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.4038843 0.4201501
# sample estimates:
#     cor 
# 0.41205

# preceeding_top_ig vs preceeding_ari
# Pearson's product-moment correlation
# 
# data:  x_vector and y_vector
# t = 91.999, df = 40030, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.4096537 0.4258264
# sample estimates:
#       cor 
# 0.4177731

# preceeding_top_saliency vs preceeding_ari
# Pearson's product-moment correlation
# 
# data:  x_vector and y_vector
# t = 93.85, df = 40030, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.4166125 0.4326713
# sample estimates:
#       cor 
# 0.4246753

###
### breast results below
###

# preceeding_top_ig vs preceeding_acc
# Pearson's product-moment correlation
# 
# data:  x_vector and y_vector
# t = 110.89, df = 40030, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.4772404 0.4922285
# sample estimates:
#     cor 
# 0.48477 

# preceeding_top_saliency vs preceeding_acc
# Pearson's product-moment correlation
# 
# data:  x_vector and y_vector
# t = 88.785, df = 40030, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.3973961 0.4137650
# sample estimates:
#      cor 
# 0.405613

# preceeding_top_ig vs preceeding_ari
# Pearson's product-moment correlation
# 
# data:  x_vector and y_vector
# t = 113.47, df = 40030, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.4858665 0.5006909
# sample estimates:
#       cor 
# 0.4933146 

# preceeding_top_saliency vs preceeding_ari
# Pearson's product-moment correlation
# 
# data:  x_vector and y_vector
# t = 86.74, df = 40030, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.3894862 0.4059786
# sample estimates:
#       cor 
# 0.3977645
