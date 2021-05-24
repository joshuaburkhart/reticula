library(dplyr)

# breast, lung

IN_DIR <- "/home/burkhart/Software/reticula/data/aim2/input/"
OUT_DIR <- "/home/burkhart/Software/reticula/data/aim1/output/"

labelled_edge_weights.df <- read.table(file=paste(IN_DIR,"pathway_hierarchy_labelled_edge_weights.csv",sep=""),header = TRUE,sep = ",")
misclass_rates.df <- readRDS(file=paste(OUT_DIR,"toi_summary_df.Rds",sep=""))
vst.count.mtx.train <- readRDS(paste(OUT_DIR,"vst_count_mtx_train.Rds",sep=""))
gtex_tissue_detail.vec.train <- readRDS(paste(OUT_DIR,"gtex_tissue_detail_vec_train.Rds",sep=""))

#TIS_NAME <- "Breast - Mammary Tissue"
TIS_NAME <- "Lung"

#ig_TIS_NAME <- "IG_Breast...Mammary.Tissue"
ig_TIS_NAME <- "IG_Lung"
#sal_TIS_NAME <- "Saliency_Breast...Mammary.Tissue"
sal_TIS_NAME <- "Saliency_Lung"

# generate mean transcript expression across tissue matrix
lung_samples <- which(gtex_tissue_detail.vec.train %in% TIS_NAME)

tis_pvals.nls <- list()
for(i in 1:nrow(vst.count.mtx.train)){
  w0 <- wilcox.test(x=as.numeric(vst.count.mtx.train[i,lung_samples]),
                    y=as.numeric(vst.count.mtx.train[i,-lung_samples]),alternative = "greater")
  tis_pvals.nls[[rownames(vst.count.mtx.train)[i]]] <- w0$p.value
  print(paste("caluclated p-value ",i," of ",nrow(vst.count.mtx.train),"...",sep=""))
}

# calculate p-values
z <- unlist(tis_pvals.nls)
transcript_pval.df <- data.frame(ens_n1 = names(z),
                                 fdr = p.adjust(z,method="fdr"))
#filter edges whose preceeding node is a reaction
tis_edges.df <- labelled_edge_weights.df %>%
  dplyr::filter(!Preceeding_Reaction %in% rownames(misclass_rates.df)) %>%
  dplyr::select(Preceeding_Reaction,
                Following_Reaction,
                ig_TIS_NAME,
                sal_TIS_NAME)
# perform pathway enrichment across mean transcript expression across tissue
ALPHA <- 0.05

#phyper(q, m, n, k, lower.tail = FALSE)
# q = number of white balls drawn (number of transcripts/reactions shared between selection & pathway) 
# m = number of white balls in urn (number of transcripts/reactions in pathway)
# n = number of black balls in urn (number of transcripts/reactions not in pathway)
# k = number of balls drawn (number of transcripts/reactions in selection)

# transcript-> pathway filenames
TRANSCRIPT_TO_PTHWY_FN <- "/home/burkhart/Software/reticula/data/aim1/input/Ensembl2Reactome_All_Levels.csv" # on box.com

# replace NA fdr values with machine minimum (underflow round-up)
transcript_pval.df$fdr[is.na(transcript_pval.df$fdr)] <- .Machine$double.xmin

transcript_2_pthwy.df <- read.table(file=TRANSCRIPT_TO_PTHWY_FN,sep="\t",header = FALSE,
                                    colClasses = c("character","character","NULL"))
colnames(transcript_2_pthwy.df) <- c("EnsemblID","Pathway")

# remove pathways not present in both reaction and transcript annotation files
# reactions have 1998 unique pathway annoations
# transcripts have 20831 unique pathway annotations
# reactions and transcripts share 1873 unique pathway annottions
shared_pathways <- intersect(as.character(tis_edges.df$Preceeding_Reaction),transcript_2_pthwy.df$Pathway)
transcript_2_pthwy.df <- transcript_2_pthwy.df %>% dplyr::filter(Pathway %in% shared_pathways)

# ensure above functions as expected
assertthat::are_equal(length(unique(shared_pathways)),
                      length(unique(transcript_2_pthwy.df$Pathway)))

# select intersection of significant reactions/transcripts and those with pathway annotations
shared_transcripts <- intersect(transcript_pval.df$ens_n1,transcript_2_pthwy.df$EnsemblID)

# filter selections by shared transcripts/reactions for consistency
transcript_pval.df.shared <- transcript_pval.df %>% dplyr::filter(ens_n1 %in% shared_transcripts)

# filter reactions/transcripts mapped to pathways
transcript_2_pthwy.df.shared <- transcript_2_pthwy.df %>% dplyr::filter(EnsemblID %in% shared_transcripts)
shared_pathways <- unique(transcript_2_pthwy.df.shared$Pathway)

# select significant reactions/transcipts
significant_transcripts.df <- transcript_pval.df.shared %>%
  dplyr::filter(fdr < ALPHA)

# build pathway lists
transcript_pathway_list <- list()
shared_pathway_2_n_transcripts.nls <- list()
for(i in 1:nrow(transcript_2_pthwy.df.shared)){
  # print(i) #debugging
  transcript <- transcript_2_pthwy.df.shared[i,1]
  pthwy <- transcript_2_pthwy.df.shared[i,2]
  if(is.null(transcript_pathway_list[[pthwy]])){
    transcript_pathway_list[[pthwy]] <- c(transcript)
    shared_pathway_2_n_transcripts.nls[[pthwy]] <- 1
  } else if(!(transcript %in% transcript_pathway_list[[pthwy]])) {
    transcript_pathway_list[[pthwy]] <- c(transcript_pathway_list[[pthwy]],transcript)
    shared_pathway_2_n_transcripts.nls[[pthwy]] <- shared_pathway_2_n_transcripts.nls[[pthwy]] + 1
  }
  else {
    print(paste("ERROR: Duplicate transcript detected on row ",i,
                ": transcript = ",transcript,
                ": pathway = ",pthwy,". Ignoring...",sep=""))
  }
}

print(length(shared_pathways))

transcript_pathway_enrichment <- list()
for(pthwy in shared_pathways){
  transcripts_in_pathway <- transcript_pathway_list[[pthwy]]
  q = length(intersect(significant_transcripts.df$ens_n1,transcripts_in_pathway))
  m = length(transcripts_in_pathway)
  n = length(unique(transcript_2_pthwy.df.shared$EnsemblID)) - length(transcripts_in_pathway)
  k = nrow(significant_transcripts.df)
  transcript_pathway_enrichment[[pthwy]] <- stats::phyper(q-1,m,n,k,lower.tail = FALSE)
}
# correlate pathway enrichment across tissue with edges whose preceeding node is a pathway
tis_edges.df <- tis_edges.df %>% dplyr::filter(as.character(Preceeding_Reaction) %in% names(transcript_pathway_enrichment))

tis_edges.df$transcript_enrichment_pval <- 1.0
for(i in 1:nrow(tis_edges.df)){
  print(i)
  preceeding_pthwy <- as.character(tis_edges.df[i,"Preceeding_Reaction"])
  print(preceeding_pthwy %in% names(transcript_pathway_enrichment))
  print(preceeding_pthwy)
  print(transcript_pathway_enrichment[[preceeding_pthwy]])
  tis_edges.df[i,"transcript_enrichment_pval"] <- transcript_pathway_enrichment[[preceeding_pthwy]]
}

tis_edges.df$transcript_enrichment_fdr <- p.adjust(tis_edges.df$transcript_enrichment_pval,method="fdr")

x_vector <- tis_edges.df[,ig_TIS_NAME]
y_vector <- -log10(tis_edges.df$transcript_enrichment_fdr)

plot(x=x_vector,y=y_vector)
abline(lm(y_vector~x_vector),col="red")
lines(lowess(x_vector,y_vector),col="blue")
cor.test(x_vector,y_vector)

#filter edges whose preceeding node is a pathway
tis_edges.df <- labelled_edge_weights.df %>%
  dplyr::filter(Preceeding_Reaction %in% rownames(misclass_rates.df)) %>%
  dplyr::select(Preceeding_Reaction,
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
### pathway hierarchy lung results below
###

# preceeding_top_ig vs preceeding_acc
#Pearson's product-moment correlation
#data:  x_vector and y_vector
#t = 24.321, df = 11137, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.2068665 0.2421358
#sample estimates:
#      cor 
#0.2245747 

# preceeding_top_saliency vs preceeding_acc
#Pearson's product-moment correlation
#data:  x_vector and y_vector
#t = 24.404, df = 11137, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.2075963 0.2428535
#sample estimates:
#      cor 
#0.2252986 

# preceeding_top_ig vs preceeding_ari
#Pearson's product-moment correlation
#data:  x_vector and y_vector
#t = 18.044, df = 11137, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.1504318 0.1865191
#sample estimates:
#      cor 
#0.1685319 

# preceeding_top_saliency vs preceeding_ari
#Pearson's product-moment correlation
#data:  x_vector and y_vector
#t = 18.447, df = 11137, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.1541069 0.1901479
#sample estimates:
#     cor 
#0.172185 

###
### pathway hierarchy breast results below
###

# preceeding_top_ig vs preceeding_acc
#Pearson's product-moment correlation
#data:  x_vector and y_vector
#t = 21.836, df = 11137, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.1847501 0.2203675
#sample estimates:
#      cor 
#0.2026258 

# preceeding_top_saliency vs preceeding_acc
#	Pearson's product-moment correlation
#data:  x_vector and y_vector
#t = 22.958, df = 11137, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  0.1947737 0.2302377
#sample estimates:
#  cor 
#0.2125757 

# preceeding_top_ig vs preceeding_ari
# Pearson's product-moment correlation
#data:  x_vector and y_vector
#t = 19.733, df = 11137, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.1657971 0.2016846
#sample estimates:
#      cor 
#0.1838021 

# preceeding_top_saliency vs preceeding_ari
#	Pearson's product-moment correlation
#data:  x_vector and y_vector
#t = 20.466, df = 11137, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  0.1724203 0.2082163
#sample estimates:
#  cor 
#0.1903816 
