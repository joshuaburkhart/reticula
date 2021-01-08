set.seed(88888888)

library(magrittr)
library(dplyr)
library(stats)

start_time <- Sys.time()

ALPHA <- 0.05

#phyper(q, m, n, k)
# q = number of white balls drawn (number of transcripts/reactions shared between selection & pathway) 
# m = number of white balls in urn (number of transcripts/reactions in pathway)
# n = number of black balls in urn (number of transcripts/reactions not in pathway)
# k = number of balls drawn (number of transcripts/reactions in selection)

# transcript/reaction p-value filenames
REACTION_PVAL_FN <- "/home/burkhart/Software/reticula/data/aim1/output/combined_w_fisher.csv" # on box.com
TRANSCRIPT_PVAL_FN <- "/home/burkhart/Software/reticula/data/aim1/output/ens_combined_w_fisher.csv" # on box.com

# transcript/reaction -> pathway filenames
REACTION_TO_PTHWY_FN <- "/home/burkhart/Software/reticula/data/aim1/input/ReactionToPathway_Rel_71_122820.csv" # on box.com
TRANSCRIPT_TO_PTHWY_FN <- "/home/burkhart/Software/reticula/data/aim1/input/Ensembl2Reactome_All_Levels.csv" # on box.com

# load as dataframes
reaction_pval.df <- read.table(file=REACTION_PVAL_FN,sep=",",header = TRUE,
                               colClasses = c("NULL","character","NULL","NULL","NULL","NULL","numeric"))
transcript_pval.df <- read.table(file=TRANSCRIPT_PVAL_FN,sep=",",header = TRUE,
                                 colClasses = c("NULL","character","NULL","NULL","NULL","NULL","numeric"))

# replace NA fdr values with machine minimum (underflow round-up)
reaction_pval.df$fdr[is.na(reaction_pval.df$fdr)] <- .Machine$double.xmin
transcript_pval.df$fdr[is.na(transcript_pval.df$fdr)] <- .Machine$double.xmin

reaction_2_pthwy.df <- read.table(file=REACTION_TO_PTHWY_FN,sep="\t",header = TRUE,
                                  colClasses = c("character","character"))
transcript_2_pthwy.df <- read.table(file=TRANSCRIPT_TO_PTHWY_FN,sep="\t",header = FALSE,
                                    colClasses = c("character","character","NULL"))
colnames(transcript_2_pthwy.df) <- c("EnsemblID","Pathway")

# select intersection of significant transcripts/reactions and those with pathway annotations
shared_reactions <- intersect(reaction_pval.df$rxn_n1,reaction_2_pthwy.df$ReactionlikeEvent)
shared_transcripts <- intersect(transcript_pval.df$ens_n1,transcript_2_pthwy.df$EnsemblID)

# filter selections by shared transcripts/reactions for consistency
reaction_pval.df.shared <- reaction_pval.df %>% dplyr::filter(rxn_n1 %in% shared_reactions)
transcript_pval.df.shared <- transcript_pval.df %>% dplyr::filter(ens_n1 %in% shared_transcripts)

reaction_2_pthwy.df.shared <- reaction_2_pthwy.df %>% dplyr::filter(ReactionlikeEvent %in% shared_reactions)
transcript_2_pthwy.df.shared <- transcript_2_pthwy.df %>% dplyr::filter(EnsemblID %in% shared_transcripts)

# select significant reactions/transcipts
significant_reactions.df <- reaction_pval.df.shared %>% dplyr::filter(fdr < ALPHA)
significant_transcripts.df <- transcript_pval.df.shared %>% dplyr::filter(fdr < ALPHA)

# build pathway lists
reaction_pathway_list <- list()
for(i in 1:nrow(reaction_2_pthwy.df.shared)){
  # print(i) #debugging
  reaction <- reaction_2_pthwy.df.shared[i,1]
  pthwy <- reaction_2_pthwy.df.shared[i,2]
  if(is.null(reaction_pathway_list[[pthwy]])){
    reaction_pathway_list[[pthwy]] <- c(reaction)
  } else {
    reaction_pathway_list[[pthwy]] <- c(reaction_pathway_list[[pthwy]],reaction)
  }
  if(length(reaction_pathway_list[[pthwy]]) != length(unique(reaction_pathway_list[[pthwy]]))){
    print(paste("ERROR: Duplicate reaction detected on row ",i,
                ": reaction = ",reaction,
                ": pathway = ",pthwy,sep=""))
  }
}

transcript_pathway_list <- list()
for(i in 1:nrow(transcript_2_pthwy.df.shared)){
  # print(i) #debugging
  transcript <- transcript_2_pthwy.df.shared[i,1]
  pthwy <- transcript_2_pthwy.df.shared[i,2]
  if(is.null(transcript_pathway_list[[pthwy]])){
    transcript_pathway_list[[pthwy]] <- c(transcript)
  } else {
    transcript_pathway_list[[pthwy]] <- c(transcript_pathway_list[[pthwy]],transcript)
  }
  if(length(transcript_pathway_list[[pthwy]]) != length(unique(transcript_pathway_list[[pthwy]]))){
    print(paste("ERROR: Duplicate transcript detected on row ",i,
                ": transcript = ",transcript,
                ": pathway = ",pthwy,sep=""))
  }
}

# calculate enrichment for each reaction/transcript pathway
reaction_pathway_enrichment <- list()
for(pthwy in names(reaction_pathway_list)){
  #print(pthwy) #debugging
  reactions_in_pathway <- reaction_pathway_list[[pthwy]]
   q = length(intersect(significant_reactions.df$rxn_n1,reactions_in_pathway))
   m = length(reactions_in_pathway)
   n = length(reaction_2_pthwy.df.shared) - length(reactions_in_pathway)
   k = length(significant_reactions.df)
   reaction_pathway_enrichment[[pthwy]] <- stats::phyper(q,m,n,k)
}
saveRDS(reaction_pathway_enrichment,
        file="/home/burkhart/Software/reticula/data/aim1/output/reaction_pathway_enrichment.rds")

transcript_pathway_enrichment <- list()
for(pthwy in names(transcript_pathway_list)){
  #print(pthwy) #debugging
  transcripts_in_pathway <- transcript_pathway_list[[pthwy]]
  q = length(intersect(significant_transcripts.df$rxn_n1,transcripts_in_pathway))
  m = length(transcripts_in_pathway)
  n = length(transcript_2_pthwy.df.shared) - length(transcripts_in_pathway)
  k = length(significant_transcripts.df)
  transcript_pathway_enrichment[[pthwy]] <- stats::phyper(q,m,n,k)
}
saveRDS(transcript_pathway_enrichment,
        file="/home/burkhart/Software/reticula/data/aim1/output/transcript_pathway_enrichment.rds")

end_time <- Sys.time()
print(paste(
  "Start: ",
  start_time,
  " End: ",
  end_time,
  " Difference: ",
  end_time - start_time
))

