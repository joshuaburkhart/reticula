set.seed(88888888)

library(magrittr)
library(ggplot2)
library(plotly)
library(dplyr)
library(stats)

start_time <- Sys.time()

ALPHA <- 0.05

#phyper(q, m, n, k, lower.tail = FALSE)
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

# remove pathways not present in both reaction and transcript annotation files
# reactions have 1998 unique pathway annoations
# transcripts have 20831 unique pathway annotations
# reactions and transcripts share 1873 unique pathway annottions
shared_pathways <- intersect(reaction_2_pthwy.df$Pathway,transcript_2_pthwy.df$Pathway)
reaction_2_pthwy.df <- reaction_2_pthwy.df %>% dplyr::filter(Pathway %in% shared_pathways)
transcript_2_pthwy.df <- transcript_2_pthwy.df %>% dplyr::filter(Pathway %in% shared_pathways)

# ensure above functions as expected
assertthat::are_equal(length(unique(reaction_2_pthwy.df$Pathway)),
                      length(unique(transcript_2_pthwy.df$Pathway)))

# select intersection of significant reactions/transcripts and those with pathway annotations
shared_reactions <- intersect(reaction_pval.df$rxn_n1,reaction_2_pthwy.df$ReactionlikeEvent)
shared_transcripts <- intersect(transcript_pval.df$ens_n1,transcript_2_pthwy.df$EnsemblID)

# filter selections by shared transcripts/reactions for consistency
reaction_pval.df.shared <- reaction_pval.df %>% dplyr::filter(rxn_n1 %in% shared_reactions)
transcript_pval.df.shared <- transcript_pval.df %>% dplyr::filter(ens_n1 %in% shared_transcripts)

# filter reactions/transcripts mapped to pathways
reaction_2_pthwy.df.shared <- reaction_2_pthwy.df %>% dplyr::filter(ReactionlikeEvent %in% shared_reactions)
transcript_2_pthwy.df.shared <- transcript_2_pthwy.df %>% dplyr::filter(EnsemblID %in% shared_transcripts)
shared_pathways <- intersect(reaction_2_pthwy.df.shared$Pathway,transcript_2_pthwy.df.shared$Pathway)

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
  } else if(!(reaction %in% reaction_pathway_list[[pthwy]])){
    reaction_pathway_list[[pthwy]] <- c(reaction_pathway_list[[pthwy]],reaction)
  }
  else {
    print(paste("ERROR: Duplicate reaction detected on row ",i,
                ": reaction = ",reaction,
                ": pathway = ",pthwy,". Ignoring...",sep=""))
  }
}

transcript_pathway_list <- list()
for(i in 1:nrow(transcript_2_pthwy.df.shared)){
  # print(i) #debugging
  transcript <- transcript_2_pthwy.df.shared[i,1]
  pthwy <- transcript_2_pthwy.df.shared[i,2]
  if(is.null(transcript_pathway_list[[pthwy]])){
    transcript_pathway_list[[pthwy]] <- c(transcript)
  } else if(!(transcript %in% transcript_pathway_list[[pthwy]])) {
    transcript_pathway_list[[pthwy]] <- c(transcript_pathway_list[[pthwy]],transcript)
  }
  else {
    print(paste("ERROR: Duplicate transcript detected on row ",i,
                ": transcript = ",transcript,
                ": pathway = ",pthwy,". Ignoring...",sep=""))
  }
}

# calculate enrichment for each reaction/transcript pathway
reaction_pathway_enrichment <- list()
for(pthwy in shared_pathways){
  #print(pthwy) #debugging
  reactions_in_pathway <- reaction_pathway_list[[pthwy]]
   q = length(intersect(significant_reactions.df$rxn_n1,reactions_in_pathway))
   m = length(reactions_in_pathway)
   n = length(unique(reaction_2_pthwy.df.shared$ReactionlikeEvent)) - length(reactions_in_pathway)
   k = nrow(significant_reactions.df)
   reaction_pathway_enrichment[[pthwy]] <- stats::phyper(q-1,m,n,k,lower.tail = FALSE)
}
saveRDS(reaction_pathway_enrichment,
        file="/home/burkhart/Software/reticula/data/aim1/output/reaction_pathway_enrichment.rds")

transcript_pathway_enrichment <- list()
for(pthwy in shared_pathways){
  #print(pthwy) #debugging
  transcripts_in_pathway <- transcript_pathway_list[[pthwy]]
  q = length(intersect(significant_transcripts.df$ens_n1,transcripts_in_pathway))
  m = length(transcripts_in_pathway)
  n = length(unique(transcript_2_pthwy.df.shared$EnsemblID)) - length(transcripts_in_pathway)
  k = nrow(significant_transcripts.df)
  transcript_pathway_enrichment[[pthwy]] <- stats::phyper(q-1,m,n,k,lower.tail = FALSE)
}
saveRDS(transcript_pathway_enrichment,
        file="/home/burkhart/Software/reticula/data/aim1/output/transcript_pathway_enrichment.rds")

# ensure both enrichment results are identically ordered and match shared_pathways
assertthat::are_equal(names(reaction_pathway_enrichment),
                      names(transcript_pathway_enrichment),
                      shared_pathways)

reaction_and_transcript_pathway_enrichment.df <- data.frame(Pathway = shared_pathways,
                                                            ReactionwisePathwayEnrichmentPVal = unlist(reaction_pathway_enrichment),
                                                            TranscriptwisePathwayEnrichmentPVal = unlist(transcript_pathway_enrichment))
library(metap)
reaction_and_transcript_pathway_enrichment.df <- reaction_and_transcript_pathway_enrichment.df %>%
  dplyr::rowwise() %>%
  dplyr::mutate(CombinedP = as.numeric((metap::sumlog(c(ReactionwisePathwayEnrichmentPVal,
                                                        TranscriptwisePathwayEnrichmentPVal)) %>% .[3])))#MaxP = max(ReactionwisePathwayEnrichmentPVal,TranscriptwisePathwayEnrichmentPVal))
reaction_and_transcript_pathway_enrichment.df$CombinedFDR <- p.adjust(reaction_and_transcript_pathway_enrichment.df$CombinedP,
                                                              method = "fdr")
reaction_and_transcript_pathway_enrichment.df$ReactionwisePathwayEnrichmentFDR <- p.adjust(reaction_and_transcript_pathway_enrichment.df$ReactionwisePathwayEnrichmentPVal,
                                                                      method = "fdr")
reaction_and_transcript_pathway_enrichment.df$TranscriptwisePathwyEnrichmentFDR <- p.adjust(reaction_and_transcript_pathway_enrichment.df$TranscriptwisePathwayEnrichmentPVal,
                                                                                           method = "fdr")
reaction_and_transcript_pathway_enrichment.df %>% write.csv(file="/home/burkhart/Software/reticula/data/aim1/output/reaction_and_transcript_pathway_enrichment_df.csv")

# horizontal and vertical lines set at significance threshold defined above
p <- ggplot(reaction_and_transcript_pathway_enrichment.df,
       aes(-log(ReactionwisePathwayEnrichmentPVal),
           -log(TranscriptwisePathwayEnrichmentPVal),
           label=Pathway)) + geom_point()
p <- p + geom_hline(yintercept=-log(ALPHA))
p <- p + geom_vline(xintercept=-log(ALPHA))
ggplotly(p)

end_time <- Sys.time()
print(paste(
  "Start: ",
  start_time,
  " End: ",
  end_time,
  " Difference: ",
  end_time - start_time
))

## add a double-check routine here for pathway component enrichemnt and then test both tails

# searching reaction pca files
fns <- c("full_rxn_pca_results_nls0-1000.Rds",
         "full_rxn_pca_results_nls1000-2000.Rds",
         "full_rxn_pca_results_nls2000-3000.Rds",
         "full_rxn_pca_results_nls3000-4000.Rds",
         "full_rxn_pca_results_nls4000-5000.Rds",
         "full_rxn_pca_results_nls5000-6000.Rds",
         "full_rxn_pca_results_nls6000-7000.Rds",
         "full_rxn_pca_results_nls7000-8000.Rds",
         "full_rxn_pca_results_nls8000-9000.Rds",
         "full_rxn_pca_results_nls9000-10000.Rds",
         "full_rxn_pca_results_nls.Rds")

printFileIdxWRxn <- function(rxn,fns){
  found <- -1
  for(i in seq(1:length(fns))){
    if(found > -1){
      break
    }
    print(paste("Searching file ",i,"...",sep=""))
    cur_fn <- fns[i]
    cur_obj <- readRDS(file=paste("/home/burkhart/Software/reticula/data/aim1/output/",cur_fn,sep=""))
    if(!(is.null(cur_obj[[rxn]]))){
      print(paste("PCA for ",rxn," found in file ",i,".",sep=""))
      found <- i
    }
    cur_obj <- list()
    gc()
  }
  return(found)
}

reaction_2_pthwy.df.shared %>%
  dplyr::filter(Pathway == "R-HSA-381070") ->
  pthwy_reactions.df

reaction_pval.df.shared %>%
  dplyr::filter(rxn_n1 %in% pthwy_reactions.df$ReactionlikeEvent)

transcript_2_pthwy.df.shared %>%
  dplyr::filter(Pathway == "R-HSA-381070") ->
  pthwy_transcripts.df

transcript_pval.df.shared %>%
  dplyr::filter(ens_n1 %in% pthwy_transcripts.df$EnsemblID)

roi <- "R-HSA-1791092"
rxn2ensembls.nls[[roi]]
                    
i <- printFileIdxWRxn(roi,fns)
pca_data <- readRDS(file=paste("/home/burkhart/Software/reticula/data/aim1/output/",fns[i],sep=""))
pca_obj <- pca_data[[roi]]

tissue_group_labels <- numeric()
tissue_group_labels[high_prolif_samples] <- 2
tissue_group_labels[med_prolif_samples] <- 1
tissue_group_labels[low_prolif_samples] <- 3
tissue_group_labels[which(is.na(tissue_group_labels))] <- -1

pca.df <- data.frame(pc1 = pca_obj$x[c(high_prolif_samples,med_prolif_samples,low_prolif_samples),1],
                     pc2 = pca_obj$x[c(high_prolif_samples,med_prolif_samples,low_prolif_samples),2],
                     grp = tissue_group_labels[c(high_prolif_samples,med_prolif_samples,low_prolif_samples)])

pca.df <- pca.df %>% dplyr::arrange(desc(grp))

pca.d <- data.frame(
  PC1 = pca.df$pc1,
  PC2 = pca.df$pc2,
  Section = pca.df$grp
)

plot_ly(
  x = pca.d$PC1,
  y = pca.d$PC2,
  type = "scatter",
  mode = "markers",
  color = pca.d$Section,
  size = 1
)
