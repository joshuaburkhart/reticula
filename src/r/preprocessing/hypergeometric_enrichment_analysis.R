set.seed(88888888)

library(magrittr)
library(ggplot2)
library(plotly)
library(dplyr)
library(stats)

start_time <- Sys.time()

OUT_DIR <- "/home/burkhart/Software/reticula/data/aim1/output/"
ALPHA <- 0.05

#phyper(q, m, n, k, lower.tail = FALSE)
# q = number of white balls drawn (number of transcripts/reactions shared between selection & pathway) 
# m = number of white balls in urn (number of transcripts/reactions in pathway)
# n = number of black balls in urn (number of transcripts/reactions not in pathway)
# k = number of balls drawn (number of transcripts/reactions in selection)

# transcript/reaction p-value filenames
REACTION_PVAL_FN <- paste(OUT_DIR,"combined_w_fisher.csv",sep="") # on box.com
TRANSCRIPT_PVAL_FN <- paste(OUT_DIR,"ens_combined_w_fisher.csv",sep="") # on box.com

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
        file=paste(OUT_DIR,"reaction_pathway_enrichment.rds",sep=""))

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
        file=paste(OUT_DIR,"transcript_pathway_enrichment.rds",sep=""))

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
reaction_and_transcript_pathway_enrichment.df %>% write.csv(file=paste(OUT_DIR,"reaction_and_transcript_pathway_enrichment_df.csv",sep=""))

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

rxn_id_2_result_file_idx.nls <- readRDS(file=paste(OUT_DIR,"rxn_id_2_result_file_idx_nls.Rds",sep=""))

printFileIdxWRxn <- function(rxn,fns){
    idx <- rxn_id_2_result_file_idx.nls[[rxn]]
    print(paste("Searching file ",i,"...",sep=""))
    cur_obj <- readRDS(file=paste(OUT_DIR,fns[idx],sep=""))
    if(!(is.null(cur_obj[[rxn]]))){
      print(paste("PCA for ",rxn," found in file ",idx,".",sep=""))
    }else{
      print(paste("Error: ",rxn," not found in file ",idx,"!",sep=""))
    }
    cur_obj <- NULL
    gc()
    return(idx)
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

rxn2ensembls.nls <- readRDS(paste(OUT_DIR, "rxn2ensembls_nls.Rds", sep = ""))
roi <- "R-HSA-1791092"
rxn2ensembls.nls[[roi]]
                    
i <- printFileIdxWRxn(roi,fns)
pca_data <- readRDS(file=paste(OUT_DIR,fns[i],sep=""))
pca_obj <- pca_data[[roi]]

high_prolif <- c("Stomach",
                 "Colon - Sigmoid",
                 "Colon - Transverse",
                 "Small Intestine - Terminal Ileum",
                 "Spleen",
                 "Esophagus - Mucosa",
                 "Esophagus - Gastroesophageal Junction")
med_prolif <- c("Bladder",
                "Skin - Not Sun Exposed (Suprapubic)",
                "Lung",
                "Liver",
                "Pancreas",
                "Kidney - Cortex")
low_prolif <- c("Adipose - Visceral (Omentum)",
                "Adipose - Subcutaneous",
                "Thyroid",
                "Muscle - Skeletal",
                "Heart - Left Ventricle",
                "Heart - Atrial Appendage",
                "Adrenal Gland")

gtex_tissue_detail.vec.train <- readRDS(paste(OUT_DIR,"gtex_tissue_detail_vec_train.Rds",sep=""))

high_prolif_samples <- which(gtex_tissue_detail.vec.train %in% high_prolif)
med_prolif_samples <- which(gtex_tissue_detail.vec.train %in% med_prolif)
low_prolif_samples <- which(gtex_tissue_detail.vec.train %in% low_prolif)

tissue_group_labels <- character()
tissue_group_labels[high_prolif_samples] <- "High"
tissue_group_labels[med_prolif_samples] <- "Med"
tissue_group_labels[low_prolif_samples] <- "Low"
tissue_group_labels[which(is.na(tissue_group_labels))] <- ""

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

# investigate transcripts with significant wilcoxon p-values
ens_ids <- c("ENSG00000105852",
             "ENSG00000137673",
             "ENSG00000149294",
             "ENSG00000169242",
             "ENSG00000180785",
             "ENSG00000151224",
             "ENSG00000164265",
             "ENSG00000115641",
             "ENSG00000115840",
             "ENSG00000140093",
             "ENSG00000140093",
             "ENSG00000119421",
             "ENSG00000115705",
             "ENSG00000160963",
             "ENSG00000140307",
             "ENSG00000166347",
             "ENSG00000078668",
             "ENSG00000173436",
             "ENSG00000104812",
             "ENSG00000134109",
             "ENSG00000178802",
             "ENSG00000181019",
             "ENSG00000166136",
             "ENSG00000188157",
             "ENSG00000117834",
             "ENSG00000164919",
             "ENSG00000197122",
             "ENSG00000185269",
             "ENSG00000162688",
             "ENSG00000165195",
             "ENSG00000130821",
             "ENSG00000164975",
             "ENSG00000132313",
             "ENSG00000077157",
             "ENSG00000258227",
             "ENSG00000110195",
             "ENSG00000123570",
             "ENSG00000136631",
             "ENSG00000160870",
             "ENSG00000234906",
             "ENSG00000168090")

vst.count.mtx.train <- readRDS(paste(OUT_DIR,"vst_count_mtx_train.Rds",sep=""))

for(ens_id in ens_ids){
  samples_by_count.df <- data.frame(
    expression_value = vst.count.mtx.train[ens_id, c(high_prolif_samples,
                                                                med_prolif_samples,
                                                                low_prolif_samples)] %>% t() %>% .[,1],
    tissue_group = factor(tissue_group_labels[c(high_prolif_samples,
                                                med_prolif_samples,
                                                low_prolif_samples)],
                          levels = c("Low","Med","High")))
  
  ggplot(samples_by_count.df, aes(x=expression_value,
                                  y=tissue_group,
                                  color=tissue_group)) +
    geom_violin() +
    coord_flip() +
    geom_boxplot(width  =0.15) +
    labs(title=paste(ens_id," expression across tissue groups"),
         x="Normalized Expression Value",
         y = "Mean Tissue Proliferation Rate")
  ggsave(paste(OUT_DIR,ens_id,"_expr_v_grps.png",sep=""),device = png())  
  dev.off()
}
