set.seed(88888888)

library(VennDiagram)
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
REACTION_PVAL_FN <- paste(OUT_DIR,"rxn_combined_w_fisher.csv",sep="") # on box.com
TRANSCRIPT_PVAL_FN <- paste(OUT_DIR,"ens_combined_w_fisher.csv",sep="") # on box.com

# transcript/reaction -> pathway filenames
REACTION_TO_PTHWY_FN <- "/home/burkhart/Software/reticula/data/aim1/input/ReactionToPathway_Rel_71_122820.csv" # on box.com
TRANSCRIPT_TO_PTHWY_FN <- "/home/burkhart/Software/reticula/data/aim1/input/Ensembl2Reactome_All_Levels.csv" # on box.com

# load as dataframes
reaction_pval.df <- read.table(file=REACTION_PVAL_FN,sep=",",header = TRUE,
                               colClasses = c("NULL","character","NULL","character","numeric","NULL","NULL","NULL","NULL","numeric"))
transcript_pval.df <- read.table(file=TRANSCRIPT_PVAL_FN,sep=",",header = TRUE,
                                 colClasses = c("NULL","character","NULL","character","numeric","NULL","NULL","NULL","numeric"))

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

# calculate transcripts per reaction per pathway
transcripts_in_rxn <- readRDS(file=paste(OUT_DIR,"transcript_in_rxn_nls.Rds",sep=""))
mean_transcripts_per_reaction_per_pathay <- list()
for(pathway in shared_pathways){
  reactions_in_pathway <- reaction_2_pthwy.df.shared %>%
    dplyr::filter(pathway == Pathway) %>% .$ReactionlikeEvent
  sum_transcripts_per_reaction_in_pathway = 0
  for(reaction in reactions_in_pathway){
    sum_transcripts_per_reaction_in_pathway = sum_transcripts_per_reaction_in_pathway + transcripts_in_rxn[[reaction]]
  }
  mean_transcripts_per_reaction_per_pathay[[pathway]] <- sum_transcripts_per_reaction_in_pathway / length(reactions_in_pathway)
}

# select significant reactions/transcipts
significant_reactions.df <- reaction_pval.df.shared %>%
  dplyr::filter(fdr < ALPHA) %>%
  dplyr::filter(direction == "positive" | direction == "negative")
significant_transcripts.df <- transcript_pval.df.shared %>%
  dplyr::filter(fdr < ALPHA) %>%
  dplyr::filter(direction == "positive" | direction == "negative")

# build pathway lists
reaction_pathway_list <- list()
shared_pathway_2_n_reactions.nls <- list()
for(i in 1:nrow(reaction_2_pthwy.df.shared)){
  # print(i) #debugging
  reaction <- reaction_2_pthwy.df.shared[i,1]
  pthwy <- reaction_2_pthwy.df.shared[i,2]
  if(is.null(reaction_pathway_list[[pthwy]])){
    reaction_pathway_list[[pthwy]] <- c(reaction)
    shared_pathway_2_n_reactions.nls[[pthwy]] <- 1
  } else if(!(reaction %in% reaction_pathway_list[[pthwy]])){
    reaction_pathway_list[[pthwy]] <- c(reaction_pathway_list[[pthwy]],reaction)
    shared_pathway_2_n_reactions.nls[[pthwy]] <- shared_pathway_2_n_reactions.nls[[pthwy]] + 1
  }
  else {
    print(paste("ERROR: Duplicate reaction detected on row ",i,
                ": reaction = ",reaction,
                ": pathway = ",pthwy,". Ignoring...",sep=""))
  }
}

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

# store reations/transcripts per pathway
saveRDS(shared_pathway_2_n_reactions.nls,file=paste(OUT_DIR,"shared_pathway_2_n_reactions_nls.Rds",sep=""))
saveRDS(shared_pathway_2_n_transcripts.nls,file=paste(OUT_DIR,"shared_pathway_2_n_transcripts_nls.Rds",sep=""))

# calculate enrichment for each reaction/transcript pathway
reaction_pathway_enrichment <- list()
transcript_pathway_enrichment <- list()
for(pthwy in shared_pathways){
  #print(pthwy) #debugging
  reactions_in_pathway <- reaction_pathway_list[[pthwy]]
  q = length(intersect(significant_reactions.df$rxn_n1,reactions_in_pathway))
  m = length(reactions_in_pathway)
  n = length(unique(reaction_2_pthwy.df.shared$ReactionlikeEvent)) - length(reactions_in_pathway)
  k = nrow(significant_reactions.df)
  reaction_pathway_enrichment[[pthwy]] <- stats::phyper(q-1,m,n,k,lower.tail = FALSE)

  transcripts_in_pathway <- transcript_pathway_list[[pthwy]]
  q = length(intersect(significant_transcripts.df$ens_n1,transcripts_in_pathway))
  m = length(transcripts_in_pathway)
  n = length(unique(transcript_2_pthwy.df.shared$EnsemblID)) - length(transcripts_in_pathway)
  k = nrow(significant_transcripts.df)
  transcript_pathway_enrichment[[pthwy]] <- stats::phyper(q-1,m,n,k,lower.tail = FALSE)
}
saveRDS(reaction_pathway_enrichment,
        file=paste(OUT_DIR,"reaction_pathway_enrichment.rds",sep=""))
saveRDS(transcript_pathway_enrichment,
        file=paste(OUT_DIR,"transcript_pathway_enrichment.rds",sep=""))

# ensure both enrichment results are identically ordered and match shared_pathways
assertthat::are_equal(names(reaction_pathway_enrichment),
                      names(transcript_pathway_enrichment),
                      shared_pathways)

common_pathway_2_n_reactions <- shared_pathway_2_n_reactions.nls[shared_pathways]
common_pathway_2_n_transcripts <- shared_pathway_2_n_transcripts.nls[shared_pathways]

reaction_and_transcript_pathway_enrichment.df <- data.frame(Pathway = shared_pathways,
                                                            ReactionwisePathwayEnrichmentPVal = unlist(reaction_pathway_enrichment),
                                                            TranscriptwisePathwayEnrichmentPVal = unlist(transcript_pathway_enrichment),
                                                            ReactionsInPathway = unlist(common_pathway_2_n_reactions),
                                                            TranscriptsInPathway = unlist(common_pathway_2_n_transcripts),
                                                            MeanTranscriptsPerReactionInPathway = unlist(mean_transcripts_per_reaction_per_pathay))
library(metap)
reaction_and_transcript_pathway_enrichment.df <- reaction_and_transcript_pathway_enrichment.df %>%
  dplyr::rowwise() %>%
  dplyr::mutate(CombinedP = as.numeric((metap::sumlog(c(ReactionwisePathwayEnrichmentPVal,
                                                        TranscriptwisePathwayEnrichmentPVal)) %>% .[3])))#MaxP = max(ReactionwisePathwayEnrichmentPVal,TranscriptwisePathwayEnrichmentPVal))
reaction_and_transcript_pathway_enrichment.df$CombinedFDR <- p.adjust(reaction_and_transcript_pathway_enrichment.df$CombinedP,
                                                              method = "fdr")
reaction_and_transcript_pathway_enrichment.df$ReactionwisePathwayEnrichmentFDR <- p.adjust(reaction_and_transcript_pathway_enrichment.df$ReactionwisePathwayEnrichmentPVal,
                                                                      method = "fdr")
reaction_and_transcript_pathway_enrichment.df$TranscriptwisePathwayEnrichmentFDR <- p.adjust(reaction_and_transcript_pathway_enrichment.df$TranscriptwisePathwayEnrichmentPVal,
                                                                                           method = "fdr")
reaction_and_transcript_pathway_enrichment.df %>% write.csv(file=paste(OUT_DIR,"reaction_and_transcript_pathway_enrichment_df.csv",sep=""))

pathways_significantly_enriched_by_transcripts <- reaction_and_transcript_pathway_enrichment.df %>%
  dplyr::filter(TranscriptwisePathwayEnrichmentFDR < ALPHA)
pathways_significantly_enriched_by_reactions <- reaction_and_transcript_pathway_enrichment.df %>%
  dplyr::filter(ReactionwisePathwayEnrichmentFDR < ALPHA)
pathways_significantly_enriched_by_both <- reaction_and_transcript_pathway_enrichment.df %>%
  dplyr::filter(TranscriptwisePathwayEnrichmentFDR < ALPHA) %>%
  dplyr::filter(ReactionwisePathwayEnrichmentFDR < ALPHA)
pathways_not_enriched_by_transcripts <- reaction_and_transcript_pathway_enrichment.df %>%
  dplyr::filter(TranscriptwisePathwayEnrichmentFDR >= ALPHA)

#svg(filename=paste(OUT_DIR,"PathwayEnrichmentVennDiagram.svg",sep=""),
#    width=15,
#    height=15,
#    pointsize=12)

q <- nrow(pathways_significantly_enriched_by_both)
m <- nrow(pathways_significantly_enriched_by_transcripts)
n <- nrow(pathways_not_enriched_by_transcripts)
k <- nrow(pathways_significantly_enriched_by_reactions)

g <- VennDiagram::draw.pairwise.venn(area1=m,
                              area2=k,
                              category = c("Pathways Enrichmed by Transcripts","Pathways Enriched by Reactions"),
                              cross.area=q,
                              fill=c("red","gold"))
require(gridExtra)
grid.arrange(gTree(children=g),
             top="Reaction vs Transcript Pathway Enrichment",
             bottom=paste("Overlap p-value = ",
                          phyper(q-1,m,n,k, lower.tail = FALSE),
                          sep=""))

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
    print(paste("Searching file ",idx,"...",sep=""))
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
  dplyr::filter(Pathway == "R-HSA-2029480") ->
  pthwy_reactions.df

p <- reaction_pval.df.shared
p <- p %>%
  dplyr::arrange(rxn_n1 %in% pthwy_reactions.df$ReactionlikeEvent) %>%
  dplyr::filter(direction != "none")
ggplot(p,
       aes(x = difference,
           y = -log10(fdr),
             colour = rxn_n1 %in% pthwy_reactions.df$ReactionlikeEvent,
             label = rxn_n1)) +
  geom_point(size = 2) +
  scale_colour_manual(values = c("grey","#F8766D")) +
  expand_limits(x = c(-6,12)) +
  geom_label(check_overlap=TRUE,
            hjust = -.05,
            vjust = .4,
            data=subset(p, fdr < ALPHA & rxn_n1 %in% pthwy_reactions.df$ReactionlikeEvent)) +
   geom_hline(yintercept=-log10(ALPHA)) +
  theme_bw() +
  theme(legend.position = "none")

#vst.count.mtx.train <- readRDS(paste(OUT_DIR,"vst_count_mtx_train.Rds",sep=""))
#rxn_pca.nls <- readRDS(paste(OUT_DIR,"rxn_pca_nls.Rds",sep=""))

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

i <- printFileIdxWRxn("R-HSA-2029458",fns)
pca_data <- readRDS(file=paste(OUT_DIR,fns[i],sep=""))
pca_obj <- pca_data[["R-HSA-2029458"]]

pca.df <- data.frame(pc1 = pca_obj$x[c(high_prolif_samples,med_prolif_samples,low_prolif_samples),1],
                     pc2 = pca_obj$x[c(high_prolif_samples,med_prolif_samples,low_prolif_samples),2],
                     pc3 = pca_obj$x[c(high_prolif_samples,med_prolif_samples,low_prolif_samples),3],
                     grp = tissue_group_labels[c(high_prolif_samples,med_prolif_samples,low_prolif_samples)])

pca.df <- pca.df %>% dplyr::arrange(desc(grp))

pca.d <- data.frame(
  PC1 = pca.df$pc1,
  PC2 = pca.df$pc2,
  PC3 = pca.df$pc3,
  ProlifGrp = pca.df$grp
)

p <- plot_ly(
  x = pca.d$PC1,
  y = pca.d$PC2,
  z = pca.d$PC3,
  type = "scatter3d",
  mode = "markers",
  color = pca.d$ProlifGrp,
  size = 1
)

p

transcript_2_pthwy.df.shared %>%
  dplyr::filter(Pathway == "R-HSA-2029480") ->
  pthwy_transcripts.df

p <- transcript_pval.df.shared

p <- p %>%
  dplyr::arrange(ens_n1 %in% pthwy_transcripts.df$EnsemblID) %>%
  dplyr::filter(direction != "none")
ggplot(p,
       aes(x = difference,
           y = -log10(fdr),
           colour = ens_n1 %in% pthwy_transcripts.df$EnsemblID,
           label = ens_n1)) +
  geom_point(size = 2) +
  scale_colour_manual(values = c("grey","#F8766D")) +
  expand_limits(x = c(-6,12)) +
  geom_label(check_overlap=TRUE,
             hjust = -.05,
             vjust = .4,
             data=subset(p, fdr < ALPHA & ens_n1 %in% pthwy_transcripts.df$EnsemblID)) +
  geom_hline(yintercept=-log10(ALPHA)) +
  theme_bw() +
  theme(legend.position = "none")

rxn2ensembls.nls <- readRDS(paste(OUT_DIR, "rxn2ensembls_nls.Rds", sep = ""))

rois <- c("R-HSA-2029458")

library(orca)
for(roi in rois){
  rxn2ensembls.nls[[roi]]
  i <- printFileIdxWRxn(roi,fns)
  pca_data <- readRDS(file=paste(OUT_DIR,fns[i],sep=""))
  pca_obj <- pca_data[[roi]]
  
  pca.df <- data.frame(pc1 = pca_obj$x[c(high_prolif_samples,med_prolif_samples,low_prolif_samples),1],
                       pc2 = pca_obj$x[c(high_prolif_samples,med_prolif_samples,low_prolif_samples),2],
                       grp = tissue_group_labels[c(high_prolif_samples,med_prolif_samples,low_prolif_samples)])
  
  pca.df <- pca.df %>% dplyr::arrange(desc(grp))
  
  pca.d <- data.frame(
    PC1 = pca.df$pc1,
    PC2 = pca.df$pc2,
    Section = pca.df$grp
  )
  
  p <- plot_ly(
    x = pca.d$PC1,
    y = pca.d$PC2,
    type = "scatter",
    mode = "markers",
    color = pca.d$Section,
    size = 1
  )
  orca(p,file=paste(roi,"_pca.png",sep="")) #for some reason this crashes when OUT_DIR is specified... beyond the scope of this project
}

# investigate transcripts with significant wilcoxon p-values, update with significant positive/negative association directions
ens_ids <- c("ENSG00000197122")

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
         y = "Mean Tissue Proliferation Rate") +
    theme_bw() +
    theme(legend.position = "none") +
    scale_color_manual(values = c("#ff7f0e","#1f77b4","#2ca02c"))
  
  kruskal.test(expression_value ~ tissue_group, data = samples_by_count.df)
  pairwise.wilcox.test(samples_by_count.df$expression_value,
                       samples_by_count.df$tissue_group,
                       p.adjust.method = "fdr")
  
  ggsave(paste(OUT_DIR,ens_id,"_expr_v_grps.png",sep=""),device = png())  
  dev.off()
}

# measure reactions/transcripts per pathway distributions across cartesion quadrants
reactions_and_transcripts_per_pathway_quadrants.df <- reaction_and_transcript_pathway_enrichment.df %>%
  dplyr::mutate(quadrant = ifelse(ReactionwisePathwayEnrichmentFDR < ALPHA &
                                    TranscriptwisePathwayEnrichmentFDR < ALPHA, "1",
                                  ifelse(ReactionwisePathwayEnrichmentFDR >= ALPHA &
                                           TranscriptwisePathwayEnrichmentFDR < ALPHA, "2",
                                         ifelse(ReactionwisePathwayEnrichmentFDR >= ALPHA &
                                                  TranscriptwisePathwayEnrichmentFDR >= ALPHA, "3",
                                                ifelse(ReactionwisePathwayEnrichmentFDR < ALPHA &
                                                         TranscriptwisePathwayEnrichmentFDR >= ALPHA, "4",
                                                       "-1")))))

# horizontal and vertical lines set at significance threshold defined above
p <- ggplot(reactions_and_transcripts_per_pathway_quadrants.df,
            aes(-log10(ReactionwisePathwayEnrichmentFDR),
                -log10(TranscriptwisePathwayEnrichmentFDR),
                color = quadrant,
                label=Pathway)) +
  geom_text(check_overlap=TRUE,
            hjust = -.05,
            vjust = .4,
            angle = -65,
             data=subset(reactions_and_transcripts_per_pathway_quadrants.df, quadrant == 1)) +
  geom_point(size=2) +
  scale_y_continuous(trans = "log10") +
  scale_x_continuous(trans = "log10") +
  theme_bw() +
  theme(legend.position = "none")
p <- p + geom_hline(yintercept=-log10(ALPHA))
p <- p + geom_vline(xintercept=-log10(ALPHA))
p

-log10(ALPHA)

ggplot(reactions_and_transcripts_per_pathway_quadrants.df, aes(x=log10(ReactionsInPathway),
                                y=quadrant,
                                color=quadrant)) +
  geom_violin() +
  coord_flip() +
  geom_boxplot(width  =0.15) +
  labs(title=paste("Reactions per pathway across sections"),
       x="Log10 reactions per pathway",
       y = "Section") +
theme_bw() +
  theme(legend.position = "none")

kruskal.test(ReactionsInPathway ~ quadrant, data = reactions_and_transcripts_per_pathway_quadrants.df)
pairwise.wilcox.test(reactions_and_transcripts_per_pathway_quadrants.df$ReactionsInPathway,
                     reactions_and_transcripts_per_pathway_quadrants.df$quadrant,p.adjust.method = "fdr")


ggsave(paste(OUT_DIR,"Reactions_per_pathway_across_quadrants.png",sep=""),device = png())  
dev.off()

ggplot(reactions_and_transcripts_per_pathway_quadrants.df, aes(x=log10(TranscriptsInPathway),
                                                               y=quadrant,
                                                               color=quadrant)) +
  geom_violin() +
  coord_flip() +
  geom_boxplot(width  =0.15) +
  labs(title=paste("Transcripts per pathway across quadrants"),
       x="Log10 transcripts per pathway",
       y = "Section") +
  theme_bw() +
  theme(legend.position = "none")

kruskal.test(TranscriptsInPathway ~ quadrant, data = reactions_and_transcripts_per_pathway_quadrants.df)
pairwise.wilcox.test(reactions_and_transcripts_per_pathway_quadrants.df$TranscriptsInPathway,
                     reactions_and_transcripts_per_pathway_quadrants.df$quadrant,p.adjust.method = "fdr")

ggsave(paste(OUT_DIR,"Transcripts_per_pathway_across_quadrants.png",sep=""),device = png())  
dev.off()

ggplot(reactions_and_transcripts_per_pathway_quadrants.df, aes(x=log10(MeanTranscriptsPerReactionInPathway),
                                                               y=quadrant,
                                                               color=quadrant)) +
  geom_violin() +
  coord_flip() +
  geom_boxplot(width  =0.15) +
  labs(title=paste("Mean transcripts per reaction per pathway across quadrants"),
       x="Log10 mean transcripts per reaction per pathway",
       y = "Section") +
theme_bw() +
  theme(legend.position = "none")

ggsave(paste(OUT_DIR,"Mean_transcripts_per_reaction_per_pathway_across_quadrants.png",sep=""),device = png())  
dev.off()

kruskal.test(MeanTranscriptsPerReactionInPathway ~ quadrant, data = reactions_and_transcripts_per_pathway_quadrants.df)
pairwise.wilcox.test(reactions_and_transcripts_per_pathway_quadrants.df$MeanTranscriptsPerReactionInPathway,
                     reactions_and_transcripts_per_pathway_quadrants.df$quadrant,p.adjust.method = "fdr")
