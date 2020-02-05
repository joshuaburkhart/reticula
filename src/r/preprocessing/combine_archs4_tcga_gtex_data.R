library(dplyr)
library(biomaRt)
library(viridis)
library(magrittr)
library(pheatmap)
library(RColorBrewer)
library(SummarizedExperiment)

ensembl_dataset <- useEnsembl(biomart="ensembl",
                              dataset = "hsapiens_gene_ensembl",
                              mirror="uswest")

#ARCHS4_DATA_DIR <- "/home/users/burkhajo/WuLab/WuLabLustreDir/reticula/input/ARCHS4/"
ARCHS4_DATA_DIR <- "/Users/burkhajo/Software/reticula/data/aim1/input/ARCHS4/"
ARCHS4_ESGS_FIL <- "Esophagus_expression_matrix.tsv"
ARCHS4_COLN_FIL <- "Colon_expression_matrix.tsv"
ARCHS4_LIVR_FIL <- "Liver_expression_matrix.tsv"
ARCHS4_BRST_FIL <- "Breast_expression_matrix.tsv"
ARCHS4_KDNY_FIL <- "Kidney_expression_matrix.tsv"

archs4.esgs.df <- read.table(paste(ARCHS4_DATA_DIR,ARCHS4_ESGS_FIL,sep=""),header = TRUE,row.names = 1)
archs4.coln.df <- read.table(paste(ARCHS4_DATA_DIR,ARCHS4_COLN_FIL,sep=""),header = TRUE,row.names = 1)
archs4.livr.df <- read.table(paste(ARCHS4_DATA_DIR,ARCHS4_LIVR_FIL,sep=""),header = TRUE,row.names = 1)
archs4.brst.df <- read.table(paste(ARCHS4_DATA_DIR,ARCHS4_BRST_FIL,sep=""),header = TRUE,row.names = 1)
archs4.kdny.df <- read.table(paste(ARCHS4_DATA_DIR,ARCHS4_KDNY_FIL,sep=""),header = TRUE,row.names = 1)

archs4.df <- archs4.esgs.df %>%
  dplyr::bind_cols(archs4.coln.df) %>%
  dplyr::bind_cols(archs4.livr.df) %>%
  dplyr::bind_cols(archs4.brst.df) %>%
  dplyr::bind_cols(archs4.kdny.df)

rownames(archs4.df) <- rownames(archs4.kdny.df) # all archs4 dataframe rows match

archs4.ensembl_gene_map <- getBM(attributes = c("ensembl_gene_id_version","external_gene_name"),
                                 filters="external_gene_name",
                                 values=row.names(archs4.df),
                                 mart=ensembl_dataset)

#TCGA_DATA_DIR <- "/home/users/burkhajo/WuLab/WuLabLustreDir/reticula/input/recount2/recount2edTCGA/"
TCGA_DATA_DIR <- "/Users/burkhajo/Software/reticula/data/aim1/input/recount2/recount2edTCGA/"
TCGA_ESGS_FIL <- "rse_gene_esophagus.Rdata"
TCGA_COLN_FIL <- "rse_gene_colorectal.Rdata"
TCGA_LIVR_FIL <- "rse_gene_liver.Rdata"
TCGA_BRST_FIL <- "rse_gene_breast.Rdata"
TCGA_KDNY_FIL <- "rse_gene_kidney.Rdata"

load(paste(TCGA_DATA_DIR,TCGA_ESGS_FIL,sep=""))
tcga.esgs.cols <- rse_gene %>% colData()
tcga.esgs.df <- rse_gene %>% SummarizedExperiment::assay() %>% .[,tcga.esgs.cols$gdc_cases.samples.sample_type == "Solid Tissue Normal"] %>% as.data.frame()
load(paste(TCGA_DATA_DIR,TCGA_COLN_FIL,sep=""))
tcga.coln.cols <- rse_gene %>% colData()
tcga.coln.df <- rse_gene %>% SummarizedExperiment::assay() %>% .[,tcga.coln.cols$gdc_cases.samples.sample_type == "Solid Tissue Normal"] %>% as.data.frame()
load(paste(TCGA_DATA_DIR,TCGA_LIVR_FIL,sep=""))
tcga.livr.cols <- rse_gene %>% colData()
tcga.livr.df <- rse_gene %>% SummarizedExperiment::assay() %>% .[,tcga.livr.cols$gdc_cases.samples.sample_type == "Solid Tissue Normal"] %>% as.data.frame()
load(paste(TCGA_DATA_DIR,TCGA_BRST_FIL,sep=""))
tcga.brst.cols <- rse_gene %>% colData()
tcga.brst.df <- rse_gene %>% SummarizedExperiment::assay() %>% .[,tcga.brst.cols$gdc_cases.samples.sample_type == "Solid Tissue Normal"] %>% as.data.frame()
load(paste(TCGA_DATA_DIR,TCGA_KDNY_FIL,sep=""))
tcga.kdny.cols <- rse_gene %>% colData()
tcga.kdny.df <- rse_gene %>% SummarizedExperiment::assay() %>% .[,tcga.kdny.cols$gdc_cases.samples.sample_type == "Solid Tissue Normal"] %>% as.data.frame()

tcga.df <- tcga.esgs.df %>%
  dplyr::bind_cols(tcga.coln.df) %>%
  dplyr::bind_cols(tcga.livr.df) %>%
  dplyr::bind_cols(tcga.brst.df) %>%
  dplyr::bind_cols(tcga.kdny.df)

rownames(tcga.df) <- rownames(tcga.kdny.df) # all tcga dataframe rows match

tcga.ensembl_gene_map <- getBM(attributes = c("ensembl_gene_id_version","external_gene_name"),
                               filters="ensembl_gene_id_version",
                               values=row.names(tcga.df),
                               mart=ensembl_dataset)

#GTEx_DATA_DIR <- "/home/users/burkhajo/WuLab/WuLabLustreDir/reticula/input/recount2/recount2edGTEx/"
GTEx_DATA_DIR <- "/Users/burkhajo/Software/reticula/data/aim1/input/recount2/recount2edGTEx/"
GTEx_ESGS_FIL <- "rse_gene_esophagus.Rdata"
GTEx_COLN_FIL <- "rse_gene_colon.Rdata"
GTEx_LIVR_FIL <- "rse_gene_liver.Rdata"
GTEx_BRST_FIL <- "rse_gene_breast.Rdata"
GTEx_KDNY_FIL <- "rse_gene_kidney.Rdata"

load(paste(GTEx_DATA_DIR,GTEx_ESGS_FIL,sep=""))
gtex.esgs.df <- rse_gene %>% SummarizedExperiment::assay() %>% as.data.frame()
load(paste(GTEx_DATA_DIR,GTEx_COLN_FIL,sep=""))
gtex.coln.df <- rse_gene %>% SummarizedExperiment::assay() %>% as.data.frame()
load(paste(GTEx_DATA_DIR,GTEx_LIVR_FIL,sep=""))
gtex.livr.df <- rse_gene %>% SummarizedExperiment::assay() %>% as.data.frame()
load(paste(GTEx_DATA_DIR,GTEx_BRST_FIL,sep=""))
gtex.brst.df <- rse_gene %>% SummarizedExperiment::assay() %>% as.data.frame()
load(paste(GTEx_DATA_DIR,GTEx_KDNY_FIL,sep=""))
gtex.kdny.df <- rse_gene %>% SummarizedExperiment::assay() %>% as.data.frame()

gtex.df <- gtex.esgs.df %>%
  dplyr::bind_cols(gtex.coln.df) %>%
  dplyr::bind_cols(gtex.livr.df) %>%
  dplyr::bind_cols(gtex.brst.df) %>%
  dplyr::bind_cols(gtex.kdny.df)

rownames(gtex.df) <- rownames(gtex.kdny.df) # all gtex dataframe rows match

gtex.ensembl_gene_map <- getBM(attributes = c("ensembl_gene_id_version","external_gene_name"),
                               filters="ensembl_gene_id_version",
                               values=row.names(gtex.df),
                               mart=ensembl_dataset)

shared.ensembl.gene.ids <- intersect(archs4.ensembl_gene_map$ensembl_gene_id_version,
                                     intersect(tcga.ensembl_gene_map$ensembl_gene_id_version,
                                               gtex.ensembl_gene_map$ensembl_gene_id_version))

#12560 rows shared among ensembl_gene_id_version columns of archs4 & tcga/gtex

archs4.gene.ids <- archs4.ensembl_gene_map %>%
  dplyr::filter(ensembl_gene_id_version %in% shared.ensembl.gene.ids)

archs4.filtered.df <- archs4.df[archs4.gene.ids$external_gene_name,]
rownames(archs4.filtered.df) <- archs4.gene.ids$ensembl_gene_id_version

archs4.filtered.df %>% class()
archs4.filtered.df %>% dim()
#[1] 12560  4573

tcga.filtered.df <- tcga.df[shared.ensembl.gene.ids,]

tcga.filtered.df %>% class()
tcga.filtered.df %>% dim()
#[1] 12560   355

gtex.filtered.df <- gtex.df[shared.ensembl.gene.ids,]

gtex.filtered.df %>% class()
gtex.filtered.df %>% dim()
#[1] 12560  1556

combined.df <- archs4.filtered.df %>%
  dplyr::bind_cols(tcga.filtered.df) %>%
  dplyr::bind_cols(gtex.filtered.df)

combined.df %>% class()
combined.df %>% dim()
#[1] 12560  6484

saveRDS(combined.df,file="~/combined_df.Rds")

datasource.vec <- c(rep("ARCHS4",ncol(archs4.df)),
                    rep("TCGA",ncol(tcga.df)),
                    rep("GTEx",ncol(gtex.df)))

saveRDS(datasource.vec,file="~/datasource_vec.Rds")

archs4.tissue.vec <- c(rep("ESGS",ncol(archs4.esgs.df)),
                       rep("COLN",ncol(archs4.coln.df)),
                       rep("LIVR",ncol(archs4.livr.df)),
                       rep("BRST",ncol(archs4.brst.df)),
                       rep("KDNY",ncol(archs4.kdny.df)))

tcga.tissue.vec <- c(rep("ESGS",ncol(tcga.esgs.df)),
                     rep("COLN",ncol(tcga.coln.df)),
                     rep("LIVR",ncol(tcga.livr.df)),
                     rep("BRST",ncol(tcga.brst.df)),
                     rep("KDNY",ncol(tcga.kdny.df)))

gtex.tissue.vec <- c(rep("ESGS",ncol(gtex.esgs.df)),
                     rep("COLN",ncol(gtex.coln.df)),
                     rep("LIVR",ncol(gtex.livr.df)),
                     rep("BRST",ncol(gtex.brst.df)),
                     rep("KDNY",ncol(gtex.kdny.df)))

tissue.vec <- c(archs4.tissue.vec,
                tcga.tissue.vec,
                gtex.tissue.vec)

saveRDS(tissue.vec,file="~/tissue_vec.Rds")

colFun <- colorRampPalette(RColorBrewer::brewer.pal(10,"Paired"))

annotation_col.df <- data.frame(Datasource = datasource.vec,
                                Tissue = tissue.vec)
rownames(annotation_col.df) <- colnames(combined.df)

Datasource <- colFun(length(unique(datasource.vec)))
Tissue <- colFun(length(unique(tissue.vec)))

names(Datasource) <- annotation_col.df$Datasource %>% unique()
names(Tissue) <- annotation_col.df$Tissue %>% unique()

anno_colors <- list(Datasource = Datasource,
                    Tissue = Tissue)

# combined.df columns ordered by datasource
pheatmap::pheatmap(mat = combined.df,
                   show_rownames = FALSE,
                   show_colnames = FALSE,
                   scale = "row",
                   color = viridis::inferno(100),
                   annotation_col = annotation_col.df,
                   annotation_colors = anno_colors,
                   cluster_cols = FALSE, # clustering takes a long time
                   cluster_rows = FALSE, # clustering takes a long time
                   silent = TRUE,
                   filename = "~/combined_pheatmap_datasource.png",
                   width=8,
                   height=8,
                   fontsize = 8
)

reorder.df <- data.frame(Sample=colnames(combined.df),
                         Datasource=datasource.vec,
                         Tissue=tissue.vec)
ordered.tissues.df <- reorder.df %>%
  dplyr::arrange(Tissue)

combined.df <- combined.df[,ordered.tissues.df$Sample]

# combined.df columns ordered by tissue
pheatmap::pheatmap(mat = combined.df,
                   show_rownames = FALSE,
                   show_colnames = FALSE,
                   scale = "row",
                   color = viridis::inferno(100),
                   annotation_col = annotation_col.df,
                   annotation_colors = anno_colors,
                   cluster_cols = FALSE, # clustering takes a long time
                   cluster_rows = FALSE, # clustering takes a long time
                   silent = TRUE,
                   filename = "~/combined_pheatmap_tissue.png",
                   width=8,
                   height=8,
                   fontsize = 8
)
