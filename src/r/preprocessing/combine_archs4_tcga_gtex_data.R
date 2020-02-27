library(dplyr)
library(biomaRt)
library(magrittr)
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

archs4.df <- archs4.df[rowSums(archs4.df) > 0, colSums(archs4.df) > 0] %>% na.omit(.)

archs4.zero.na.dropped.transcripts <- setdiff(rownames(archs4.kdny.df),rownames(archs4.df))

archs4.ensembl_gene_map <- getBM(attributes = c("ensembl_gene_id_version","external_gene_name"),
                                 filters="external_gene_name",
                                 values=row.names(archs4.df),
                                 mart=ensembl_dataset)

archs4.map.dropped.transcripts <- setdiff(rownames(archs4.df),archs4.ensembl_gene_map$external_gene_name)

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

tcga.df <- tcga.df[rowSums(tcga.df) > 0,colSums(tcga.df) > 0] %>% na.omit(.)

tcga.zero.na.dropped.transcripts <- setdiff(rownames(tcga.kdny.df),rownames(tcga.df))

tcga.ensembl_gene_map <- getBM(attributes = c("ensembl_gene_id_version","external_gene_name"),
                               filters="ensembl_gene_id_version",
                               values=row.names(tcga.df),
                               mart=ensembl_dataset)

tcga.map.dropped.transcripts <- setdiff(rownames(tcga.df),tcga.ensembl_gene_map$ensembl_gene_id_version)

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

gtex.df <- gtex.df[rowSums(gtex.df) > 0, colSums(gtex.df) > 0] %>% na.omit(.)

gtex.zero.na.dropped.transcripts <- setdiff(rownames(gtex.kdny.df),rownames(gtex.df))

gtex.ensembl_gene_map <- getBM(attributes = c("ensembl_gene_id_version","external_gene_name"),
                               filters="ensembl_gene_id_version",
                               values=row.names(gtex.df),
                               mart=ensembl_dataset)

gtex.map.dropped.transcripts <- setdiff(rownames(gtex.df),gtex.ensembl_gene_map$ensembl_gene_id_version)

union.ensembl.gene.ids <- union(archs4.ensembl_gene_map$ensembl_gene_id_version,
                                union(tcga.ensembl_gene_map$ensembl_gene_id_version,
                                      gtex.ensembl_gene_map$ensembl_gene_id_version))

shared.ensembl.gene.ids <- intersect(archs4.ensembl_gene_map$ensembl_gene_id_version,
                                     intersect(tcga.ensembl_gene_map$ensembl_gene_id_version,
                                               gtex.ensembl_gene_map$ensembl_gene_id_version))

ensembl.mapped.dropped.transcript.df <- data.frame(ARCHS4=as.numeric(union.ensembl.gene.ids %in% archs4.ensembl_gene_map$ensembl_gene_id_version),
                                                      TCGA=as.numeric(union.ensembl.gene.ids %in% tcga.ensembl_gene_map$ensembl_gene_id_version),
                                                      GTEX=as.numeric(union.ensembl.gene.ids %in% gtex.ensembl_gene_map$ensembl_gene_id_version),
                                                   row.names = union.ensembl.gene.ids)

ensembl.mapped.dropped.transcript.df <- ensembl.mapped.dropped.transcript.df[rowSums(ensembl.mapped.dropped.transcript.df) < 3,]

write.csv(ensembl.mapped.dropped.transcript.df,
          file="~/dropped_transcripts_df.csv")

write.table(rownames(ensembl.mapped.dropped.transcript.df) %>% gsub("\\..+","",.),
            file="~/dropped_transcripts.csv",
            row.names=FALSE,
            col.names = FALSE,
            quote=FALSE)
write.table(archs4.zero.na.dropped.transcripts %>% gsub("\\..+","",.),
            file="~/archs4_zero_na_dropped_transcripts.csv",
            row.names = FALSE,
            col.names = FALSE,
            quote=FALSE)
write.table(archs4.map.dropped.transcripts %>% gsub("\\..+","",.),
            file="~/archs4_map_dropped_transcripts.csv",
            row.names=FALSE,
            col.names=FALSE,
            quote=FALSE)
write.table(tcga.zero.na.dropped.transcripts %>% gsub("\\..+","",.),
            file="~/tcga_zero_na_dropped_transcripts.csv",
            row.names = FALSE,
            col.names = FALSE,
            quote=FALSE)
write.table(tcga.map.dropped.transcripts %>% gsub("\\..+","",.),
            file="~/tcga_map_dropped_transcripts.csv",
            row.names=FALSE,
            col.names=FALSE,
            quote=FALSE)
write.table(gtex.zero.na.dropped.transcripts %>% gsub("\\..+","",.),
            file="~/gtex_zero_na_dropped_transcripts.csv",
            row.names = FALSE,
            col.names = FALSE,
            quote=FALSE)
write.table(gtex.map.dropped.transcripts %>% gsub("\\..+","",.),
            file="~/gtex_map_dropped_transcripts.csv",
            row.names=FALSE,
            col.names=FALSE,
            quote=FALSE)

archs4.gene.ids <- archs4.ensembl_gene_map %>%
  dplyr::filter(ensembl_gene_id_version %in% shared.ensembl.gene.ids)

archs4.filtered.df <- archs4.df[archs4.gene.ids$external_gene_name,]
rownames(archs4.filtered.df) <- archs4.gene.ids$ensembl_gene_id_version

archs4.filtered.df %>% class()
archs4.filtered.df %>% dim()

tcga.filtered.df <- tcga.df[shared.ensembl.gene.ids,]

tcga.filtered.df %>% class()
tcga.filtered.df %>% dim()

gtex.filtered.df <- gtex.df[shared.ensembl.gene.ids,]

gtex.filtered.df %>% class()
gtex.filtered.df %>% dim()

combined.df <- archs4.filtered.df %>%
  dplyr::bind_cols(tcga.filtered.df) %>%
  dplyr::bind_cols(gtex.filtered.df)

rownames(combined.df) <- rownames(archs4.filtered.df)

combined.df <- combined.df[rowSums(is.na(combined.df)) == 0,]

combined.df %>% class()
combined.df %>% dim()

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
