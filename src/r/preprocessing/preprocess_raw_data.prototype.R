library(CePa)
library(dplyr)
library(DESeq2)
library(biomaRt)
library(ggplot2)
library(magrittr)
library(SummarizedExperiment)

cls <- function(){
  cat("\014")
}

ensembl_dataset <- useEnsembl(biomart="ensembl",
                              dataset = "hsapiens_gene_ensembl",
                              mirror="uswest")

ARCHS4_DATA_DIR <- "/Users/burkhajo/Software/reticula/data/aim1/input/ARCHS4/"
ARCHS4_DATA_FIL <- "Ovary_expression_matrix.tsv"
archs4.df <- read.table(paste(ARCHS4_DATA_DIR,ARCHS4_DATA_FIL,sep=""),header = TRUE,row.names = 1)

#archs4.df with 35238 rows

archs4.ensembl_gene_map <- getBM(attributes = c("ensembl_gene_id_version","external_gene_name"),
                          filters="external_gene_name",
                          values=row.names(archs4.df),
                          mart=ensembl_dataset)

#archs4.ensembl_gene_map with 31028 rows

TCGA_DATA_DIR <- "/Users/burkhajo/Software/reticula/data/aim1/input/recount2/recount2edTCGA/"
TCGA_DATA_FIL <- "rse_gene_breast.Rdata"
load(paste(TCGA_DATA_DIR,TCGA_DATA_FIL,sep=""))
tcga.df <- rse_gene %>% colData()
tcga.df$gdc_cases.samples.sample_type %>% as.factor() %>% summary()
fz <- rse_gene %>% assay() %>% .[,tcga.df$gdc_cases.samples.sample_type == "Solid Tissue Normal"] %>% as.data.frame()

#fz with 58037 rows

tcga.ensembl_gene_map <- getBM(attributes = c("ensembl_gene_id_version","external_gene_name"),
                          filters="ensembl_gene_id_version",
                          values=row.names(fz),
                          mart=ensembl_dataset)

#tcga.ensembl_gene_map with 38321 rows




GTEx_DATA_DIR <- "/Users/burkhajo/Software/reticula/data/aim1/input/recount2/recount2edGTEx/"
GTEx_DATA_FIL <- "rse_gene_esophagus.Rdata"
load(paste(GTEx_DATA_DIR,GTEx_DATA_FIL,sep=""))
gtex.df <- rse_gene %>% colData()
gz <- rse_gene %>% assay() %>% as.data.frame()

#gz with 58037 rows

gtex.ensembl_gene_map <- getBM(attributes = c("ensembl_gene_id_version","external_gene_name"),
                               filters="ensembl_gene_id_version",
                               values=row.names(gz),
                               mart=ensembl_dataset)

#gtex.ensembl_gene_map with 38321 rows





INPUT_DATA_DIR <- "/home/users/burkhajo/WuLab/WuLabLustreDir/reticula/input/"

GTEX_COUNTS_FN <- "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.bz2"
GTEX_S_ATTR_FN <- "GTEx_Analysis_v8_Annotations_SampleAttributesDS.csv"
GTEX_TIS_OI_FN <- "20_types_6_subtypes.csv"
GTEX_RDS_MX_FN <- "gtex_matrix.rds"
DESEQ2_NCTS_FN <- "normalized_counts.rds"
GGPLT_TISSU_FN <- "gtex_tissue_counts.png"

# read third line of the humongous GTEx counts file
sample.names <- read.table(paste(INPUT_DATA_DIR,GTEX_COUNTS_FN,sep=""),
                           nrows=1,
                           skip=2)

sample.attrs <- read.csv(paste(INPUT_DATA_DIR,GTEX_S_ATTR_FN,sep=""),
                         header = TRUE) %>%
  dplyr::select(SAMPID,
                SMTS,
                SMTSD)

tis.o.intrst <- read.csv(paste(INPUT_DATA_DIR,GTEX_TIS_OI_FN,sep=""),
                         header = TRUE)

samples.o.intrst <- sample.attrs %>%
  dplyr::filter(SMTSD %in% tis.o.intrst$Subtype)

soi <- samples.o.intrst$SAMPID %>% as.character()

# initiate colClassVec with character columns for Name and Description
colClassVec <- c("character","character")

for(i in 3:(ncol(sample.names))){
  if(mod(i,1000) == 0){
    print(paste(i," elements in colClassVec..."))
  }
  element <- "NULL"
  if(sample.names[1,i] %in% soi){
    element <- "integer"
  }
  colClassVec <- colClassVec %>% append(element)
}

count.data <- read.table(paste(INPUT_DATA_DIR,GTEX_COUNTS_FN,sep=""),
                         header=TRUE,
                         skip=2,
                         colClasses = colClassVec)
saveRDS(count.data,
        file=paste(INPUT_DATA_DIR,GTEX_RDS_MX_FN,sep=""))

#design matrix... coerce tissue subtypes into factors/integers and store as RDS files
samples.o.intrst$SAMPID <- gsub('-','.',samples.o.intrst$SAMPID)
samples.o.intrst <- samples.o.intrst %>%
  dplyr::filter(SAMPID %in% colnames(count.data))

samples.o.intrst$CLASS <- as.numeric(samples.o.intrst$SMTSD)
rownames(samples.o.intrst) <- samples.o.intrst$SAMPID

ggplot(data.frame(samples.o.intrst$SMTS),
       aes(x=samples.o.intrst$SMTS)) +
  geom_bar()
ggsave(filename = paste(INPUT_DATA_DIR,GGPLT_TISSU_FN,sep=""),
       width = 5,
       height = 5)

rownames(count.data) <- count.data$Name

count.data <- count.data %>%
  dplyr::select(rownames(samples.o.intrst))

#test DESeq2 on count data/design mtx subsets
dds <- DESeq2::DESeqDataSetFromMatrix(countData = count.data,
                                      colData = samples.o.intrst,
                                      design = ~ SMTS)
dds <- DESeq2::DESeq(dds)
dds <- estimateSizeFactors(dds)
normalized_counts <- DESeq2::counts(dds,normalized = TRUE)
saveRDS(normalized_counts,
        file=paste(INPUT_DATA_DIR,DESEQ2_NCTS_FN,sep=""))

#rewrite this file as a template and integrate into reticula workflow

