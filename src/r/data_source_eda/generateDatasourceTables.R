library(dplyr)
library(biomaRt)
library(magrittr)
library(SummarizedExperiment)

ensembl_dataset <- useEnsembl(biomart="ensembl",
                              dataset = "hsapiens_gene_ensembl",
                              mirror="uswest")

# load gtex data and make individual tables
#GTEx_DATA_DIR <- "/home/users/burkhajo/WuLab/WuLabLustreDir/reticula/input/recount2/recount2edGTEx/"
GTEx_DATA_DIR <- "/Users/burkhajo/Software/reticula/data/aim1/input/recount2/recount2edGTEx/"
GTEx_adipose_tissue_FIL <- "rse_gene_adipose_tissue.Rdata"
GTEx_adrenal_gland_FIL <- "rse_gene_adrenal_gland.Rdata"
GTEx_bladder_FIL <- "rse_gene_bladder.Rdata"
GTEx_blood_vessel_FIL <- "rse_gene_blood_vessel.Rdata"
GTEx_blood_FIL <- "rse_gene_blood.Rdata"
GTEx_bone_marrow_FIL <- "rse_gene_bone_marrow.Rdata"
GTEx_brain_FIL <- "rse_gene_brain.Rdata"
GTEx_breast_FIL <- "rse_gene_breast.Rdata"
GTEx_cervix_uteri_FIL <- "rse_gene_cervix_uteri.Rdata"
GTEx_colon_FIL <- "rse_gene_colon.Rdata"
GTEx_esophagus_FIL <- "rse_gene_esophagus.Rdata"
GTEx_fallopian_tube_FIL <- "rse_gene_fallopian_tube.Rdata"
GTEx_heart_FIL <- "rse_gene_heart.Rdata"
GTEx_kidney_FIL <- "rse_gene_kidney.Rdata"
GTEx_liver_FIL <- "rse_gene_liver.Rdata"
GTEx_lung_FIL <- "rse_gene_lung.Rdata"
GTEx_muscle_FIL <- "rse_gene_muscle.Rdata"
GTEx_nerve_FIL <- "rse_gene_nerve.Rdata"
GTEx_ovary_FIL <- "rse_gene_ovary.Rdata"
GTEx_pancreas_FIL <- "rse_gene_pancreas.Rdata"
GTEx_pituitary_FIL <- "rse_gene_pituitary.Rdata"
GTEx_prostate_FIL <- "rse_gene_prostate.Rdata"
GTEx_salivary_gland_FIL <- "rse_gene_salivary_gland.Rdata"
GTEx_skin_FIL <- "rse_gene_skin.Rdata"
GTEx_small_intestine_FIL <- "rse_gene_small_intestine.Rdata"
GTEx_spleen_FIL <- "rse_gene_spleen.Rdata"
GTEx_stomach_FIL <- "rse_gene_stomach.Rdata"
GTEx_testis_FIL <- "rse_gene_testis.Rdata"
GTEx_thyroid_FIL <- "rse_gene_thyroid.Rdata"
GTEx_uterus_FIL <- "rse_gene_uterus.Rdata"
GTEx_vagina_FIL <- "rse_gene_vagina.Rdata"

load(paste(GTEx_DATA_DIR,GTEx_adipose_tissue_FIL,sep=""))
gtex.adipose_tissue.cols <- rse_gene %>% colData()
gtex.adipose_tissue.df <- rse_gene %>% SummarizedExperiment::assay() %>% as.data.frame()

load(paste(GTEx_DATA_DIR,GTEx_adrenal_gland_FIL,sep=""))
gtex.adrenal_gland.cols <- rse_gene %>% colData()
gtex.adrenal_gland.df <- rse_gene %>% SummarizedExperiment::assay() %>% as.data.frame()

load(paste(GTEx_DATA_DIR,GTEx_bladder_FIL,sep=""))
gtex.bladder.cols <- rse_gene %>% colData()
gtex.bladder.df <- rse_gene %>% SummarizedExperiment::assay() %>% as.data.frame()

load(paste(GTEx_DATA_DIR,GTEx_blood_vessel_FIL,sep=""))
gtex.blood_vessel.cols <- rse_gene %>% colData()
gtex.blood_vessel.df <- rse_gene %>% SummarizedExperiment::assay() %>% as.data.frame()

load(paste(GTEx_DATA_DIR,GTEx_blood_FIL,sep=""))
gtex.blood.cols <- rse_gene %>% colData()
gtex.blood.df <- rse_gene %>% SummarizedExperiment::assay() %>% as.data.frame()

load(paste(GTEx_DATA_DIR,GTEx_bone_marrow_FIL,sep=""))
gtex.bone_marrow.cols <- rse_gene %>% colData()
gtex.bone_marrow.df <- rse_gene %>% SummarizedExperiment::assay() %>% as.data.frame()

load(paste(GTEx_DATA_DIR,GTEx_brain_FIL,sep=""))
gtex.brain.cols <- rse_gene %>% colData()
gtex.brain.df <- rse_gene %>% SummarizedExperiment::assay() %>% as.data.frame()

load(paste(GTEx_DATA_DIR,GTEx_breast_FIL,sep=""))
gtex.breast.cols <- rse_gene %>% colData()
gtex.breast.df <- rse_gene %>% SummarizedExperiment::assay() %>% as.data.frame()

load(paste(GTEx_DATA_DIR,GTEx_cervix_uteri_FIL,sep=""))
gtex.cervix_uteri.cols <- rse_gene %>% colData()
gtex.cervix_uteri.df <- rse_gene %>% SummarizedExperiment::assay() %>% as.data.frame()

load(paste(GTEx_DATA_DIR,GTEx_colon_FIL,sep=""))
gtex.colon.cols <- rse_gene %>% colData()
gtex.colon.df <- rse_gene %>% SummarizedExperiment::assay() %>% as.data.frame()

load(paste(GTEx_DATA_DIR,GTEx_esophagus_FIL,sep=""))
gtex.esophagus.cols <- rse_gene %>% colData()
gtex.esophagus.df <- rse_gene %>% SummarizedExperiment::assay() %>% as.data.frame()

load(paste(GTEx_DATA_DIR,GTEx_fallopian_tube_FIL,sep=""))
gtex.fallopian_tube.cols <- rse_gene %>% colData()
gtex.fallopian_tube.df <- rse_gene %>% SummarizedExperiment::assay() %>% as.data.frame()

load(paste(GTEx_DATA_DIR,GTEx_heart_FIL,sep=""))
gtex.heart.cols <- rse_gene %>% colData()
gtex.heart.df <- rse_gene %>% SummarizedExperiment::assay() %>% as.data.frame()

load(paste(GTEx_DATA_DIR,GTEx_kidney_FIL,sep=""))
gtex.kidney.cols <- rse_gene %>% colData()
gtex.kidney.df <- rse_gene %>% SummarizedExperiment::assay() %>% as.data.frame()

load(paste(GTEx_DATA_DIR,GTEx_liver_FIL,sep=""))
gtex.liver.cols <- rse_gene %>% colData()
gtex.liver.df <- rse_gene %>% SummarizedExperiment::assay() %>% as.data.frame()

load(paste(GTEx_DATA_DIR,GTEx_lung_FIL,sep=""))
gtex.lung.cols <- rse_gene %>% colData()
gtex.lung.df <- rse_gene %>% SummarizedExperiment::assay() %>% as.data.frame()

load(paste(GTEx_DATA_DIR,GTEx_muscle_FIL,sep=""))
gtex.muscle.cols <- rse_gene %>% colData()
gtex.muscle.df <- rse_gene %>% SummarizedExperiment::assay() %>% as.data.frame()

load(paste(GTEx_DATA_DIR,GTEx_nerve_FIL,sep=""))
gtex.nerve.cols <- rse_gene %>% colData()
gtex.nerve.df <- rse_gene %>% SummarizedExperiment::assay() %>% as.data.frame()

load(paste(GTEx_DATA_DIR,GTEx_ovary_FIL,sep=""))
gtex.ovary.cols <- rse_gene %>% colData()
gtex.ovary.df <- rse_gene %>% SummarizedExperiment::assay() %>% as.data.frame()

load(paste(GTEx_DATA_DIR,GTEx_pancreas_FIL,sep=""))
gtex.pancreas.cols <- rse_gene %>% colData()
gtex.pancreas.df <- rse_gene %>% SummarizedExperiment::assay() %>% as.data.frame()

load(paste(GTEx_DATA_DIR,GTEx_pituitary_FIL,sep=""))
gtex.pituitary.cols <- rse_gene %>% colData()
gtex.pituitary.df <- rse_gene %>% SummarizedExperiment::assay() %>% as.data.frame()

load(paste(GTEx_DATA_DIR,GTEx_prostate_FIL,sep=""))
gtex.prostate.cols <- rse_gene %>% colData()
gtex.prostate.df <- rse_gene %>% SummarizedExperiment::assay() %>% as.data.frame()

load(paste(GTEx_DATA_DIR,GTEx_salivary_gland_FIL,sep=""))
gtex.salivary_gland.cols <- rse_gene %>% colData()
gtex.salivary_gland.df <- rse_gene %>% SummarizedExperiment::assay() %>% as.data.frame()

load(paste(GTEx_DATA_DIR,GTEx_skin_FIL,sep=""))
gtex.skin.cols <- rse_gene %>% colData()
gtex.skin.df <- rse_gene %>% SummarizedExperiment::assay() %>% as.data.frame()

load(paste(GTEx_DATA_DIR,GTEx_small_intestine_FIL,sep=""))
gtex.small_intestine.cols <- rse_gene %>% colData()
gtex.small_intestine.df <- rse_gene %>% SummarizedExperiment::assay() %>% as.data.frame()

load(paste(GTEx_DATA_DIR,GTEx_spleen_FIL,sep=""))
gtex.spleen.cols <- rse_gene %>% colData()
gtex.spleen.df <- rse_gene %>% SummarizedExperiment::assay() %>% as.data.frame()

load(paste(GTEx_DATA_DIR,GTEx_stomach_FIL,sep=""))
gtex.stomach.cols <- rse_gene %>% colData()
gtex.stomach.df <- rse_gene %>% SummarizedExperiment::assay() %>% as.data.frame()

load(paste(GTEx_DATA_DIR,GTEx_testis_FIL,sep=""))
gtex.testis.cols <- rse_gene %>% colData()
gtex.testis.df <- rse_gene %>% SummarizedExperiment::assay() %>% as.data.frame()

load(paste(GTEx_DATA_DIR,GTEx_thyroid_FIL,sep=""))
gtex.thyroid.cols <- rse_gene %>% colData()
gtex.thyroid.df <- rse_gene %>% SummarizedExperiment::assay() %>% as.data.frame()

load(paste(GTEx_DATA_DIR,GTEx_uterus_FIL,sep=""))
gtex.uterus.cols <- rse_gene %>% colData()
gtex.uterus.df <- rse_gene %>% SummarizedExperiment::assay() %>% as.data.frame()

load(paste(GTEx_DATA_DIR,GTEx_vagina_FIL,sep=""))
gtex.vagina.cols <- rse_gene %>% colData()
gtex.vagina.df <- rse_gene %>% SummarizedExperiment::assay() %>% as.data.frame()

gtex.tissue.vec <- c(gtex.adipose_tissue.cols$smts,
                     gtex.adrenal_gland.cols$smts,
                     gtex.bladder.cols$smts,
                     gtex.blood_vessel.cols$smts,
                     gtex.blood.cols$smts,
                     gtex.bone_marrow.cols$smts,
                     gtex.brain.cols$smts,
                     gtex.breast.cols$smts,
                     gtex.cervix_uteri.cols$smts,
                     gtex.colon.cols$smts,
                     gtex.esophagus.cols$smts,
                     gtex.fallopian_tube.cols$smts,
                     gtex.heart.cols$smts,
                     gtex.kidney.cols$smts,
                     gtex.liver.cols$smts,
                     gtex.lung.cols$smts,
                     gtex.muscle.cols$smts,
                     gtex.nerve.cols$smts,
                     gtex.ovary.cols$smts,
                     gtex.pancreas.cols$smts,
                     gtex.pituitary.cols$smts,
                     gtex.prostate.cols$smts,
                     gtex.salivary_gland.cols$smts,
                     gtex.skin.cols$smts,
                     gtex.small_intestine.cols$smts,
                     gtex.spleen.cols$smts,
                     gtex.stomach.cols$smts,
                     gtex.testis.cols$smts,
                     gtex.thyroid.cols$smts,
                     gtex.uterus.cols$smts,
                     gtex.vagina.cols$smts)

gtex.tissue.vec %>% table()

gtex.tissue.detail.vec <- c(gtex.adipose_tissue.cols$smtsd,
                            gtex.adrenal_gland.cols$smtsd,
                            gtex.bladder.cols$smtsd,
                            gtex.blood_vessel.cols$smtsd,
                            gtex.blood.cols$smtsd,
                            gtex.bone_marrow.cols$smtsd,
                            gtex.brain.cols$smtsd,
                            gtex.breast.cols$smtsd,
                            gtex.cervix_uteri.cols$smtsd,
                            gtex.colon.cols$smtsd,
                            gtex.esophagus.cols$smtsd,
                            gtex.fallopian_tube.cols$smtsd,
                            gtex.heart.cols$smtsd,
                            gtex.kidney.cols$smtsd,
                            gtex.liver.cols$smtsd,
                            gtex.lung.cols$smtsd,
                            gtex.muscle.cols$smtsd,
                            gtex.nerve.cols$smtsd,
                            gtex.ovary.cols$smtsd,
                            gtex.pancreas.cols$smtsd,
                            gtex.pituitary.cols$smtsd,
                            gtex.prostate.cols$smtsd,
                            gtex.salivary_gland.cols$smtsd,
                            gtex.skin.cols$smtsd,
                            gtex.small_intestine.cols$smtsd,
                            gtex.spleen.cols$smtsd,
                            gtex.stomach.cols$smtsd,
                            gtex.testis.cols$smtsd,
                            gtex.thyroid.cols$smtsd,
                            gtex.uterus.cols$smtsd,
                            gtex.vagina.cols$smtsd)

gtex.tissue.detail.vec %>% table()

gtex.df <- gtex.adipose_tissue.df %>%
  dplyr::bind_cols(gtex.adrenal_gland.df) %>%
  dplyr::bind_cols(gtex.bladder.df) %>%
  dplyr::bind_cols(gtex.blood_vessel.df) %>%
  dplyr::bind_cols(gtex.blood.df) %>%
  dplyr::bind_cols(gtex.bone_marrow.df) %>%
  dplyr::bind_cols(gtex.brain.df) %>%
  dplyr::bind_cols(gtex.breast.df) %>%
  dplyr::bind_cols(gtex.cervix_uteri.df) %>%
  dplyr::bind_cols(gtex.colon.df) %>%
  dplyr::bind_cols(gtex.esophagus.df) %>%
  dplyr::bind_cols(gtex.fallopian_tube.df) %>%
  dplyr::bind_cols(gtex.heart.df) %>%
  dplyr::bind_cols(gtex.kidney.df) %>%
  dplyr::bind_cols(gtex.liver.df) %>%
  dplyr::bind_cols(gtex.lung.df) %>%
  dplyr::bind_cols(gtex.muscle.df) %>%
  dplyr::bind_cols(gtex.nerve.df) %>%
  dplyr::bind_cols(gtex.ovary.df) %>%
  dplyr::bind_cols(gtex.pancreas.df) %>%
  dplyr::bind_cols(gtex.pituitary.df) %>%
  dplyr::bind_cols(gtex.prostate.df) %>%
  dplyr::bind_cols(gtex.salivary_gland.df) %>%
  dplyr::bind_cols(gtex.skin.df) %>%
  dplyr::bind_cols(gtex.small_intestine.df) %>%
  dplyr::bind_cols(gtex.spleen.df) %>%
  dplyr::bind_cols(gtex.stomach.df) %>%
  dplyr::bind_cols(gtex.testis.df) %>%
  dplyr::bind_cols(gtex.thyroid.df) %>%
  dplyr::bind_cols(gtex.uterus.df) %>%
  dplyr::bind_cols(gtex.vagina.df)

rownames(gtex.df) <- rownames(gtex.adipose_tissue.df) #%>% gsub("\\..+","",.) # all gtex dataframe rows match, stripping ensembl version

gtex.df <- gtex.df[rowSums(gtex.df) > 0, colSums(gtex.df) > 0] %>% na.omit(.)

gtex.zero.na.dropped.transcripts <- setdiff(rownames(gtex.adipose_tissue.df),rownames(gtex.df))

gtex.ensembl_gene_map <- getBM(attributes = c("ensembl_gene_id_version","external_gene_name"),
                               filters="ensembl_gene_id_version",
                               values=row.names(gtex.df),
                               mart=ensembl_dataset)

gtex.map.dropped.transcripts <- setdiff(rownames(gtex.df),gtex.ensembl_gene_map$ensembl_gene_id_version)

# load tcga data and make combined tables
#TCGA_DATA_DIR <- "/home/users/burkhajo/WuLab/WuLabLustreDir/reticula/input/recount2/recount2edTCGA/"
TCGA_DATA_DIR <- "/Users/burkhajo/Software/reticula/data/aim1/input/recount2/recount2edTCGA/"
TCGA_adrenal_gland_FIL <- "rse_gene_adrenal_gland.Rdata"
TCGA_bladder_FIL <- "rse_gene_bladder.Rdata"
TCGA_bone_marrow_FIL <- "rse_gene_bone_marrow.Rdata"
TCGA_brain_FIL <- "rse_gene_brain.Rdata"
TCGA_breast_FIL <- "rse_gene_breast.Rdata"
TCGA_cervix_FIL <- "rse_gene_cervix.Rdata"
TCGA_colorectal_FIL <- "rse_gene_colorectal.Rdata"
TCGA_esophagus_FIL <- "rse_gene_esophagus.Rdata"
TCGA_kidney_FIL <- "rse_gene_kidney.Rdata"
TCGA_liver_FIL <- "rse_gene_liver.Rdata"
TCGA_lung_FIL <- "rse_gene_lung.Rdata"
TCGA_ovary_FIL <- "rse_gene_ovary.Rdata"
TCGA_pancreas_FIL <- "rse_gene_pancreas.Rdata"
TCGA_prostate_FIL <- "rse_gene_prostate.Rdata"
TCGA_skin_FIL <- "rse_gene_skin.Rdata"
TCGA_stomach_FIL <- "rse_gene_stomach.Rdata"
TCGA_testis_FIL <- "rse_gene_testis.Rdata"
TCGA_thyroid_FIL <- "rse_gene_thyroid.Rdata"
TCGA_uterus_FIL <- "rse_gene_uterus.Rdata"


load(paste(TCGA_DATA_DIR,TCGA_adrenal_gland_FIL,sep=""))
tcga.adrenal_gland.cols <- rse_gene %>% colData()
tcga.adrenal_gland.df <- rse_gene %>% SummarizedExperiment::assay() %>% .[,tcga.adrenal_gland.cols$gdc_cases.samples.sample_type == "Solid Tissue Normal"] %>% as.data.frame()

load(paste(TCGA_DATA_DIR,TCGA_bladder_FIL,sep=""))
tcga.bladder.cols <- rse_gene %>% colData()
tcga.bladder.df <- rse_gene %>% SummarizedExperiment::assay() %>% .[,tcga.bladder.cols$gdc_cases.samples.sample_type == "Solid Tissue Normal"] %>% as.data.frame()

load(paste(TCGA_DATA_DIR,TCGA_bone_marrow_FIL,sep=""))
tcga.bone_marrow.cols <- rse_gene %>% colData()
tcga.bone_marrow.df <- rse_gene %>% SummarizedExperiment::assay() %>% .[,tcga.bone_marrow.cols$gdc_cases.samples.sample_type == "Solid Tissue Normal"] %>% as.data.frame()

load(paste(TCGA_DATA_DIR,TCGA_brain_FIL,sep=""))
tcga.brain.cols <- rse_gene %>% colData()
tcga.brain.df <- rse_gene %>% SummarizedExperiment::assay() %>% .[,tcga.brain.cols$gdc_cases.samples.sample_type == "Solid Tissue Normal"] %>% as.data.frame()

load(paste(TCGA_DATA_DIR,TCGA_breast_FIL,sep=""))
tcga.breast.cols <- rse_gene %>% colData()
tcga.breast.df <- rse_gene %>% SummarizedExperiment::assay() %>% .[,tcga.breast.cols$gdc_cases.samples.sample_type == "Solid Tissue Normal"] %>% as.data.frame()

load(paste(TCGA_DATA_DIR,TCGA_cervix_FIL,sep=""))
tcga.cervix.cols <- rse_gene %>% colData()
tcga.cervix.df <- rse_gene %>% SummarizedExperiment::assay() %>% .[,tcga.cervix.cols$gdc_cases.samples.sample_type == "Solid Tissue Normal"] %>% as.data.frame()

load(paste(TCGA_DATA_DIR,TCGA_colorectal_FIL,sep=""))
tcga.colorectal.cols <- rse_gene %>% colData()
tcga.colorectal.df <- rse_gene %>% SummarizedExperiment::assay() %>% .[,tcga.colorectal.cols$gdc_cases.samples.sample_type == "Solid Tissue Normal"] %>% as.data.frame()

load(paste(TCGA_DATA_DIR,TCGA_esophagus_FIL,sep=""))
tcga.esophagus.cols <- rse_gene %>% colData()
tcga.esophagus.df <- rse_gene %>% SummarizedExperiment::assay() %>% .[,tcga.esophagus.cols$gdc_cases.samples.sample_type == "Solid Tissue Normal"] %>% as.data.frame()

load(paste(TCGA_DATA_DIR,TCGA_kidney_FIL,sep=""))
tcga.kidney.cols <- rse_gene %>% colData()
tcga.kidney.df <- rse_gene %>% SummarizedExperiment::assay() %>% .[,tcga.kidney.cols$gdc_cases.samples.sample_type == "Solid Tissue Normal"] %>% as.data.frame()

load(paste(TCGA_DATA_DIR,TCGA_liver_FIL,sep=""))
tcga.liver.cols <- rse_gene %>% colData()
tcga.liver.df <- rse_gene %>% SummarizedExperiment::assay() %>% .[,tcga.liver.cols$gdc_cases.samples.sample_type == "Solid Tissue Normal"] %>% as.data.frame()

load(paste(TCGA_DATA_DIR,TCGA_lung_FIL,sep=""))
tcga.lung.cols <- rse_gene %>% colData()
tcga.lung.df <- rse_gene %>% SummarizedExperiment::assay() %>% .[,tcga.lung.cols$gdc_cases.samples.sample_type == "Solid Tissue Normal"] %>% as.data.frame()

load(paste(TCGA_DATA_DIR,TCGA_ovary_FIL,sep=""))
tcga.ovary.cols <- rse_gene %>% colData()
tcga.ovary.df <- rse_gene %>% SummarizedExperiment::assay() %>% .[,tcga.ovary.cols$gdc_cases.samples.sample_type == "Solid Tissue Normal"] %>% as.data.frame()

load(paste(TCGA_DATA_DIR,TCGA_pancreas_FIL,sep=""))
tcga.pancreas.cols <- rse_gene %>% colData()
tcga.pancreas.df <- rse_gene %>% SummarizedExperiment::assay() %>% .[,tcga.pancreas.cols$gdc_cases.samples.sample_type == "Solid Tissue Normal"] %>% as.data.frame()

load(paste(TCGA_DATA_DIR,TCGA_prostate_FIL,sep=""))
tcga.prostate.cols <- rse_gene %>% colData()
tcga.prostate.df <- rse_gene %>% SummarizedExperiment::assay() %>% .[,tcga.prostate.cols$gdc_cases.samples.sample_type == "Solid Tissue Normal"] %>% as.data.frame()

load(paste(TCGA_DATA_DIR,TCGA_skin_FIL,sep=""))
tcga.skin.cols <- rse_gene %>% colData()
tcga.skin.df <- rse_gene %>% SummarizedExperiment::assay() %>% .[,tcga.skin.cols$gdc_cases.samples.sample_type == "Solid Tissue Normal"] %>% as.data.frame()

load(paste(TCGA_DATA_DIR,TCGA_stomach_FIL,sep=""))
tcga.stomach.cols <- rse_gene %>% colData()
tcga.stomach.df <- rse_gene %>% SummarizedExperiment::assay() %>% .[,tcga.stomach.cols$gdc_cases.samples.sample_type == "Solid Tissue Normal"] %>% as.data.frame()

load(paste(TCGA_DATA_DIR,TCGA_testis_FIL,sep=""))
tcga.testis.cols <- rse_gene %>% colData()
tcga.testis.df <- rse_gene %>% SummarizedExperiment::assay() %>% .[,tcga.testis.cols$gdc_cases.samples.sample_type == "Solid Tissue Normal"] %>% as.data.frame()

load(paste(TCGA_DATA_DIR,TCGA_thyroid_FIL,sep=""))
tcga.thyroid.cols <- rse_gene %>% colData()
tcga.thyroid.df <- rse_gene %>% SummarizedExperiment::assay() %>% .[,tcga.thyroid.cols$gdc_cases.samples.sample_type == "Solid Tissue Normal"] %>% as.data.frame()

load(paste(TCGA_DATA_DIR,TCGA_uterus_FIL,sep=""))
tcga.uterus.cols <- rse_gene %>% colData()
tcga.uterus.df <- rse_gene %>% SummarizedExperiment::assay() %>% .[,tcga.uterus.cols$gdc_cases.samples.sample_type == "Solid Tissue Normal"] %>% as.data.frame()

tcga.tissue.vec <- c(tcga.adrenal_gland.cols[tcga.adrenal_gland.cols$gdc_cases.samples.sample_type == "Solid Tissue Normal",]$gdc_cases.project.primary_site,
                     tcga.bladder.cols[tcga.bladder.cols$gdc_cases.samples.sample_type == "Solid Tissue Normal",]$gdc_cases.project.primary_site,
                     tcga.bone_marrow.cols[tcga.bone_marrow.cols$gdc_cases.samples.sample_type == "Solid Tissue Normal",]$gdc_cases.project.primary_site,
                     tcga.brain.cols[tcga.brain.cols$gdc_cases.samples.sample_type == "Solid Tissue Normal",]$gdc_cases.project.primary_site,
                     tcga.breast.cols[tcga.breast.cols$gdc_cases.samples.sample_type == "Solid Tissue Normal",]$gdc_cases.project.primary_site,
                     tcga.cervix.cols[tcga.cervix.cols$gdc_cases.samples.sample_type == "Solid Tissue Normal",]$gdc_cases.project.primary_site,
                     tcga.colorectal.cols[tcga.colorectal.cols$gdc_cases.samples.sample_type == "Solid Tissue Normal",]$gdc_cases.project.primary_site,
                     tcga.esophagus.cols[tcga.esophagus.cols$gdc_cases.samples.sample_type == "Solid Tissue Normal",]$gdc_cases.project.primary_site,
                     tcga.kidney.cols[tcga.kidney.cols$gdc_cases.samples.sample_type == "Solid Tissue Normal",]$gdc_cases.project.primary_site,
                     tcga.liver.cols[tcga.liver.cols$gdc_cases.samples.sample_type == "Solid Tissue Normal",]$gdc_cases.project.primary_site,
                     tcga.lung.cols[tcga.lung.cols$gdc_cases.samples.sample_type == "Solid Tissue Normal",]$gdc_cases.project.primary_site,
                     tcga.ovary.cols[tcga.ovary.cols$gdc_cases.samples.sample_type == "Solid Tissue Normal",]$gdc_cases.project.primary_site,
                     tcga.pancreas.cols[tcga.pancreas.cols$gdc_cases.samples.sample_type == "Solid Tissue Normal",]$gdc_cases.project.primary_site,
                     tcga.prostate.cols[tcga.prostate.cols$gdc_cases.samples.sample_type == "Solid Tissue Normal",]$gdc_cases.project.primary_site,
                     tcga.skin.cols[tcga.skin.cols$gdc_cases.samples.sample_type == "Solid Tissue Normal",]$gdc_cases.project.primary_site,
                     tcga.stomach.cols[tcga.stomach.cols$gdc_cases.samples.sample_type == "Solid Tissue Normal",]$gdc_cases.project.primary_site,
                     tcga.testis.cols[tcga.testis.cols$gdc_cases.samples.sample_type == "Solid Tissue Normal",]$gdc_cases.project.primary_site,
                     tcga.thyroid.cols[tcga.thyroid.cols$gdc_cases.samples.sample_type == "Solid Tissue Normal",]$gdc_cases.project.primary_site,
                     tcga.uterus.cols[tcga.uterus.cols$gdc_cases.samples.sample_type == "Solid Tissue Normal",]$gdc_cases.project.primary_site)

tcga.tissue.vec %>% table()

tcga.df <- tcga.adrenal_gland.df %>%
  dplyr::bind_cols(tcga.bladder.df) %>%
  dplyr::bind_cols(tcga.bone_marrow.df) %>%
  dplyr::bind_cols(tcga.brain.df) %>%
  dplyr::bind_cols(tcga.breast.df) %>%
  dplyr::bind_cols(tcga.cervix.df) %>%
  dplyr::bind_cols(tcga.colorectal.df) %>%
  dplyr::bind_cols(tcga.esophagus.df) %>%
  dplyr::bind_cols(tcga.kidney.df) %>%
  dplyr::bind_cols(tcga.liver.df) %>%
  dplyr::bind_cols(tcga.lung.df) %>%
  dplyr::bind_cols(tcga.ovary.df) %>%
  dplyr::bind_cols(tcga.pancreas.df) %>%
  dplyr::bind_cols(tcga.prostate.df) %>%
  dplyr::bind_cols(tcga.skin.df) %>%
  dplyr::bind_cols(tcga.stomach.df) %>%
  dplyr::bind_cols(tcga.testis.df) %>%
  dplyr::bind_cols(tcga.thyroid.df) %>%
  dplyr::bind_cols(tcga.uterus.df)

rownames(tcga.df) <- rownames(tcga.adrenal_gland.df)  #%>% gsub("\\..+","",.) # all tcga dataframe rows match, stripping ensembl version

tcga.df <- tcga.df[rowSums(tcga.df) > 0,colSums(tcga.df) > 0] %>% na.omit(.)

tcga.zero.na.dropped.transcripts <- setdiff(rownames(tcga.adrenal_gland.df),rownames(tcga.df))

tcga.ensembl_gene_map <- getBM(attributes = c("ensembl_gene_id_version","external_gene_name"),
                               filters="ensembl_gene_id_version",
                               values=row.names(tcga.df),
                               mart=ensembl_dataset)

tcga.map.dropped.transcripts <- setdiff(rownames(tcga.df),tcga.ensembl_gene_map$ensembl_gene_id_version)


# load archs4 data an dmake combined tables
#ARCHS4_DATA_DIR <- "/home/users/burkhajo/WuLab/WuLabLustreDir/reticula/input/ARCHS4/"
ARCHS4_DATA_DIR <- "/Users/burkhajo/Software/reticula/data/aim1/input/ARCHS4/"
ARCHS4_breast_FIL <- "Breast_expression_matrix.tsv"
ARCHS4_colon_FIL <- "Colon_expression_matrix.tsv"
ARCHS4_esophagus_FIL <- "Esophagus_expression_matrix.tsv"
ARCHS4_kidney_FIL <- "Kidney_expression_matrix.tsv"
ARCHS4_liver_FIL <- "Liver_expression_matrix.tsv"

archs4.breast.df <- read.table(paste(ARCHS4_DATA_DIR,ARCHS4_breast_FIL,sep=""),header = TRUE,row.names = 1)
archs4.colon.df <- read.table(paste(ARCHS4_DATA_DIR,ARCHS4_colon_FIL,sep=""),header = TRUE,row.names = 1)
archs4.esophagus.df <- read.table(paste(ARCHS4_DATA_DIR,ARCHS4_esophagus_FIL,sep=""),header = TRUE,row.names = 1)
archs4.kidney.df <- read.table(paste(ARCHS4_DATA_DIR,ARCHS4_kidney_FIL,sep=""),header = TRUE,row.names = 1)
archs4.liver.df <- read.table(paste(ARCHS4_DATA_DIR,ARCHS4_liver_FIL,sep=""),header = TRUE,row.names = 1)

archs4.tissue.vec <- c(rep("breast",ncol(archs4.breast.df)),
                       rep("colon",ncol(archs4.colon.df)),
                       rep("esophagus",ncol(archs4.esophagus.df)),
                       rep("kidney",ncol(archs4.kidney.df)),
                       rep("liver",ncol(archs4.liver.df)))

archs4.tissue.vec %>% table()

archs4.df <- archs4.breast.df %>%
  dplyr::bind_cols(archs4.colon.df) %>%
  dplyr::bind_cols(archs4.esophagus.df) %>%
  dplyr::bind_cols(archs4.kidney.df) %>%
  dplyr::bind_cols(archs4.liver.df)

rownames(archs4.df) <- rownames(archs4.breast.df) # all archs4 dataframe rows match

archs4.df <- archs4.df[rowSums(archs4.df) > 0, colSums(archs4.df) > 0] %>% na.omit(.)

archs4.zero.na.dropped.transcripts <- setdiff(rownames(archs4.breast.df),rownames(archs4.df))

archs4.ensembl_gene_map <- getBM(attributes = c("ensembl_gene_id_version","external_gene_name"),
                                 filters="external_gene_name",
                                 values=row.names(archs4.df),
                                 mart=ensembl_dataset)

archs4.map.dropped.transcripts <- setdiff(rownames(archs4.df),archs4.ensembl_gene_map$external_gene_name)

# check map dropped transcripts
union.ensembl.gene.ids <- union(archs4.ensembl_gene_map$ensembl_gene_id_version,
                                union(gtex.ensembl_gene_map$ensembl_gene_id_version,
                                      union(tcga.ensembl_gene_map$ensembl_gene_id_version)))

shared.ensembl.gene.ids <- intersect(archs4.ensembl_gene_map$ensembl_gene_id,
                                     intersect(gtex.ensembl_gene_map$ensembl_gene_id,
                                               intersect(tcga.ensembl_gene_map$ensembl_gene_id)))

ensembl.mapped.dropped.transcript.df <- data.frame(ARCHS4=as.numeric(union.ensembl.gene.ids %in% archs4.ensembl_gene_map$ensembl_gene_id),
                                                   GTEX=as.numeric(union.ensembl.gene.ids %in% gtex.ensembl_gene_map$ensembl_gene_id),
                                                   TCGA=as.numeric(union.ensembl.gene.ids %in% tcga.ensembl_gene_map$ensembl_gene_id),
                                                   row.names = union.ensembl.gene.ids)

ensembl.mapped.dropped.transcript.df <- ensembl.mapped.dropped.transcript.df[rowSums(ensembl.mapped.dropped.transcript.df) < 3,]