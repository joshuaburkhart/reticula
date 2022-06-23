library(multiROC)
library(magrittr)
library(ggplot2)
library(data.table)
library(stringr)
library(dplyr)

DATA_DIR <- "/home/jgburk/PycharmProjects/reticula/data/tcga/output/"
tcga_tune_test_dataset.df <- read.table(paste(DATA_DIR,"TCGA_Solid_Tissue_Normal_Samples_Dataset.csv",sep=""), sep = ",")
z <- table(tcga_tune_test_dataset.df)
z[,c(2,1)]

tcga_test_gnn_calls.df <- read.table(paste(DATA_DIR,"TCGA_Solid_Tissue_Normal_Samples_GNN_Calls.csv",sep=""),sep=",")
tcga_test_resnet_calls.df <- read.table(paste(DATA_DIR,"TCGA_Solid_Tissue_Normal_Samples_Resnet_Calls.csv",sep=""),sep=",")

tissue_name2code_gnn <- list()
tissue_code2name_gnn <- list()
for(i in 1:nrow(tcga_test_gnn_calls.df)){
  tis_name <- tcga_test_gnn_calls.df[i,1]
  tis_code <- tcga_test_gnn_calls.df[i,3]
  tissue_name2code_gnn[[tis_name]] <- tis_code
  tissue_code2name_gnn[[as.character(tis_code)]] <- tis_name
}

tissue_name2code_res <- list()
tissue_code2name_res <- list()
for(i in 1:nrow(tcga_test_resnet_calls.df)){
  tis_name <- tcga_test_resnet_calls.df[i,1]
  tis_code <- tcga_test_resnet_calls.df[i,3]
  tissue_name2code_res[[tis_name]] <- tis_code
  tissue_code2name_res[[as.character(tis_code)]] <- tis_name
}
stopifnot(assertthat::are_equal(tissue_name2code_gnn,tissue_name2code_res))
stopifnot(assertthat::are_equal(tissue_code2name_gnn,tissue_code2name_res))

tcga_test_gnn_calls.df$rownames <- rownames(tcga_test_gnn_calls.df)
tcga_test_resnet_calls.df$rownames <- rownames(tcga_test_resnet_calls.df)

stopifnot(assertthat::are_equal(tcga_test_gnn_calls.df$rownames,tcga_test_resnet_calls.df$rownames))

table(c(unique(tcga_test_gnn_calls.df$V3),unique(tcga_test_gnn_calls.df$V4))) # 5, 12, 17
table(c(unique(tcga_test_resnet_calls.df$V3),unique(tcga_test_resnet_calls.df$V4))) # 12, 15, 17

tissue_code2name_gnn[["5"]]
#[1] "Cervix"
tissue_code2name_gnn[["12"]]
#[1] "Pancreas"
tissue_code2name_gnn[["15"]]
#[1] "Soft Tissue"
tissue_code2name_gnn[["17"]]
#[1] "Thymus"

true_labels <- dcast(tcga_test_gnn_calls.df,rownames ~ V1, fun.aggregate = function(x) 1L, fill = 0L)
colnames(true_labels) <-c("Rownames", stringr::str_replace_all(paste(colnames(true_labels) %>% .[-1],"_true",sep="")," ","_"))

gnn_calls <- dcast(tcga_test_gnn_calls.df,rownames ~ V4, fun.aggregate = function(x) 1L, fill = 0L)
colnames(gnn_calls) <- c("Rownames", stringr::str_replace_all(paste(tissue_code2name_gnn[gnn_calls %>% colnames() %>% .[-1]] %>% unlist(),"_pred_GNN",sep="")," ","_"))

res_calls <- dcast(tcga_test_resnet_calls.df,rownames ~ V4, fun.aggregate = function(x) 1L, fill = 0L)
colnames(res_calls) <- c("Rownames", stringr::str_replace_all(paste(tissue_code2name_gnn[res_calls %>% colnames() %>% .[-1]] %>% unlist(),"_pred_Resnet",sep="")," ","_"))

gnn_true_v_calls <- data.frame(sapply(cbind(true_labels,gnn_calls[,-1]),as.numeric))
res_true_v_calls <- data.frame(sapply(cbind(true_labels,res_calls[,-1]),as.numeric))

gnn_true_v_calls$Cervix_pred_GNN <- rep(0,nrow(gnn_true_v_calls))
gnn_true_v_calls$Pancreas_pred_GNN <- rep(0,nrow(gnn_true_v_calls))
gnn_true_v_calls$Thymus_pred_GNN <- rep(0,nrow(gnn_true_v_calls))

res_true_v_calls$Pancreas_pred_Resnet <- rep(0,nrow(res_true_v_calls))
res_true_v_calls$Soft_Tissue_pred_Resnet <- rep(0,nrow(res_true_v_calls))
res_true_v_calls$Thymus_pred_Resnet <- rep(0,nrow(res_true_v_calls))

all_true_v_calls <- data.frame(sapply(cbind(true_labels,gnn_calls[,-1]),as.numeric))
all_true_v_calls <- data.frame(sapply(cbind(all_true_v_calls,res_calls[,-1]),as.numeric))

all_true_v_calls$Cervix_pred_GNN <- rep(0,nrow(all_true_v_calls))
all_true_v_calls$Pancreas_pred_GNN <- rep(0,nrow(all_true_v_calls))
all_true_v_calls$Thymus_pred_GNN <- rep(0,nrow(all_true_v_calls))
all_true_v_calls$Pancreas_pred_Resnet <- rep(0,nrow(all_true_v_calls))
all_true_v_calls$Soft_Tissue_pred_Resnet <- rep(0,nrow(all_true_v_calls))
all_true_v_calls$Thymus_pred_Resnet <- rep(0,nrow(all_true_v_calls))

roc_result_gnn <- multiROC::multi_roc(gnn_true_v_calls[,-1],force_diag = TRUE)
roc_result_res <- multiROC::multi_roc(res_true_v_calls[,-1],force_diag = TRUE)

roc_result_all <- multiROC::multi_roc(all_true_v_calls[,-1],force_diag = TRUE)

plot_roc_gnn_df <- plot_roc_data(roc_result_gnn) %>% dplyr::filter(Group != "Micro" & Group != "Macro")
plot_roc_res_df <- plot_roc_data(roc_result_res) %>% dplyr::filter(Group != "Micro" & Group != "Macro")
plot_roc_all_df <- plot_roc_data(roc_result_all) %>% dplyr::filter(Group != "Micro" & Group != "Macro")

# from https://github.com/WandeRum/multiROC

ggplot(plot_roc_all_df, aes(x = 1-Specificity, y=Sensitivity)) +
  geom_path(aes(color = Group, linetype=Method), size=1.5) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
               colour='grey', linetype = 'dotdash') +
  theme_bw() + 
  #theme(plot.title = element_text(hjust = 0.5), 
  #      legend.justification=c(1, 0), legend.position=c(0.95, .05),
  #      legend.title=element_blank(), 
  #      legend.background = element_rect(fill=NULL, size=0.5, 
  #                                       linetype="solid", colour ="black")) +
  ggtitle("GNN vs Resnet MultiROC")

unlist(roc_result_gnn$AUC)
unlist(roc_result_res$AUC)

mclust::adjustedRandIndex(tcga_test_gnn_calls.df$V3,tcga_test_gnn_calls.df$V4)
#[1] 0.7909318
mclust::adjustedRandIndex(tcga_test_resnet_calls.df$V3,tcga_test_resnet_calls.df$V4)
#[1] 0.7781712