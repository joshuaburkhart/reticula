library(multiROC)
library(magrittr)
library(ggplot2)
library(data.table)
library(stringr)
library(dplyr)

DATA_DIR <- "/home/jgburk/PycharmProjects/reticula/data/tcga/output/"

generate_misclass_chord <- function(n_elements,nm_elements,tcga_test_calls_df,tissue_code2name,plot_name){
  misclass_df <- data.frame(matrix(data=0,ncol=n_elements,nrow = n_elements))
  rownames(misclass_df) <- nm_elements
  colnames(misclass_df) <- nm_elements
  
  for(i in 0:n_elements) {
    tis_name <- tissue_code2name[[as.character(i)]]
    if (length(tis_name) > 0) {
      tis_calls <-
        sapply(tcga_test_calls_df[which(tcga_test_calls_df$V1 == tis_name), ] %>% .$V4, function(x)
          tissue_code2name[[as.character(x)]])
      tis_miscalls <- tis_calls[which(tis_calls != tis_name)]
      n_miscalls <- length(tis_miscalls)
      if (n_miscalls > 0) {
        for (j in 1:n_miscalls) {
          miscall_tis_name <- tis_miscalls[j]
          misclass_df[tis_name, miscall_tis_name] <-
            misclass_df[tis_name, miscall_tis_name] + 1
        }
      } else{
        print(paste("INFO: No GNN miscalls for ", tis_name, sep = ""))
      }
    } else{
      print(paste("ERROR: No tissue name for i=", i, sep = ""))
    }
  }
  
  # from https://www.data-to-viz.com/graph/chord.html#code
  
  library(tidyverse)
  library(viridis)
  library(patchwork)
  library(hrbrthemes)
  library(circlize)
  library(chorddiag)  #devtools::install_github("mattflor/chorddiag")
  
  # I need a long format
  data_long <- misclass_df %>%
    rownames_to_column %>%
    gather(key = 'key', value = 'value', -rowname)
  
  # parameters
  circos.clear()
  circos.par(start.degree = 90, gap.degree = 4, track.margin = c(-0.1, 0.1), points.overflow.warning = FALSE)
  par(mar = rep(0, 4))
  
  # color palette
  mycolor <- viridis(10, alpha = 1, begin = 0, end = 1, option = "D")
  mycolor <- mycolor[sample(1:10)]
  
  svg(paste("~/",plot_name,sep=""), width = 11, height = 10)
  # Base plot
  chordDiagram(
    x = data_long, 
    grid.col = c("#F7756C",  # Adrenal Gland
                 "#E9842C",  # Bile Duct
                 "#D69101",  # Bladder
                 "#BC9D01",  # Brain
                 "#9CA700",  # Breast
                 "#6FB000",  # Cervix
                 "#00B712",  # Colorectal
                 "#00BC60",  # Esophagus
                 "#00C08E",  # Head and Neck
                 "#18C6BB",  # Kidney
                 "#08BFD5",  # Liver
                 "#00B5EE",  # Lung
                 "#00A6FF",  # Pancreas
                 "#7F96FF",  # Prostate
                 "#BC82FF",  # Soft Tissue
                 "#E26EF7",  # Stomach
                 "#F762DE",  # Thymus
                 "#FF62BF",  # Thyroid
                 "#FF6C9B"), # Uterus
    transparency = 0.25,
    directional = 1,
    direction.type = c("arrows", "diffHeight"), 
    diffHeight  = -0.04,
    annotationTrack = "grid", 
    annotationTrackHeight = c(0.05, 0.1),
    link.arr.type = "big.arrow", 
    link.sort = TRUE, 
    link.largest.ontop = TRUE)
  
  # Add text and axis
  skip_section_counter <<- 0
  circos.trackPlotRegion(
    track.index = 1, 
    bg.border = NA, 
    panel.fun = function(x, y) {
      
      xlim = get.cell.meta.data("xlim")
      sector.index = get.cell.meta.data("sector.index")
      
      # Add names to the sector. 
      circos.text(
        x = mean(xlim), 
        y = 4,                     # Section label distance outside edge (3.2 default)
        labels = sector.index, 
        facing = "downward", 
        cex = 1                    # Font size (0.8 default)
      )
      
      section_idx <- get.cell.meta.data("sector.numeric.index")
      
      n_section_miscalls <- rowSums(misclass_df) %>% .[section_idx]
      n_section_incoming <- colSums(misclass_df) %>% .[section_idx]
      
      if(n_section_miscalls == 0 && n_section_incoming == 0){
        skip_section_counter <<- skip_section_counter + 1
        print(paste("skip_section_counter incremented to ",skip_section_counter,sep=""))
      }
      
      row_sum_idx <- section_idx + skip_section_counter
      
      n_section_miscalls <- rowSums(misclass_df) %>% .[row_sum_idx]
      
      n_ticks <- n_section_miscalls
      
      print(paste("Adding ",n_ticks," from row index ",row_sum_idx," for section ",section_idx,
                  " '",get.cell.meta.data("sector.index"),"'",sep=""))
      
      # Add graduation on axis
      circos.axis(
        h = "top", 
        major.at = seq(from = 0, to = n_ticks,
                       by = 1), 
        minor.ticks = 1,
        labels.cex = 0.5,               # Tick Label Font Size (0.5 default)
        #major.tick.percentage = 0.5,
        #major.tick.length = 1,           # Tick Length (unset default)
        #lwd = 5,                         # Tick width (unset default)
        labels.niceFacing = FALSE
      )
    }
  )
  dev.off()
}

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

stopifnot(assertthat::are_equal(tcga_test_gnn_calls.df$V1 %>% unique() %>% length(),
                                tcga_test_resnet_calls.df$V1 %>% unique() %>% length()))

num_elements <- tcga_test_gnn_calls.df$V1 %>% unique() %>% length()
name_elements <- tcga_test_gnn_calls.df$V1 %>% unique()

generate_misclass_chord(n_elements = num_elements,
                        nm_elements = name_elements,
                        tcga_test_calls_df = tcga_test_gnn_calls.df,
                        tissue_code2name = tissue_code2name_gnn,
                        plot_name = "Resnet_Misclass_Chord.svg")

num_elements <- tcga_test_resnet_calls.df$V1 %>% unique() %>% length()
name_elements <- tcga_test_resnet_calls.df$V1 %>% unique()

generate_misclass_chord(n_elements = num_elements,
                        nm_elements = name_elements,
                        tcga_test_calls_df = tcga_test_resnet_calls.df,
                        tissue_code2name = tissue_code2name_res,
                        "GNN_Misclass_Chord.svg")
