library(ggplot2)
library(stringr)
DATA_DIR <- "/home/jgburk/PycharmProjects/reticula/data/gtex/output/accuracy_logs/"

fixDecimalDelim <- function(x){
  replaced <- unlist(stringr::str_replace_all(x,"0\\.",",0."))
  replaced <- unlist(stringr::str_replace_all(replaced,"1\\.",",1."))
  if(sum(is.na(replaced))>0){
    print(replaced)
  }
  tokenized <- unlist(stringr::str_split(replaced,","))
  if(sum(is.na(tokenized))>0){
    print(tokenized)
  }
  new_numeric <- numeric()
  for(t_idx in 2:length(tokenized)){
    new_numeric <- c(new_numeric,as.numeric(tokenized[t_idx]))
    if(sum(is.na(new_numeric))>0){
      print(t_idx)
      print(tokenized[t_idx])
      print(new_numeric)
      stop()
    }
  }
  return(new_numeric)
}

# 10-fold CV GTEx GNN
fold1_gnn_acc <- read.table(paste(DATA_DIR,"graph_classification_acc_fold_1.txt",sep=""),sep = ",")$V2
fold2_gnn_acc <- read.table(paste(DATA_DIR,"graph_classification_acc_fold_2.txt",sep=""),sep = ",")$V2
fold3_gnn_acc <- read.table(paste(DATA_DIR,"graph_classification_acc_fold_3.txt",sep=""),sep = ",")$V2
fold4_gnn_acc <- read.table(paste(DATA_DIR,"graph_classification_acc_fold_4.txt",sep=""),sep = ",")$V2
fold5_gnn_acc <- read.table(paste(DATA_DIR,"graph_classification_acc_fold_5.txt",sep=""),sep = ",")$V2
fold6_gnn_acc <- read.table(paste(DATA_DIR,"graph_classification_acc_fold_6.txt",sep=""),sep = ",")$V2
fold7_gnn_acc <- read.table(paste(DATA_DIR,"graph_classification_acc_fold_7.txt",sep=""),sep = ",")$V2
fold8_gnn_acc <- read.table(paste(DATA_DIR,"graph_classification_acc_fold_8.txt",sep=""),sep = ",")$V2
fold9_gnn_acc <- read.table(paste(DATA_DIR,"graph_classification_acc_fold_9.txt",sep=""),sep = ",")$V2
fold10_gnn_acc <- read.table(paste(DATA_DIR,"graph_classification_acc_fold_10.txt",sep=""),sep = ",")$V2
# fullly-trained GNN
gnn_acc <- read.table(paste(DATA_DIR,"graph_classification_acc_real_edges_full_dataset.txt",sep=""),sep=",")$V1
# 10-rewired fully-trained GNN
rewire1_gnn_acc <- read.table(paste(DATA_DIR,"graph_classification_acc_rewire_full_dataset.txt",sep=""))$V1
rewire2_gnn_acc <- fixDecimalDelim(read.table(paste(DATA_DIR,"graph_classification_acc_rewired2_full_dataset.txt",sep=""))$V1)
rewire3_gnn_acc <- fixDecimalDelim(read.table(paste(DATA_DIR,"graph_classification_acc_rewired3_full_dataset.txt",sep=""))$V1)
rewire4_gnn_acc <- fixDecimalDelim(read.table(paste(DATA_DIR,"graph_classification_acc_rewired4_full_dataset.txt",sep=""))$V1)
rewire5_gnn_acc <- fixDecimalDelim(read.table(paste(DATA_DIR,"graph_classification_acc_rewired5_full_dataset.txt",sep=""))$V1)
rewire6_gnn_acc <- fixDecimalDelim(read.table(paste(DATA_DIR,"graph_classification_acc_rewired6_full_dataset.txt",sep=""))$V1)
rewire7_gnn_acc <- fixDecimalDelim(read.table(paste(DATA_DIR,"graph_classification_acc_rewired7_full_dataset.txt",sep=""))$V1)
rewire8_gnn_acc <- fixDecimalDelim(read.table(paste(DATA_DIR,"graph_classification_acc_rewired8_full_dataset.txt",sep=""))$V1)
rewire9_gnn_acc <- fixDecimalDelim(read.table(paste(DATA_DIR,"graph_classification_acc_rewired9_full_dataset.txt",sep=""))$V1)
rewire10_gnn_acc <- fixDecimalDelim(read.table(paste(DATA_DIR,"graph_classification_acc_rewired10_full_dataset.txt",sep=""))$V1)
# Shuffled features fully-trained GNN
shuffled_features_gnn_acc <- fixDecimalDelim(read.table(paste(DATA_DIR,"graph_classification_acc_shuffled_features_full_dataset.txt",sep=""))$V1)
# Shuffled targets fully-trained GNN
shuffled_targets_gnn_acc <- fixDecimalDelim(read.table(paste(DATA_DIR,"graph_classification_acc_shuffled_targets_full_dataset.txt",sep=""))$V1)
# 10-fold CV Resnet
fold1_resnet_acc <- read.table(paste(DATA_DIR,"resnet_resnet_classification_acc_fold_1.txt",sep=""),sep = ",")$V2
fold2_resnet_acc <- read.table(paste(DATA_DIR,"resnet_resnet_classification_acc_fold_2.txt",sep=""),sep = ",")$V2
fold3_resnet_acc <- read.table(paste(DATA_DIR,"resnet_resnet_classification_acc_fold_3.txt",sep=""),sep = ",")$V2
fold4_resnet_acc <- read.table(paste(DATA_DIR,"resnet_resnet_classification_acc_fold_4.txt",sep=""),sep = ",")$V2
fold5_resnet_acc <- read.table(paste(DATA_DIR,"resnet_resnet_classification_acc_fold_5.txt",sep=""),sep = ",")$V2
fold6_resnet_acc <- read.table(paste(DATA_DIR,"resnet_resnet_classification_acc_fold_6.txt",sep=""),sep = ",")$V2
fold7_resnet_acc <- read.table(paste(DATA_DIR,"resnet_resnet_classification_acc_fold_7.txt",sep=""),sep = ",")$V2
fold8_resnet_acc <- read.table(paste(DATA_DIR,"resnet_resnet_classification_acc_fold_8.txt",sep=""),sep = ",")$V2
fold9_resnet_acc <- read.table(paste(DATA_DIR,"resnet_resnet_classification_acc_fold_9.txt",sep=""),sep = ",")$V2
fold10_resnet_acc <- read.table(paste(DATA_DIR,"resnet_resnet_classification_acc_fold_10.txt",sep=""),sep = ",")$V2
# fully-trained Resnet
res_acc <- fixDecimalDelim(read.table(paste(DATA_DIR,"resnet_resnet_classification_acc_fold_full_dataset.txt",sep=""))$V1)
# summation features fully-trained GNN
summation_features_gnn_acc <- read.table(paste(DATA_DIR,"summation_graph_classification_acc_full_dataset.txt",sep=""))$V1

acc_df <- data.frame(f1_gnn = fold1_gnn_acc,
                     f2_gnn = fold2_gnn_acc,
                     f3_gnn = fold3_gnn_acc,
                     f4_gnn = fold4_gnn_acc,
                     f5_gnn = fold5_gnn_acc,
                     f6_gnn = fold6_gnn_acc,
                     f7_gnn = fold7_gnn_acc,
                     f8_gnn = fold8_gnn_acc,
                     f9_gnn = fold9_gnn_acc,
                     f10_gnn = fold10_gnn_acc,
                     f_gnn = gnn_acc,
                     r1_gnn = rewire1_gnn_acc,
                     r2_gnn = rewire2_gnn_acc,
                     r3_gnn = rewire3_gnn_acc,
                     r4_gnn = rewire4_gnn_acc,
                     r5_gnn = rewire5_gnn_acc,
                     r6_gnn = rewire6_gnn_acc,
                     r7_gnn = rewire7_gnn_acc,
                     r8_gnn = rewire8_gnn_acc,
                     r9_gnn = rewire9_gnn_acc,
                     r10_gnn = rewire10_gnn_acc,
                     sf_gnn = shuffled_features_gnn_acc,
                     st_gnn = shuffled_targets_gnn_acc,
                     f1_res = fold1_resnet_acc,
                     f2_res = fold2_resnet_acc,
                     f3_res = fold3_resnet_acc,
                     f4_res = fold4_resnet_acc,
                     f5_res = fold5_resnet_acc,
                     f6_res = fold6_resnet_acc,
                     f7_res = fold7_resnet_acc,
                     f8_res = fold8_resnet_acc,
                     f9_res = fold9_resnet_acc,
                     f10_res = fold10_resnet_acc,
                     f_res = res_acc,
                     sum_res = summation_features_gnn_acc)

#calculate gnn fold mean
acc_df <- acc_df %>% dplyr::mutate(mean_f_gnn = rowMeans(
  dplyr::select(acc_df,c(f1_gnn,
                         f2_gnn,
                         f3_gnn,
                         f4_gnn,
                         f5_gnn,
                         f6_gnn,
                         f7_gnn,
                         f8_gnn,
                         f9_gnn,
                         f10_gnn))
  ))

#calculate rewire mean
acc_df <- acc_df %>% dplyr::mutate(mean_r_gnn = rowMeans(
  dplyr::select(acc_df,c(r1_gnn,
                         r2_gnn,
                         r3_gnn,
                         r4_gnn,
                         r5_gnn,
                         r6_gnn,
                         r7_gnn,
                         r8_gnn,
                         r9_gnn,
                         r10_gnn))
))

#calculate resnet mean  
acc_df <- acc_df %>% dplyr::mutate(mean_f_res = rowMeans(
  dplyr::select(acc_df,c(f1_res,
                         f2_res,
                         f3_res,
                         f4_res,
                         f5_res,
                         f6_res,
                         f7_res,
                         f8_res,
                         f9_res,
                         f10_res))
))

acc_df$rownames <- as.numeric(rownames(acc_df))

final_scores <- subset(acc_df,rownames == max(acc_df$rownames))

control_dark_color <- "#6E1717"
control_light_color <- "#EA9999"
resnet_dark_color <- "#27659B"
resnet_light_color <- "#CFE2F3"
gnn_dark_color <- "#523E7E"
gnn_light_color <- "#D9D2E9"

ggplot(acc_df,aes(x=rownames)) +
  geom_line(aes(y=r1_gnn),color=control_light_color,alpha=0.5) +
  geom_line(aes(y=r2_gnn),color=control_light_color,alpha=0.5) +
  geom_line(aes(y=r3_gnn),color=control_light_color,alpha=0.5) +
  geom_line(aes(y=r4_gnn),color=control_light_color,alpha=0.5) +
  geom_line(aes(y=r5_gnn),color=control_light_color,alpha=0.5) +
  geom_line(aes(y=r6_gnn),color=control_light_color,alpha=0.5) +
  geom_line(aes(y=r7_gnn),color=control_light_color,alpha=0.5) +
  geom_line(aes(y=r8_gnn),color=control_light_color,alpha=0.5) +
  geom_line(aes(y=r9_gnn),color=control_light_color,alpha=0.5) +
  geom_line(aes(y=r10_gnn),color=control_light_color,alpha=0.5) +
  geom_line(aes(y=f1_gnn),color=gnn_light_color,alpha=0.5) +
  geom_line(aes(y=f2_gnn),color=gnn_light_color,alpha=0.5) +
  geom_line(aes(y=f3_gnn),color=gnn_light_color,alpha=0.5) +
  geom_line(aes(y=f4_gnn),color=gnn_light_color,alpha=0.5) +
  geom_line(aes(y=f5_gnn),color=gnn_light_color,alpha=0.5) +
  geom_line(aes(y=f6_gnn),color=gnn_light_color,alpha=0.5) +
  geom_line(aes(y=f7_gnn),color=gnn_light_color,alpha=0.5) +
  geom_line(aes(y=f8_gnn),color=gnn_light_color,alpha=0.5) +
  geom_line(aes(y=f9_gnn),color=gnn_light_color,alpha=0.5) +
  geom_line(aes(y=f10_gnn),color=gnn_light_color,alpha=0.5) +
  geom_line(aes(y=f1_res),color=resnet_light_color,alpha=0.5) +
  geom_line(aes(y=f2_res),color=resnet_light_color,alpha=0.5) +
  geom_line(aes(y=f3_res),color=resnet_light_color,alpha=0.5) +
  geom_line(aes(y=f4_res),color=resnet_light_color,alpha=0.5) +
  geom_line(aes(y=f5_res),color=resnet_light_color,alpha=0.5) +
  geom_line(aes(y=f6_res),color=resnet_light_color,alpha=0.5) +
  geom_line(aes(y=f7_res),color=resnet_light_color,alpha=0.5) +
  geom_line(aes(y=f8_res),color=resnet_light_color,alpha=0.5) +
  geom_line(aes(y=f9_res),color=resnet_light_color,alpha=0.5) +
  geom_line(aes(y=f10_res),color=resnet_light_color,alpha=0.5) +
  geom_smooth(aes(y=mean_f_res),color=resnet_dark_color) +
  geom_smooth(aes(y=sf_gnn),color=control_dark_color) +
  geom_smooth(aes(y=st_gnn),color=control_dark_color) +
  geom_smooth(aes(y=mean_f_gnn),color=gnn_dark_color) +
  geom_smooth(aes(y=mean_r_gnn),color=control_dark_color) +
  geom_smooth(aes(y=sum_res),color=control_dark_color) +
  #geom_smooth(aes(y=f_gnn),color=gnn_dark_color) +
  #geom_smooth(aes(y=f_res),color=resnet_dark_color) +
  # geom_text(data=final_scores, aes(x=rownames,y=mean_f_res),
  #           label=paste("Resnet18 CV (",
  #                       round(final_scores$mean_f_res,4),")",sep=""),
  #           vjust=-.7, hjust=1) +
  # geom_text(data=final_scores,aes(x=rownames,y=sf_gnn),
  #           label=paste("Shuffled Features (",
  #                       round(final_scores$sf_gnn,4),")",sep=""),
  #           vjust=2.7,hjust=1) +
  # geom_text(data=final_scores,aes(x=rownames,y=st_gnn),
  #           label=paste("Shuffled Targets (",
  #                       round(final_scores$st_gnn,4),")",sep=""),
  #           vjust=-.7,hjust=1) +
  # geom_text(data=final_scores, aes(x=rownames,y=mean_f_gnn),
  #           label=paste("Reaction Network CV (",
  #                       round(final_scores$mean_f_gnn,4),")",sep=""),
  #           vjust=-.7, hjust=1) +
  # geom_text(data=final_scores, aes(x=rownames,y=mean_r_gnn),
  #           label=paste("Rewired Network (",
  #                       round(final_scores$mean_r_gnn,4),")",sep=""),
  #           vjust=2.5, hjust=1) +
  # geom_text(data=final_scores, aes(x=rownames,y=sum_res),
  #           label=paste("Summation (",
  #                       round(final_scores$sum_res,4),")",sep=""),
  #           vjust=3.2,hjust=1) +
  #geom_text(data=final_scores, aes(x=rownames,y=f_gnn),
  #          label=paste("Reaction Network (",
  #                      round(final_scores$f_gnn,4),")",sep=""),
  #          vjust=-.7, hjust=1) +
  #geom_text(data=final_scores, aes(x=rownames,y=f_res),
  #          label=paste("Resnet18 (",
  #                      round(final_scores$f_res,4),")",sep=""),
  #          vjust=-.7, hjust=1) +
  theme_minimal()

