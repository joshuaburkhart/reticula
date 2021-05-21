#resnet 10fold CV vs reaction network 10 foldCV vs pathway hierarchy 10fold CV with lines for random tissue label and reaction pc1 shuffles
library(magrittr)
library(dplyr)
library(ggplot2)

OUT_DIR <- "/home/burkhart/Software/reticula/data/aim2/output/"

# resnet 10fold CV
resnet_10foldCV_acc_fns <- c("resnet_resnet_classification_acc_fold_1.txt",
                             "resnet_resnet_classification_acc_fold_2.txt",
                             "resnet_resnet_classification_acc_fold_3.txt",
                             "resnet_resnet_classification_acc_fold_4.txt",
                             "resnet_resnet_classification_acc_fold_5.txt",
                             "resnet_resnet_classification_acc_fold_6.txt",
                             "resnet_resnet_classification_acc_fold_7.txt",
                             "resnet_resnet_classification_acc_fold_8.txt",
                             "resnet_resnet_classification_acc_fold_9.txt",
                             "resnet_resnet_classification_acc_fold_10.txt")

# reaction network 10foldCV
rxnnet_10foldCV_acc_fns <- c("graph_classification_acc_fold_1.txt",
                             "graph_classification_acc_fold_2.txt",
                             "graph_classification_acc_fold_3.txt",
                             "graph_classification_acc_fold_4.txt",
                             "graph_classification_acc_fold_5.txt",
                             "graph_classification_acc_fold_6.txt",
                             "graph_classification_acc_fold_7.txt",
                             "graph_classification_acc_fold_8.txt",
                             "graph_classification_acc_fold_9.txt",
                             "graph_classification_acc_fold_10.txt")

# pathway hierarchy 10foldCV

ph_10foldCV_acc_fns <- c("pathway_hierarchy_graph_classification_acc_fold_1.txt",
                         "pathway_hierarchy_graph_classification_acc_fold_2.txt",
                         "pathway_hierarchy_graph_classification_acc_fold_3.txt",
                         "pathway_hierarchy_graph_classification_acc_fold_4.txt",
                         "pathway_hierarchy_graph_classification_acc_fold_5.txt",
                         "pathway_hierarchy_graph_classification_acc_fold_6.txt",
                         "pathway_hierarchy_graph_classification_acc_fold_7.txt",
                         "pathway_hierarchy_graph_classification_acc_fold_8.txt",
                         "pathway_hierarchy_graph_classification_acc_fold_9.txt",
                         "pathway_hierarchy_graph_classification_acc_fold_10.txt")

performance_controls.df <- data.frame(
  res1 = read.csv(file=paste(OUT_DIR,resnet_10foldCV_acc_fns[1],sep=""),header = FALSE) %>% .[2],
  res2 = read.csv(file=paste(OUT_DIR,resnet_10foldCV_acc_fns[2],sep=""),header = FALSE) %>% .[2],
  res3 = read.csv(file=paste(OUT_DIR,resnet_10foldCV_acc_fns[3],sep=""),header = FALSE) %>% .[2],
  res4 = read.csv(file=paste(OUT_DIR,resnet_10foldCV_acc_fns[4],sep=""),header = FALSE) %>% .[2],
  res5 = read.csv(file=paste(OUT_DIR,resnet_10foldCV_acc_fns[5],sep=""),header = FALSE) %>% .[2],
  res6 = read.csv(file=paste(OUT_DIR,resnet_10foldCV_acc_fns[6],sep=""),header = FALSE) %>% .[2],
  res7 = read.csv(file=paste(OUT_DIR,resnet_10foldCV_acc_fns[7],sep=""),header = FALSE) %>% .[2],
  res8 = read.csv(file=paste(OUT_DIR,resnet_10foldCV_acc_fns[8],sep=""),header = FALSE) %>% .[2],
  res9 = read.csv(file=paste(OUT_DIR,resnet_10foldCV_acc_fns[9],sep=""),header = FALSE) %>% .[2],
  res10 = read.csv(file=paste(OUT_DIR,resnet_10foldCV_acc_fns[10],sep=""),header = FALSE) %>% .[2],
  rxn1 = read.csv(file=paste(OUT_DIR,rxnnet_10foldCV_acc_fns[1],sep=""),header = FALSE) %>% .[2],
  rxn2 = read.csv(file=paste(OUT_DIR,rxnnet_10foldCV_acc_fns[2],sep=""),header = FALSE) %>% .[2],
  rxn3 = read.csv(file=paste(OUT_DIR,rxnnet_10foldCV_acc_fns[3],sep=""),header = FALSE) %>% .[2],
  rxn4 = read.csv(file=paste(OUT_DIR,rxnnet_10foldCV_acc_fns[4],sep=""),header = FALSE) %>% .[2],
  rxn5 = read.csv(file=paste(OUT_DIR,rxnnet_10foldCV_acc_fns[5],sep=""),header = FALSE) %>% .[2],
  rxn6 = read.csv(file=paste(OUT_DIR,rxnnet_10foldCV_acc_fns[6],sep=""),header = FALSE) %>% .[2],
  rxn7 = read.csv(file=paste(OUT_DIR,rxnnet_10foldCV_acc_fns[7],sep=""),header = FALSE) %>% .[2],
  rxn8 = read.csv(file=paste(OUT_DIR,rxnnet_10foldCV_acc_fns[8],sep=""),header = FALSE) %>% .[2],
  rxn9 = read.csv(file=paste(OUT_DIR,rxnnet_10foldCV_acc_fns[9],sep=""),header = FALSE) %>% .[2],
  rxn10 = read.csv(file=paste(OUT_DIR,rxnnet_10foldCV_acc_fns[10],sep=""),header = FALSE) %>% .[2],
  ph1 = read.csv(file=paste(OUT_DIR,ph_10foldCV_acc_fns[1],sep=""),header = FALSE) %>% .[2] %>% na.omit(),
  ph2 = read.csv(file=paste(OUT_DIR,ph_10foldCV_acc_fns[2],sep=""),header = FALSE) %>% .[2] %>% na.omit(),
  ph3 = read.csv(file=paste(OUT_DIR,ph_10foldCV_acc_fns[3],sep=""),header = FALSE) %>% .[2] %>% na.omit(),
  ph4 = read.csv(file=paste(OUT_DIR,ph_10foldCV_acc_fns[4],sep=""),header = FALSE) %>% .[2] %>% na.omit(),
  ph5 = read.csv(file=paste(OUT_DIR,ph_10foldCV_acc_fns[5],sep=""),header = FALSE) %>% .[2] %>% na.omit(),
  ph6 = read.csv(file=paste(OUT_DIR,ph_10foldCV_acc_fns[6],sep=""),header = FALSE) %>% .[2] %>% na.omit(),
  ph7 = read.csv(file=paste(OUT_DIR,ph_10foldCV_acc_fns[7],sep=""),header = FALSE) %>% .[2] %>% na.omit(),
  ph8 = read.csv(file=paste(OUT_DIR,ph_10foldCV_acc_fns[8],sep=""),header = FALSE) %>% .[2] %>% na.omit(),
  ph9 = read.csv(file=paste(OUT_DIR,ph_10foldCV_acc_fns[9],sep=""),header = FALSE) %>% .[2] %>% na.omit(),
  ph10 = read.csv(file=paste(OUT_DIR,ph_10foldCV_acc_fns[10],sep=""),header = FALSE) %>% .[2] %>% na.omit(),
  shuffled_features_ph = read.csv(file=paste(OUT_DIR,"pathway_hierarchy_graph_classification_acc_shuffle_features_fold_full_dataset.txt",sep=""),header = FALSE),
  shuffled_targets_ph = read.csv(file=paste(OUT_DIR,"pathway_hierarchy_graph_classification_acc_shuffled_targets_fold_full_dataset.txt",sep=""),header = FALSE),
  shuffled_features_rxn = read.csv(file=paste(OUT_DIR,"graph_classification_acc_shuffled_features_full_dataset.txt",sep=""),header = FALSE),
  shuffled_targets_rxn = read.csv(file=paste(OUT_DIR,"graph_classification_acc_shuffled_targets_full_dataset.txt",sep=""),header = FALSE)
)

performance_controls.df <- performance_controls.df %>% dplyr::mutate(res_mean = rowMeans(select(.,c(1,2,3,4,5,6,7,8,9,10))),
                                                                     rxn_mean = rowMeans(select(.,c(11,12,13,14,15,16,17,18,19,20))),
                                                                     ph_mean = rowMeans(select(.,c(21,22,23,24,25,26,27,28,29,30)))) %>%
  dplyr::select(res_mean,rxn_mean,ph_mean,31,32,33,34)

colnames(performance_controls.df) <- c("ResNet18","Reaction Network (RXN)","Pathway Hierarchy (PH)","PH Shuffled Features","PH Shuffled Targets","RXN Shuffled Features","RXN Shuffled Targets")

performance_controls.df$epoch <- seq(1:nrow(performance_controls.df))

plot(performance_controls.df$ResNet18,type = "l",col="blue",xlim=c(1,500),ylim=c(0,1))
lines(performance_controls.df$`Reaction Network (RXN)`,col="purple")
lines(performance_controls.df$`Pathway Hierarchy (PH)`,col="green")
lines(performance_controls.df$`PH Shuffled Features`,col="dark green")
lines(performance_controls.df$`PH Shuffled Targets`,col="light green")
lines(performance_controls.df$`RXN Shuffled Features`,col="darkorchid4")
lines(performance_controls.df$`RXN Shuffled Targets`,col="mediumpurple")


# 10fold CV reaction network vs 10 degree-preserving randomizations

# 10foldCV reaction network

# reaction network 10foldCV
rxnnet_10foldCV_acc_fns <- c("graph_classification_acc_fold_1.txt",
                             "graph_classification_acc_fold_2.txt",
                             "graph_classification_acc_fold_3.txt",
                             "graph_classification_acc_fold_4.txt",
                             "graph_classification_acc_fold_5.txt",
                             "graph_classification_acc_fold_6.txt",
                             "graph_classification_acc_fold_7.txt",
                             "graph_classification_acc_fold_8.txt",
                             "graph_classification_acc_fold_9.txt",
                             "graph_classification_acc_fold_10.txt")

# rewirings

rewired_10foldCV_acc_fns <- c("graph_classification_acc_rewire_full_dataset.txt",
                         "graph_classification_acc_rewired2_full_dataset.txt",
                         "graph_classification_acc_rewired3_full_dataset.txt",
                         "graph_classification_acc_rewired4_full_dataset.txt",
                         "graph_classification_acc_rewired5_full_dataset.txt",
                         "graph_classification_acc_rewired6_full_dataset.txt",
                         "graph_classification_acc_rewired7_full_dataset.txt",
                         "graph_classification_acc_rewired8_full_dataset.txt",
                         "graph_classification_acc_rewired9_full_dataset.txt",
                         "graph_classification_acc_rewired10_full_dataset.txt")

structure_assessment.df <- data.frame(
  rxn1 = read.csv(file=paste(OUT_DIR,rxnnet_10foldCV_acc_fns[1],sep=""),header = FALSE) %>% .[2],
  rxn2 = read.csv(file=paste(OUT_DIR,rxnnet_10foldCV_acc_fns[2],sep=""),header = FALSE) %>% .[2],
  rxn3 = read.csv(file=paste(OUT_DIR,rxnnet_10foldCV_acc_fns[3],sep=""),header = FALSE) %>% .[2],
  rxn4 = read.csv(file=paste(OUT_DIR,rxnnet_10foldCV_acc_fns[4],sep=""),header = FALSE) %>% .[2],
  rxn5 = read.csv(file=paste(OUT_DIR,rxnnet_10foldCV_acc_fns[5],sep=""),header = FALSE) %>% .[2],
  rxn6 = read.csv(file=paste(OUT_DIR,rxnnet_10foldCV_acc_fns[6],sep=""),header = FALSE) %>% .[2],
  rxn7 = read.csv(file=paste(OUT_DIR,rxnnet_10foldCV_acc_fns[7],sep=""),header = FALSE) %>% .[2],
  rxn8 = read.csv(file=paste(OUT_DIR,rxnnet_10foldCV_acc_fns[8],sep=""),header = FALSE) %>% .[2],
  rxn9 = read.csv(file=paste(OUT_DIR,rxnnet_10foldCV_acc_fns[9],sep=""),header = FALSE) %>% .[2],
  rxn10 = read.csv(file=paste(OUT_DIR,rxnnet_10foldCV_acc_fns[10],sep=""),header = FALSE) %>% .[2],
  rew1 = read.csv(file=paste(OUT_DIR,rewired_10foldCV_acc_fns[1],sep=""),header = FALSE),
  rew2 = read.csv(file=paste(OUT_DIR,rewired_10foldCV_acc_fns[2],sep=""),header = FALSE),
  rew3 = read.csv(file=paste(OUT_DIR,rewired_10foldCV_acc_fns[3],sep=""),header = FALSE),
  rew4 = read.csv(file=paste(OUT_DIR,rewired_10foldCV_acc_fns[4],sep=""),header = FALSE),
  rew5 = read.csv(file=paste(OUT_DIR,rewired_10foldCV_acc_fns[5],sep=""),header = FALSE),
  rew6 = read.csv(file=paste(OUT_DIR,rewired_10foldCV_acc_fns[6],sep=""),header = FALSE),
  rew7 = read.csv(file=paste(OUT_DIR,rewired_10foldCV_acc_fns[7],sep=""),header = FALSE),
  rew8 = read.csv(file=paste(OUT_DIR,rewired_10foldCV_acc_fns[8],sep=""),header = FALSE),
  rew9 = read.csv(file=paste(OUT_DIR,rewired_10foldCV_acc_fns[9],sep=""),header = FALSE),
  rew10 = read.csv(file=paste(OUT_DIR,rewired_10foldCV_acc_fns[10],sep=""),header = FALSE)
)
  
real <- structure_assessment.df %>% .[500,1:10] %>% as.numeric()
rewire <- structure_assessment.df %>% .[500,11:20] %>% as.numeric()

plotting_df <- data.frame(accuracy = c(real,rewire),
                          category = c(rep("real",10),rep("rewire",10)))
  
  ggplot(plotting_df, aes(x=accuracy,
                                  y=category,
                                  color=category)) +
    geom_violin() +
    coord_flip() +
    geom_boxplot(width  =0.15) +
    theme_bw() +
    theme(legend.position = "none") +
    scale_color_manual(values = c("purple","red"))
  
  wilcox.test(real,rewire,alternative = "greater")