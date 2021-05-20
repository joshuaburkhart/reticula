#resnet 10fold CV vs reaction network 10 foldCV vs pathway hierarchy 10fold CV with lines for random tissue label and reaction pc1 shuffles



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