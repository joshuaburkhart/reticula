library(viridis)
library(pheatmap)
library(RColorBrewer)

combined.df <- readRDS("~/combined_df.Rds")
tissue.vec <- readRDS("~/tissue_vec.Rds")
datasource.vec <- readRDS("~/datasource_vec.Rds")

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

reorder.df <- data.frame(Sample=colnames(combined.df),
                         Datasource=datasource.vec,
                         Tissue=tissue.vec)
ordered.tissues.df <- reorder.df %>%
  dplyr::arrange(Tissue)

combined.df <- combined.df[,as.character(ordered.tissues.df$Sample)]

# combined.df columns ordered by tissue
pheatmap::pheatmap(mat = combined.df,
                   show_rownames = FALSE,
                   show_colnames = FALSE,
                   scale = "row",
                   color = inferno(100),
                   annotation_col = annotation_col.df,
                   annotation_colors = anno_colors,
                   cluster_cols = FALSE, # clustering takes a long time
                   cluster_rows = TRUE, # clustering takes a long time
                   silent = TRUE,
                   filename = "~/combined_pheatmap_rowscale_rowcluster_tissue.c.png",
                   width=8,
                   height=8,
                   fontsize = 8
)