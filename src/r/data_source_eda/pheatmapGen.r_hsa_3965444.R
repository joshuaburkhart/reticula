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

reorder.df <- data.frame(Sample=colnames(z),
                         Datasource=datasource.vec,
                         Tissue=tissue.vec)
ordered.tissues.df <- reorder.df %>%
  dplyr::arrange(Tissue)

z1 <- z[,as.character(ordered.tissues.df$Sample)]


pheatmap::pheatmap(mat = z1 - mean(z1),
                   show_rownames = TRUE,
                   show_colnames = FALSE,
                   scale = "row",
                   color = inferno(100),
                   annotation_col = annotation_col.df,
                   annotation_colors = anno_colors,
                   cluster_cols = FALSE, # clustering takes a long time
                   cluster_rows = TRUE, # clustering takes a long time
                   silent = TRUE,
                   filename = "~/3965444_pheatmap_tissue.png",
                   width=8,
                   height=8,
                   fontsize = 8
)
