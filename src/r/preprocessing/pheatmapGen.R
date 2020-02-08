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

# combined.df columns ordered by clustering
pheatmap::pheatmap(mat = combined.df,
                   show_rownames = FALSE,
                   show_colnames = FALSE,
                   scale = "row",
                   color = inferno(100),
                   annotation_col = annotation_col.df,
                   annotation_colors = anno_colors,
                   cluster_cols = TRUE, # clustering takes a long time
                   cluster_rows = TRUE, # clustering takes a long time
                   silent = TRUE,
                   filename = "~/combined_pheatmap_rowscale_row_col_clustering.c.png",
                   width=8,
                   height=8,
                   fontsize = 8
)