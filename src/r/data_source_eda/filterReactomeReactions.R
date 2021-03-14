library(magrittr)
library(viridis)
library(RColorBrewer)
library(pheatmap)

vst.rbe <- readRDS("~/vst_rbe.Rds")
sample.transcripts <- rownames(vst.rbe)
sample.transcripts.no.version <- gsub("\\.[0-9]+","",
                                      sample.transcripts)
ensembl2rxns.df <- read.table("~/Downloads/Ensembl2ReactomeReactions.txt",
                           sep="\t",
                           comment.char = "")

sample.transcripts2rxns.df <- ensembl2rxns.df %>%
  dplyr::filter(V1 %in% sample.transcripts.no.version)

rxn.id.instances <- table(as.character(sample.transcripts2rxns.df$V2))

# rxn           freq  name
# R-HSA-381750	75    Olfactory Receptor - G Protein olfactory trimer complex formation
# R-HSA-975040	55    KRAB-ZNF / KAP Interaction
# R-HSA-2730833	50    Phosphorylation of TEC kinases by p-SYK
# R-HSA-2730888	50    Phosphorylation of PLC-gamma

r.hsa.381750 <- sample.transcripts2rxns.df %>%
  dplyr::filter(as.character(V2) == "R-HSA-381750")
r.hsa.975040 <- sample.transcripts2rxns.df %>%
  dplyr::filter(as.character(V2) == "R-HSA-975040")
r.hsa.2730833 <- sample.transcripts2rxns.df %>%
  dplyr::filter(as.character(V2) == "R-HSA-2730833")
r.hsa.2730888 <- sample.transcripts2rxns.df %>%
  dplyr::filter(as.character(V2) == "R-HSA-2730888")

rownames(vst.rbe) <- sample.transcripts.no.version

vst.rbe.381750 <- vst.rbe[as.character(r.hsa.381750$V1),]
vst.rbe.975040 <- vst.rbe[as.character(r.hsa.975040$V1),]
vst.rbe.2730833 <- vst.rbe[as.character(r.hsa.2730833$V1),]
vst.rbe.2730888 <- vst.rbe[as.character(r.hsa.2730888$V1),]

tissue.vec <- readRDS("~/tissue_vec.Rds")
datasource.vec <- readRDS("~/datasource_vec.Rds")

colFun <- colorRampPalette(RColorBrewer::brewer.pal(10,"Paired"))

annotation_col.df <- data.frame(Datasource = datasource.vec,
                                Tissue = tissue.vec)
rownames(annotation_col.df) <- colnames(vst.rbe)

Datasource <- colFun(length(unique(datasource.vec)))
Tissue <- colFun(length(unique(tissue.vec)))

names(Datasource) <- annotation_col.df$Datasource %>% unique()
names(Tissue) <- annotation_col.df$Tissue %>% unique()

anno_colors <- list(Datasource = Datasource,
                    Tissue = Tissue)

# columns ordered by clustering
pheatmap::pheatmap(mat = vst.rbe.381750,
                   show_rownames = TRUE,
                   show_colnames = FALSE,
                   scale = "row",
                   color = inferno(100),
                   annotation_col = annotation_col.df,
                   annotation_colors = anno_colors,
                   cluster_cols = TRUE, # clustering takes a long time
                   cluster_rows = TRUE, # clustering takes a long time
                   silent = TRUE,
                   filename = "~/vst.rbe.381750.png",
                   width=8,
                   height=8,
                   fontsize = 8
)

pheatmap::pheatmap(mat = vst.rbe.975040,
                   show_rownames = TRUE,
                   show_colnames = FALSE,
                   scale = "row",
                   color = inferno(100),
                   annotation_col = annotation_col.df,
                   annotation_colors = anno_colors,
                   cluster_cols = TRUE, # clustering takes a long time
                   cluster_rows = TRUE, # clustering takes a long time
                   silent = TRUE,
                   filename = "~/vst.rbe.975040.png",
                   width=8,
                   height=8,
                   fontsize = 8
)

pheatmap::pheatmap(mat = vst.rbe.2730833,
                   show_rownames = TRUE,
                   show_colnames = FALSE,
                   scale = "row",
                   color = inferno(100),
                   annotation_col = annotation_col.df,
                   annotation_colors = anno_colors,
                   cluster_cols = TRUE, # clustering takes a long time
                   cluster_rows = TRUE, # clustering takes a long time
                   silent = TRUE,
                   filename = "~/vst.rbe.2730833.png",
                   width=8,
                   height=8,
                   fontsize = 8
)

pheatmap::pheatmap(mat = vst.rbe.2730888,
                   show_rownames = TRUE,
                   show_colnames = FALSE,
                   scale = "row",
                   color = inferno(100),
                   annotation_col = annotation_col.df,
                   annotation_colors = anno_colors,
                   cluster_cols = TRUE, # clustering takes a long time
                   cluster_rows = TRUE, # clustering takes a long time
                   silent = TRUE,
                   filename = "~/vst.rbe.2730888.png",
                   width=8,
                   height=8,
                   fontsize = 8
)

archs4.filter <- datasource.vec == "ARCHS4"

pheatmap::pheatmap(mat = vst.rbe.381750[,archs4.filter],
                   show_rownames = TRUE,
                   show_colnames = FALSE,
                   scale = "row",
                   color = inferno(100),
                   annotation_col = annotation_col.df,
                   annotation_colors = anno_colors,
                   cluster_cols = TRUE, # clustering takes a long time
                   cluster_rows = TRUE, # clustering takes a long time
                   silent = TRUE,
                   filename = "~/vst.rbe.381750.archs4.png",
                   width=8,
                   height=8,
                   fontsize = 8
)

pheatmap::pheatmap(mat = vst.rbe.975040[,archs4.filter],
                   show_rownames = TRUE,
                   show_colnames = FALSE,
                   scale = "row",
                   color = inferno(100),
                   annotation_col = annotation_col.df,
                   annotation_colors = anno_colors,
                   cluster_cols = TRUE, # clustering takes a long time
                   cluster_rows = TRUE, # clustering takes a long time
                   silent = TRUE,
                   filename = "~/vst.rbe.975040.archs4.png",
                   width=8,
                   height=8,
                   fontsize = 8
)

pheatmap::pheatmap(mat = vst.rbe.2730833[,archs4.filter],
                   show_rownames = TRUE,
                   show_colnames = FALSE,
                   scale = "row",
                   color = inferno(100),
                   annotation_col = annotation_col.df,
                   annotation_colors = anno_colors,
                   cluster_cols = TRUE, # clustering takes a long time
                   cluster_rows = TRUE, # clustering takes a long time
                   silent = TRUE,
                   filename = "~/vst.rbe.2730833.archs4.png",
                   width=8,
                   height=8,
                   fontsize = 8
)

pheatmap::pheatmap(mat = vst.rbe.2730888[,archs4.filter],
                   show_rownames = TRUE,
                   show_colnames = FALSE,
                   scale = "row",
                   color = inferno(100),
                   annotation_col = annotation_col.df,
                   annotation_colors = anno_colors,
                   cluster_cols = TRUE, # clustering takes a long time
                   cluster_rows = TRUE, # clustering takes a long time
                   silent = TRUE,
                   filename = "~/vst.rbe.2730888.archs4.png",
                   width=8,
                   height=8,
                   fontsize = 8
)

gtex.filter <- datasource.vec == "GTEx"

pheatmap::pheatmap(mat = vst.rbe.381750[,gtex.filter],
                   show_rownames = TRUE,
                   show_colnames = FALSE,
                   scale = "row",
                   color = inferno(100),
                   annotation_col = annotation_col.df,
                   annotation_colors = anno_colors,
                   cluster_cols = TRUE, # clustering takes a long time
                   cluster_rows = TRUE, # clustering takes a long time
                   silent = TRUE,
                   filename = "~/vst.rbe.381750.gtex.png",
                   width=8,
                   height=8,
                   fontsize = 8
)

pheatmap::pheatmap(mat = vst.rbe.975040[,gtex.filter],
                   show_rownames = TRUE,
                   show_colnames = FALSE,
                   scale = "row",
                   color = inferno(100),
                   annotation_col = annotation_col.df,
                   annotation_colors = anno_colors,
                   cluster_cols = TRUE, # clustering takes a long time
                   cluster_rows = TRUE, # clustering takes a long time
                   silent = TRUE,
                   filename = "~/vst.rbe.975040.gtex.png",
                   width=8,
                   height=8,
                   fontsize = 8
)

pheatmap::pheatmap(mat = vst.rbe.2730833[,gtex.filter],
                   show_rownames = TRUE,
                   show_colnames = FALSE,
                   scale = "row",
                   color = inferno(100),
                   annotation_col = annotation_col.df,
                   annotation_colors = anno_colors,
                   cluster_cols = TRUE, # clustering takes a long time
                   cluster_rows = TRUE, # clustering takes a long time
                   silent = TRUE,
                   filename = "~/vst.rbe.2730833.gtex.png",
                   width=8,
                   height=8,
                   fontsize = 8
)

pheatmap::pheatmap(mat = vst.rbe.2730888[,gtex.filter],
                   show_rownames = TRUE,
                   show_colnames = FALSE,
                   scale = "row",
                   color = inferno(100),
                   annotation_col = annotation_col.df,
                   annotation_colors = anno_colors,
                   cluster_cols = TRUE, # clustering takes a long time
                   cluster_rows = TRUE, # clustering takes a long time
                   silent = TRUE,
                   filename = "~/vst.rbe.2730888.gtex.png",
                   width=8,
                   height=8,
                   fontsize = 8
)

tcga.filter <- datasource.vec == "TCGA"

pheatmap::pheatmap(mat = vst.rbe.381750[,tcga.filter],
                   show_rownames = TRUE,
                   show_colnames = FALSE,
                   scale = "row",
                   color = inferno(100),
                   annotation_col = annotation_col.df,
                   annotation_colors = anno_colors,
                   cluster_cols = TRUE, # clustering takes a long time
                   cluster_rows = TRUE, # clustering takes a long time
                   silent = TRUE,
                   filename = "~/vst.rbe.381750.tcga.png",
                   width=8,
                   height=8,
                   fontsize = 8
)

pheatmap::pheatmap(mat = vst.rbe.975040[,tcga.filter],
                   show_rownames = TRUE,
                   show_colnames = FALSE,
                   scale = "row",
                   color = inferno(100),
                   annotation_col = annotation_col.df,
                   annotation_colors = anno_colors,
                   cluster_cols = TRUE, # clustering takes a long time
                   cluster_rows = TRUE, # clustering takes a long time
                   silent = TRUE,
                   filename = "~/vst.rbe.975040.tcga.png",
                   width=8,
                   height=8,
                   fontsize = 8
)

pheatmap::pheatmap(mat = vst.rbe.2730833[,tcga.filter],
                   show_rownames = TRUE,
                   show_colnames = FALSE,
                   scale = "row",
                   color = inferno(100),
                   annotation_col = annotation_col.df,
                   annotation_colors = anno_colors,
                   cluster_cols = TRUE, # clustering takes a long time
                   cluster_rows = TRUE, # clustering takes a long time
                   silent = TRUE,
                   filename = "~/vst.rbe.2730833.tcga.png",
                   width=8,
                   height=8,
                   fontsize = 8
)

pheatmap::pheatmap(mat = vst.rbe.2730888[,tcga.filter],
                   show_rownames = TRUE,
                   show_colnames = FALSE,
                   scale = "row",
                   color = inferno(100),
                   annotation_col = annotation_col.df,
                   annotation_colors = anno_colors,
                   cluster_cols = TRUE, # clustering takes a long time
                   cluster_rows = TRUE, # clustering takes a long time
                   silent = TRUE,
                   filename = "~/vst.rbe.2730888.tcga.png",
                   width=8,
                   height=8,
                   fontsize = 8
)

