library(magrittr)

IN_DIR <- "/home/burkhart/Software/reticula/data/aim2/input/"

# Reaction Network

edges.df <- read.table(paste(IN_DIR,"edges.txt",sep=""),sep = " ")

unique_nodes.vec <- c(edges.df$V1,edges.df$V2) %>% unique()

edges.df$V1 <- sample(unique_nodes.vec,nrow(edges.df),replace = TRUE)
edges.df$V2 <- sample(unique_nodes.vec,nrow(edges.df),replace = TRUE)

write.table(edges.df,
            file=paste(IN_DIR,"random_rewired_edges.txt",sep=""),
            row.names = FALSE,
            col.names = FALSE)

# Pathway Hierarchy

edges.df <- read.table(paste(IN_DIR,"pathway_hierarchy_edges.txt",sep=""),sep = " ")

unique_nodes.vec <- c(edges.df$V1,edges.df$V2) %>% unique()

edges.df$V1 <- sample(unique_nodes.vec,nrow(edges.df),replace = TRUE)
edges.df$V2 <- sample(unique_nodes.vec,nrow(edges.df),replace = TRUE)

write.table(edges.df,
            file=paste(IN_DIR,"random_rewired_pathway_hierarchy_edges.txt",sep=""),
            row.names = FALSE,
            col.names = FALSE)
