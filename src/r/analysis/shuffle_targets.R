library(magrittr)

IN_DIR <- "/home/burkhart/Software/reticula/data/aim2/input/"

# Reaction Network

graph_targets.df <- read.table(paste(IN_DIR,"graph_targets.txt",sep=""),sep = " ")

graph_targets.df$V1 <- sample(graph_targets.df$V1)

write.table(graph_targets.df,
            file=paste(IN_DIR,"shuffled_graph_targets.txt",sep=""),
            row.names = FALSE,
            col.names = FALSE)

# Pathway Hierarchy

graph_targets.df <- read.table(paste(IN_DIR,"pathway_hierarchy_graph_targets.txt",sep=""),sep = " ")

graph_targets.df$V1 <- sample(graph_targets.df$V1)

write.table(graph_targets.df,
            file=paste(IN_DIR,"shuffled_pathway_hierarchy_graph_targets.txt",sep=""),
            row.names = FALSE,
            col.names = FALSE)
