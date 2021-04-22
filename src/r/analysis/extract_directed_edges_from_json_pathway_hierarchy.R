library(tidyjson)

IN_DIR <- "/home/burkhart/Software/reticula/data/aim2/input/"

json <- tidyjson::read_json(path=paste(IN_DIR,"reactome_pathway_hierarchy.json",sep=""))
