
IN_DIR <- "/home/burkhart/Software/reticula/data/aim2/input/"

#https://arxiv.org/pdf/1202.3473.pdf "the results in Sec. 3.1 furnish a workable estimate of N, if one uses <10−5."

# Reaction Network

edges.df <- read.table(paste(IN_DIR,"edges.txt",sep=""),sep = " ")

EPSILON <- 10^-7
E <- nrow(edges.df)
N_REWIRES <- ceiling(E * log(1/EPSILON))

for(REWIRE in 1:N_REWIRES){
  e_idxs <- sample(E,2)
  rands <- rbinom(2,1,0.5)
  p_rxn1 <- edges.df$V1[e_idxs[1]]
  f_rxn1 <- edges.df$V2[e_idxs[1]]
  p_rxn2 <- edges.df$V1[e_idxs[2]]
  f_rxn2 <- edges.df$V2[e_idxs[2]]
  if(rands[1] > 0){
    temp <- p_rxn1
    p_rxn1 <- f_rxn2
    f_rxn2 <- temp
  }
  if(rands[2] > 0){
    temp <- p_rxn2
    p_rxn2 <- f_rxn1
    f_rxn1 <- temp
  }
  edges.df[e_idxs[1],] <- c(p_rxn1,f_rxn2)
  edges.df[e_idxs[2],] <- c(p_rxn2,f_rxn1)
  if(mod(REWIRE,10000) == 0){
    print(paste("Reaction network rewire ",REWIRE," of ",N_REWIRES," complete...",sep=""))
  }
}

write.table(edges.df,
            file=paste(IN_DIR,"rewired_edges.txt",sep=""),
            row.names = FALSE,
            col.names = FALSE)

# Pathway Hierarchy

edges.df <- read.table(paste(IN_DIR,"pathway_hierarchy_edges.txt",sep=""),sep = " ")

EPSILON <- 10^-7
E <- nrow(edges.df)
N_REWIRES <- ceiling(E * log(1/EPSILON))

for(REWIRE in 1:N_REWIRES){
  e_idxs <- sample(E,2)
  rands <- rbinom(2,1,0.5)
  p_rxn1 <- edges.df$V1[e_idxs[1]]
  f_rxn1 <- edges.df$V2[e_idxs[1]]
  p_rxn2 <- edges.df$V1[e_idxs[2]]
  f_rxn2 <- edges.df$V2[e_idxs[2]]
  if(rands[1] > 0){
    temp <- p_rxn1
    p_rxn1 <- f_rxn2
    f_rxn2 <- temp
  }
  if(rands[2] > 0){
    temp <- p_rxn2
    p_rxn2 <- f_rxn1
    f_rxn1 <- temp
  }
  edges.df[e_idxs[1],] <- c(p_rxn1,f_rxn2)
  edges.df[e_idxs[2],] <- c(p_rxn2,f_rxn1)
  if(mod(REWIRE,10000) == 0){
    print(paste("Pathway hierarchy rewire ",REWIRE," of ",N_REWIRES," complete...",sep=""))
  }
}

write.table(edges.df,
            file=paste(IN_DIR,"rewired_pathway_hierarchy_edges.txt",sep=""),
            row.names = FALSE,
            col.names = FALSE)
