# reticula
## Graphical Abstract
![Graphical Abstract Image](Graphical%20Abstract.png)
## Program Flowchart
![Program Flowchart Image](Program%20Flowchart.png)

This repository accompanies our article: "Biology Inspired Graph Neural Network Encodes Reactome and Reveals Biochemical Reactions of Disease".

The above graphical abstract represents our information processing workflow and the program flowchart indicates which files are used for each phase of analysis where colors matched to represent each dataset.

## Data Files

### rse_gene.Rdata
Provided by Recount2 for each dataset.

### Ensembl2ReactomeReactions.txt
Provided by Reactome.

### ReactionNetwork_Rel_71_122820.txt
Provided by Reactome.

## Source Code Files

### final_dataset_processing_a.R
Requires rse_gene.Rdata and Ensembl2ReactomeReactions.txt.

### final_dataset_processing_b.R

### final_dataset_processing_c.R

### pca_and_knn_calculation.R

### prepare_data_for_pytorch_reaction_network.R
Requires ReactionNetwork_Rel_71_122820.txt

### ReactomeGraphClassification.ipynb (GTEx-specific)

### ReactomeResNetClassification.ipynb (GTEx-specific)

### ReactomeGraphClassificationValidationTCGA_local.py (TCGA-specific)

### ReactomeResnetClassificationValidationTCGA_local.py (TCGA-specific)

### ReactomeGraphClassificationValidationSRP035988_local.py (SRP035988-specific)

### ExtractModelWeightsSRP035988_local.py (SRP035988-specific)

