# reticula
## Graphical Abstract
![Graphical Abstract Image](Graphical%20Abstract.png)
## Program Flowchart
![Program Flowchart Image](Program%20Flowchart.png)

This repository accompanies our article: "Biology Inspired Graph Neural Network Encodes Reactome and Reveals Biochemical Reactions of Disease".

The above graphical abstract represents our information processing workflow and the program flowchart indicates which files are used for each phase of analysis. Colors are matched to represent each dataset.

## Data Files

### rse_gene.Rdata
Provided by Recount2 and unique to each dataset.

### Ensembl2ReactomeReactions.txt
Provided by Reactome and shared across datasets.

### ReactionNetwork_Rel_71_122820.txt
Provided by Reactome and shared across datasets.

## Source Code Files

### final_dataset_processing_a.R
#### Requires rse_gene.Rdata and Ensembl2ReactomeReactions.txt.
This file is copied with hardcoded values for each dataset:  
GTEx: reticula/src/r/process_gtex_data/  
TCGA: reticula/src/r/process_tcga_data/  
SRP035988: reticula/src/r/process_srp035988_data/  

This script parses the rse_gene.Rdata file, storing the tissue phenotype labels and filtering by identifiers found in Ensembl2ReactomeReactions.txt.  

### final_dataset_processing_b.R
This file is copied with hardcoded values for each dataset:  
GTEx: reticula/src/r/process_gtex_data/  
TCGA: reticula/src/r/process_tcga_data/  
SRP035988: reticula/src/r/process_srp035988_data/  

This script applies a scale factor and small pseudocount to an R dataframe of transcript counts and transforms it into an object compatible with DESeq2.  

### final_dataset_processing_c.R
This file is copied with hardcoded values for each dataset:  
GTEx: reticula/src/r/process_gtex_data/  
TCGA: reticula/src/r/process_tcga_data/  
SRP035988: reticula/src/r/process_srp035988_data/  

This script applies the DESeq2::vst() "Variance Stabalizing Transform" to the expression data object resulting in normalized counts. Additionally, reaction <-> gene identifier files are generated for annotation and debugging purposes.  

### pca_and_knn_calculation.R
This file is copied with hardcoded values for each dataset:  
GTEx: reticula/src/r/analysis/  
TCGA: reticula/src/r/validation/tcga/  
SRP035988: reticula/src/r/validation/srp035988/  

This script performs principal component analysis and stores the first principal component values across reactions. Additionally, k-nearest neighbor classification and plotting is performed for debugging purposes.  

### prepare_data_for_pytorch_reaction_network.R
#### Requires ReactionNetwork_Rel_71_122820.txt
This file is copied with hardcoded values for each dataset:  
GTEx: reticula/src/r/analysis/  
TCGA: reticula/src/r/validation/tcga/  
SRP035988: reticula/src/r/validation/srp035988/  

This script filters Reactome reaction network edges to match retained identifiers and stores edge, feature and target data in a form appropriate for loading into pytorch geometric.  

### ReactomeGraphClassification.ipynb (GTEx-specific)

### ReactomeResNetClassification.ipynb (GTEx-specific)

### ReactomeGraphClassificationValidationTCGA_local.py (TCGA-specific)

### ReactomeResnetClassificationValidationTCGA_local.py (TCGA-specific)

### ReactomeGraphClassificationValidationSRP035988_local.py (SRP035988-specific)

### ExtractModelWeightsSRP035988_local.py (SRP035988-specific)

