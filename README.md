# reticula [![DOI](https://zenodo.org/badge/212217385.svg)](https://zenodo.org/badge/latestdoi/212217385)
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
This file is located in reticula/src/python/drivers/aim2/

This python notebook trains the graph neural network model on the GTEx dataset for 500 epochs and stores it. Additionally, cross-validation, batch-size evaluation and figure generation logic is included for benchmarking and debugging purposes.  

### ReactomeResNetClassification.ipynb (GTEx-specific)
This file is located in reticula/src/python/drivers/aim2/

This python notebook trains the Resnet18 model on the GTEx dataset for 500 epochs and stores it. Additionally, cross-validation logic is included for benchmarking purposes.  

### ReactomeGraphClassificationValidationTCGA_local.py (TCGA-specific)
This file is located in reticula/src/python/drivers/validation/tcga/

This python script loads the graph neural network model trained on the GTEx dataset, tunes it on half of the TCGA dataset for 500 epochs and tests it on the other half, calculating the ARI.  

### ReactomeResnetClassificationValidationTCGA_local.py (TCGA-specific)
This file is located in reticula/src/python/drivers/validation/tcga/

This python script loads the Resnet18 model trained on the GTEx dataset, tunes it on half of the TCGA dataset for 500 epochs and tests it on the other half, calculating the ARI.  

### ReactomeGraphClassificationValidationSRP035988_local.py (SRP035988-specific)
This file is located in reticula/src/python/drivers/validation/srp035988/

This python script loads the graph neural network model trained on the GTEx dataset, tunes it on the SRP035988 dataset for 500 epochs and stores it.  

### ExtractModelWeightsSRP035988_local.py (SRP035988-specific)
This file is located in reticula/src/python/drivers/validation/srp035988/

This python script calculates the reaction <-> reaction edge weight of the graph neural network trained on the GTEx dataset and tuned on the SRP035988 dataset and stores it.  
