#!/bin/bash

#Go Interactive
#sh ~/go_interactive_gpus.sh

#unzip ./zsl_validation.zip -d zsl_validation/
#cd ~/zsl_validation/

module purge
module load lang/Python/3.9.5-GCCcore-10.3.0
module load system/CUDA/11.0.2
python ReactomeGraphClassificationZSLValidationGTEX_Mana.py

#Editing
#module purge
#module load tools/Vim

#Monitor jobs on hpc
#squeue -o "%all" | head -n1

#Monitor user “anirvan”
#squeue -o "%all" | grep gpu | grep anirvan | head -n1
