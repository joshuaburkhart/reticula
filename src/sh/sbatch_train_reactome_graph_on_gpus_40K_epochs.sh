#!/bin/bash
#SBATCH --job-name=train_reactome_graph
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jgburk@hawaii.edu
#SBATCH --partition=gpu
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --gres=gpu:8

module purge
module load lang/Python/3.9.5-GCCcore-10.3.0
module load system/CUDA/11.0.2

python ~/zsl_validation/ReactomeGraphClassificationZSLValidationGTEX_Mana.py
