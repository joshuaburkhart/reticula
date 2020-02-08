#!/usr/bin/bash
#
# set max wallclock time hh:mm:ss
#SBATCH --time=12:00:00
#
# set number of tasks
#SBATCH --ntasks=1
#
# set number of cores per node
#SBATCH --cpus-per-task=12
#
# set memory
#SBATCH --mem=64G
#
# set output filename
#SBATCH -o slurm-%j.out-%N
#
# mail all alerts (start, end and abortion)
#SBATCH --mail-type=ALL
#
# mailto
#SBATCH --mail-user=burkhajo@ohsu.edu
#
# set job name
#SBATCH --job-name=pheatmapGen
#
# set queue
#SBATCH -p exacloud
#
# sequentially execute commands in a session on a single node

# run the r script
srun /home/users/burkhajo/R-3.6.2/bin/Rscript ~/WuLab/WuLabLustreDir/reticula/src/pheatmapGen.R
