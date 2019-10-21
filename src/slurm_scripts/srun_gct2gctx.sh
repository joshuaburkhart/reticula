#!/usr/bin/bash
#
# set max wallclock time hh:mm:ss
#SBATCH --time=24:00:00
#
# set number of tasks
#SBATCH --ntasks=1
#
# set number of cores per node
#SBATCH --cpus-per-task=24
#
# set memory
#SBATCH --mem=32G
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
#SBATCH --job-name=gct2gctx
#
# set queue
#SBATCH -p exacloud
#
# specify node with ssd scratch disk
#SBATCH -C ssdscratch
#
# request 8GB scratch disk allocation
#SBATCH --gres disk:16G
#
# sequentially execute commands in a session on a single node

# create node-local scratch directory for processed gtex output
srun mkdir -p /mnt/scratch/burkhajo/.data /mnt/scratch/burkhajo/.container

# copy gtex data from Lustre to node-local scratch directory
srun cp ~/WuLab/WuLabLustreDir/reticula/data/input/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.bz2 \
 /mnt/scratch/burkhajo/.data/

# copy docker image from Lustre to node-local scratch directory
srun cp ~/WuLab/WuLabLustreDir/reticula/bin/savedGct2gctx.tar.gz \
 /mnt/scratch/burkhajo/.container/

# load docker image from node-local scratch directory
srun sudo /opt/acc/sbin/exadocker load --input /mnt/scratch/burkhajo/.container/savedGct2gctx.tar.gz

# run docker image on bound gtex data and then automatically remove container
#TODO: figure out how to set /reticula/data more elegantly--currently specified by dockerfile
srun sudo /opt/acc/sbin/exadocker run --rm=true --volume "/mnt/scratch/burkhajo/.data:/reticula/data" burkhajo/gct2gctx

# ensure user images are removed from exacloud
srun sudo /opt/acc/sbin/exadocker rmi "$(sudo /opt/acc/sbin/exadocker images | grep burk | tr -s ' ' | cut -d' ' -f3)"

# remove gtex data from node-local scratch directory
srun rm -f /mnt/scratch/burkhajo/.data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.bz2

# remove docker image from node-local scratch directory
srun rm -f /mnt/scratch/burkhajo/.container/savedGct2gctx.tar.gz

# copy remainder of node-local scratch directory to Lustre directory
srun cp /mnt/scratch/burkhajo/* ~/WuLab/WuLabLustreDir/reticula/data/input

# clean node-local scratch directory
srun rm -f /tmp/burkhajo/*
srun rmdir /tmp/burkhajo