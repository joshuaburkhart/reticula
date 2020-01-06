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
#SBATCH --mem=12G
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
#SBATCH --job-name=exec_python_script
#
# set queue
#SBATCH -p exacloud
#
# specify node with ssd scratch disk
#SBATCH -C ssdscratch
#
# request x GB scratch disk allocation
#SBATCH --gres disk:16
#
# sequentially execute commands in a session on a single node

# create node-local scratch directory for processed gtex output
srun mkdir -p /mnt/scratch/burkhajo/.data /mnt/scratch/burkhajo/.container

# copy gtex data from Lustre to node-local scratch directory
srun cp ~/WuLab/WuLabLustreDir/reticula/input/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.bz2 /mnt/scratch/burkhajo/.data

#bunzip
srun bzip2 -dc /mnt/scratch/burkhajo/.data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.bz2 >/mnt/scratch/burkhajo/.data/decompressed.dat

# copy docker image from Lustre to node-local scratch directory
srun cp ~/WuLab/WuLabLustreDir/reticula/docker_images/burkhajo_exec_python_script.tar.gz /mnt/scratch/burkhajo/.container/

# load docker image from node-local scratch directory
srun sudo /opt/acc/sbin/exadocker load --input /mnt/scratch/burkhajo/.container/burkhajo_exec_python_script.tar.gz

# run docker image on bound gtex data and then automatically remove container
srun sudo /opt/acc/sbin/exadocker run --rm=true --volume "/mnt/scratch/burkhajo/.data:/reticula/data" burkhajo/exec_python_script

# ensure user images are removed from exacloud
srun sudo /opt/acc/sbin/exadocker rmi "$(sudo /opt/acc/sbin/exadocker images | grep burk | tr -s ' ' | cut -d' ' -f3)"

# archive node-local output
srun tar -cfP node_data.tar /mnt/scratch/burkhajo/.data

# compress node-local output
srun bzip2 -c node_data.tar >node_data.tar.bz2

# copy output to Lustre directory
srun cp node_data.tar.bz2 ~/WuLab/WuLabLustreDir/reticula/input/

# clean node-local scratch directory
srun rm -rf /mnt/scratch/burkhajo
