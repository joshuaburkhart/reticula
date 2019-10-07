#!/usr/bin/bash
#
# set max wallclock time hh:mm:ss
#SBATCH --time=24:00:00
#
# set number of nodes
#SBATCH --nodes=1
#
# set number of tasks
#SBATCH --ntasks=1
#
# set number of cores per node
#SBATCH --cpus-per-task=24
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
# execute on cluster
srun sg - WuLab
srun sudo /opt/acc/sbin/exadocker load -i /home/users/burkhajo/WuLab/WuLabLustreDir/burkhajo/reticula/bin/savedGct2gctx.tar.gz
srun sudo /opt/acc/sbin/exadocker run burkhajo/gct2gctx --rm=true
