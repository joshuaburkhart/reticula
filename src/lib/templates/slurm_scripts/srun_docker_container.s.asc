#!/usr/bin/bash
#
# set max wallclock time hh:mm:ss
#SBATCH ${sbatchTime}
#
# set number of tasks
#SBATCH ${sbatchNTasks}
#
# set number of cores per node
#SBATCH ${sbatchCpusPerTask}
#
# set memory
#SBATCH ${sbatchMemory}
#
# set output filename
#SBATCH ${sbatchOutputFilename}
#
# mail all alerts (start, end and abortion)
#SBATCH ${sbatchMailAlerts}
#
# mailto
#SBATCH ${sbatchMailTo}
#
# set job name
#SBATCH ${sbatchJobName}
#
# set queue
#SBATCH ${sbatchQueue}
#
# specify node with ssd scratch disk
#SBATCH ${sbatchSSD}
#
# request x GB scratch disk allocation
#SBATCH ${sbatchScratchStorage}
#
# sequentially execute commands in a session on a single node

# create node-local scratch directory for processed gtex output
srun mkdir -p ${scratchDataDir} ${scratchContDir}

# copy gtex data from Lustre to node-local scratch directory
srun cp ${lusRelDatPath} ${scratchDataDir}

# copy gtex design matrix from Lustre to node-local scratch directory
srun cp ${lusRelDsnPath} ${scratchDataDir}

#bunzip data and design matrix
srun bzip2 -dc ${scratchDataDir}/${datFilNme} >${scratchDataDir}/decompressed.dat
srun bzip2 -dc ${scratchDataDir}/${dsnFilNme} >${scratchDataDir}/decompressed.dsn

# copy docker image from Lustre to node-local scratch directory
srun cp ${lusRelImgPath} ${scratchContDir}/

# load docker image from node-local scratch directory
srun ${exadocker} load --input ${scratchContDir}/${dockerImgFilNme}

# run docker image on bound gtex data and then automatically remove container
srun ${exadocker} run --rm=true --volume "${scratchDataDir}:${dockerDataDir}" ${dockerNmeSpc}

# ensure user images are removed from exacloud
srun ${exadocker} rmi "$(${exadocker} images | grep burk | tr -s ' ' | cut -d' ' -f3)"

# archive, compress and transfer node-local output to Lustre directory
srun tar cfjS ${lusRelDatDir}/node_data.tar.bz2 ${scratchDataDir}

# clean node-local scratch directory
srun rm -rf ${scratchHomeDir}
