#!/bin/bash
echo "initiating asynchronous docker image construction..."
echo "copying /home/burkhart/PycharmProjects/reticula/src/lib/drivers/aim1_genewise/process_w_deseq2.py to ./executable.py..."
cp /home/burkhart/PycharmProjects/reticula/src/lib/drivers/aim1_genewise/process_w_deseq2.py ./executable.py
echo "logging into docker with dockerNmeSpc == burkhajo/exec_python_script, dockerImgFilNme == burkhajo_exec_python_script.tar.gz..."
docker login
docker build -t burkhajo/exec_python_script .
docker save burkhajo/exec_python_script | gzip >./burkhajo_exec_python_script.tar.gz
docker push burkhajo/exec_python_script
echo "copying burkhajo_exec_python_script.tar.gz to burkhajo@exahead1.ohsu.edu:~/WuLab/WuLabLustreDir/reticula/docker_images/burkhajo_exec_python_script.tar.gz..."
scp ./burkhajo_exec_python_script.tar.gz burkhajo@exahead1.ohsu.edu:~/WuLab/WuLabLustreDir/reticula/docker_images/burkhajo_exec_python_script.tar.gz
echo "copying srun_docker_container.sh to burkhajo@exahead1.ohsu.edu:~/WuLab/WuLabLustreDir/reticula/docker_images/burkhajo_exec_python_script.tar.gz..."
scp ./srun_docker_container.sh burkhajo@exahead1.ohsu.edu:~/WuLab/WuLabLustreDir/reticula/slurm_scripts/srun_docker_container.sh

#rm ./Dockerfile
#rm ./executable.py
#rm ./burkhajo_exec_python_script.tar.gz
