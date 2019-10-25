#!/bin/bash
# I don't see a need to ever chmod+x this file as it should be executed with sbatch.

DOCKERFILE_PATH='../src/lib/aim1_genewise/dockerfiles/dockerfile_gct2gctx'

#PYTHON3_PY_PATH='/Users/burkhajo/.conda/envs/reticula/lib/python3.7/site-packages/cmapPy/pandasGEXpress/gct2gctx.py'
PYTHON3_PY_PATH='/home/burkhart/miniconda3/lib/python3.7/site-packages/cmapPy/pandasGEXpress


DOCKER_ACCOUNT='burkhajo'

DOCKER_CNT_IMG='gct2gctx'

DOCKER_NAMESPC="$DOCKER_ACCOUNT/$DOCKER_CNT_IMG"

EXACLD_ACCOUNT='burkhajo'

LUSTRE_CNT_STR='/home/users/burkhajo/WuLab/WuLabLustreDir/reticula/bin/'

#add dockerfile to context
cp $DOCKERFILE_PATH ./Dockerfile

#add executable to context
cp $PYTHON3_PY_PATH ./executable.py

#docker login for tagging, save & push
docker login
docker build -t $DOCKER_NAMESPC .
docker save $DOCKER_NAMESPC | gzip > ./saved$DOCKER_CNT_IMG.tar.gz

#push to dockerhub
docker push $DOCKER_NAMESPC

#tx to Exacloud Lustre fs
scp ./saved$DOCKER_CNT_IMG.tar.gz $EXACLD_ACCOUNT@exahead1.ohsu.edu:$LUSTRE_CNT_STR
#remove temporary files
rm ./Dockerfile
rm ./executable.py
rm ./saved$DOCKER_CNT_IMG.tar.gz
