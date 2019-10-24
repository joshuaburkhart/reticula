#!/bin/bash
# I don't see a need to ever chmod+x this file as it should be executed with sbatch.

#add dockerfile to context
cp /Users/burkhajo/PycharmProjects/reticula/src/lib/aim1_genewise/dockerfiles/dockerfile_gct2gctx ./Dockerfile

#add executable to context
cp /Users/burkhajo/.conda/envs/reticula/lib/python3.7/site-packages/cmapPy/pandasGEXpress/gct2gctx.py ./executable.py

#docker login for tagging, save & push
docker login
docker build -t burkhajo/gct2gctx .
docker save burkhajo/gct2gctx | gzip > ./savedGct2gctx.tar.gz

#push to dockerhub
docker push burkhajo/gct2gctx

#tx to Exacloud Lustre fs
scp ./savedGct2gctx.tar.gz burkhajo@exahead1.ohsu.edu:/home/users/burkhajo/WuLab/WuLabLustreDir/reticula/bin/

#remove temporary files
rm ./Dockerfile
rm ./executable.py
rm ./savedGct2gctx.tar.gz
