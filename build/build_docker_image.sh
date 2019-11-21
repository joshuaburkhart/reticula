#!/bin/bash
cp ../src/lib/templates/dockerfiles/.d.asc ./Dockerfile
cp /Users/burkhajo/PycharmProjects/reticula/src/lib/drivers/aim1_genewise/process_w_deseq2.py ./executable.py
docker login
docker build_scripts -t burkhajo/ .
docker save burkhajo/ | gzip > ./.tar.gz
docker push burkhajo/
scp ./.tar.gz burkhajo@exahead1.ohsu.edu:~/WuLab/WuLabLustreDir/reticula/bin/
rm ./Dockerfile
rm ./executable.py
rm ./.tar.gz
