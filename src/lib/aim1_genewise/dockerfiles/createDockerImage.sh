#!/bin/bash
docker login
docker build -t burkhajo/gct2gctx .
docker save burkhajo/gct2gctx | gzip > savedGct2gctx.tar.gz
docker push burkhajo/gct2gctx
scp savedGct2gctx.tar.gz burkhajo@exahead1.ohsu.edu:/home/users/burkhajo/WuLab/WuLabLustreDir/burkhajo/reticula/bin/
