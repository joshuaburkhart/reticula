#!/bin/bash

mkdir /tmp/zsl_validation
cd /tmp/zsl_validation/
cp ~/PycharmProjects/reticula/data/gtex/input/edges.txt ./
cp ~/PycharmProjects/reticula/data/gtex/input/zsl_gtex_* ./
cp ~/PycharmProjects/reticula/data/tcga/input/zsl_tcga_* ./
cp ~/PycharmProjects/reticula/src/sh/prepare_mana_env.sh ./
cp ~/PycharmProjects/reticula/src/python/drivers/zsl_validation/gtex/ReactomeGraphClassificationZSLValidationGTEX_Mana.py ./
cp ~/PycharmProjects/reticula/src/python/drivers/zsl_validation/tcga/ReactomeGraphClassificationZSLValidationTCGA_Mana.py ./
cd ../
zip -rj zsl_validation.zip zsl_validation
. ~/.profile
echo "Initiating transfer to Mana HPC. Authorization Required."
huscp zsl_validation.zip ~/
echo "Transfer complete."