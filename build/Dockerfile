##TODO: pull a smaller image after all dependencies are known
FROM ubuntu:latest

# set contact email
LABEL contact="burkhajo@ohsu.edu"

# update, add deps &
RUN apt-get update &&\
 apt-get -y install python3-pip &&\
 pip3 install numpy scipy pandas matplotlib cmapPy &&\
 mkdir -p /reticula/code /reticula/data

# add executible python script to container /reticula/code directory
ADD ./executable.py /reticula/code/

# execute python script on node-local data
ENTRYPOINT exec python3 /reticula/code/executable.py -verbose \
 -filename reticula/data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct \
 -output_filepath reticula/data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gctx
