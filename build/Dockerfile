FROM ubuntu:latest
LABEL contact="burkhajo@ohsu.edu"
RUN apt-get update &&\
 apt-get -y install python3-pip r-base &&\
 pip3 install numpy scipy pandas matplotlib rpy2 &&\
 Rscript -e 'install.packages("BiocManager"); BiocManager::install("DESeq2")' &&\
 mkdir -p /reticula/code /reticula/data
ADD ./executable.py /reticula/code/
ENTRYPOINT exec python3 /reticula/code/executable.py -verbose -file_name ~/WuLab/WuLabLustreDir/reticula/input/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.bz2
