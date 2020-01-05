FROM ubuntu:latest
LABEL contact="burkhajo@ohsu.edu"
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 &&\
 echo "deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/" > /etc/apt/sources.list.d/cran.list &&\
 apt-get update &&\
 DEBIAN_FRONTEND=noninteractive apt-get -y install python3-pip r-base tzdata &&\
 pip3 install numpy scipy pandas matplotlib rpy2 &&\
 Rscript -e 'install.packages("BiocManager"); BiocManager::install("DESeq2")' &&\
 mkdir -p /reticula/code /reticula/data
ADD ./executable.py /reticula/code/
ENTRYPOINT exec python3 /reticula/code/executable.py -verbose -file_name ~/WuLab/WuLabLustreDir/reticula/input/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.bz2
