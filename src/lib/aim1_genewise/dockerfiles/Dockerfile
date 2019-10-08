FROM debian:buster-slim
LABEL contact="burkhajo@ohsu.edu"
VOLUME ["/home/users/burkhajo/WuLab/WuLabLustreDir/burkhajo"]
RUN apt-get update && apt-get install -y python-pip && pip install --verbose pandas cmapPy
CMD ["gct2gctx.py","-f GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct","-o GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gctx"]
