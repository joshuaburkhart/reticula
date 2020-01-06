FROM ubuntu:18.04
LABEL contact="burkhajo@ohsu.edu"
ENV DEBIAN_FRONTEND noninteractive
# install R 3.6.1, which was designed by satan
RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
	    apt-utils \
		ed \
		less \
		locales \
		vim-tiny \
		wget \
		ca-certificates \
		apt-transport-https \
		gsfonts \
		gnupg2 \
	&& rm -rf /var/lib/apt/lists/*
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
	&& locale-gen en_US.utf8 \
	&& /usr/sbin/update-locale LANG=en_US.UTF-8
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8
RUN echo "deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/" > /etc/apt/sources.list.d/cran.list
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
ENV R_BASE_VERSION 3.6.1
RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
		littler \
        r-cran-littler \
		r-base=${R_BASE_VERSION}* \
		r-base-dev=${R_BASE_VERSION}* \
		r-recommended=${R_BASE_VERSION}* \
        && echo 'options(repos = c(CRAN = "https://cloud.r-project.org/"), download.file.method = "libcurl")' >> /etc/R/Rprofile.site \
        && echo 'source("/etc/R/Rprofile.site")' >> /etc/littler.r \
	&& ln -s /usr/share/doc/littler/examples/install.r /usr/local/bin/install.r \
	&& ln -s /usr/share/doc/littler/examples/install2.r /usr/local/bin/install2.r \
	&& ln -s /usr/share/doc/littler/examples/installGithub.r /usr/local/bin/installGithub.r \
	&& ln -s /usr/share/doc/littler/examples/testInstalled.r /usr/local/bin/testInstalled.r \
	&& install.r docopt \
	&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds \
	&& rm -rf /var/lib/apt/lists/*

# install python3 and other components
RUN apt-get update &&\
 apt-get update &&\
 apt-get -yf install python3-pip libxml2-dev &&\
 pip3 install numpy scipy pandas matplotlib rpy2 &&\
 Rscript -e 'install.packages("BiocManager"); BiocManager::install("DESeq2")' &&\
 mkdir -p /reticula/code /reticula/data
ADD ./executable.py /reticula/code/
ENTRYPOINT exec python3 /reticula/code/executable.py -verbose -file_name ~/WuLab/WuLabLustreDir/reticula/input/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.bz2
