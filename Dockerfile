############################################################
# Dockerfile to build HLA typing tools
# Based on Ubuntu 16.04
############################################################
# Set the base image to Ubuntu
FROM ubuntu:16.04
# File Author / Maintainer
LABEL maintainer="Ruth Nanjala rnanjala@icipe.org"
LABEL description="This is custom Docker Image for \
    HLA typing tools and their dependencies."
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH
################## BEGIN INSTALLATION ######################
# Install wget
RUN apt-get update && apt-get install -y \
    autoconf \
    build-essential \
    git \
    libncurses5-dev \
    pkg-config \
    unzip \
    wget curl \
    python python-dev \
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev &&\
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
# Install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
    tar -xvf samtools-1.9.tar.bz2 && \
    cd samtools-1.9 && \
    ./configure --prefix=/usr/local && \
    make && \
    make install && \
    cd .. && rm -rf samtools-1.9*
# Install VCFTools
RUN wget https://github.com/vcftools/vcftools/releases/download/v0.1.16/vcftools-0.1.16.tar.gz && \
    tar -xvf vcftools-0.1.16.tar.gz && \
    cd vcftools-0.1.16 && \
    ./configure && \
    make && \
    make install
# Install bcftools
RUN wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 && \
    tar -xvf bcftools-1.9.tar.bz2 && \
    cd bcftools-1.9 && \
    ./configure --prefix=/usr/local && \
    make && \
    make install && \
    cd .. && rm -rf bcftools-1.9*
WORKDIR /usr/local/bin/
# Install HLA-VBSEQ
RUN wget http://nagasakilab.csml.org/hla/HLAVBSeq.jar && \
    wget http://nagasakilab.csml.org/hla/bamNameIndex.jar && \
    wget http://nagasakilab.csml.org/hla/parse_result.pl && \
    wget http://nagasakilab.csml.org/hla/call_hla_digits.py && \
    wget http://sourceforge.net/projects/picard/files/picard-tools/1.119/picard-tools-1.119.zip && \
        unzip picard-tools-1.119.zip && \
        mv -f picard-tools-1.119/SamToFastq.jar . && \
        rm -fr picard-tools-1.119.zip picard-tools-1.119
# Install miniconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    echo "conda activate base" >> ~/.bashrc
#install csh
RUN conda clean --all --yes && \
    conda install -c conda-forge tcsh
#install svn  
RUN conda clean --all --yes && \
    conda install -c anaconda svn
#install nano
RUN conda clean --all --yes && \
    conda install -c conda-forge nano
#install java
RUN conda clean --all --yes && \
    conda install -c bioconda java-jdk 
#python 3.6
RUN conda clean --all --yes && \
    conda install -c anaconda python=3.6
# pandas 
RUN conda clean --all --yes && \
    conda install -c anaconda pandas
#perl 5.26.2
RUN conda clean --all --yes && \
    conda install -c anaconda perl=5.26.2
#minimap2
RUN conda clean --all --yes && \
    conda install -c bioconda minimap2

RUN useradd -m docker && echo "docker:docker" | chpasswd && adduser docker sudo

# RUN useradd --create-home --shell /bin/bash ubuntu && \
#   chown -R ubuntu:ubuntu /home/ubuntu
# USER ubuntu
USER docker
# CMD /bin/bash
CMD ["/bin/bash","-i"]
