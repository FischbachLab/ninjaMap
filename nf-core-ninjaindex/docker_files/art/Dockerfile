FROM nfcore/base

LABEL authors="Xiandong Meng" \
      description="Docker image containing all requirements for nf-core/ninjaindex pipeline" \
      software="art" \
      about.summary="ART is a set of simulation tools to generate synthetic next-generation sequencing reads." \
      about.home="https://www.niehs.nih.gov/research/resources/software/biostatistics/art/" \
      software.version="0.5.4-1-deb"


USER root:root

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y awscli pigz vim htop && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/*

RUN conda install -c bioconda -y picard curl
RUN conda clean -ya

RUN mkdir -p /mnt
WORKDIR /mnt

#RUN wget https://www.niehs.nih.gov/research/resources/assets/docs/artbinmountrainier2016.06.05linux64.tgz
#RUN tar zxf artbinmountrainier2016.06.05linux64.tgz

RUN curl https://www.niehs.nih.gov/research/resources/assets/docs/artbinmountrainier2016.06.05linux64.tgz --output - | tar xz
WORKDIR /mnt

ADD run_art.sh /usr/local/bin/run_art.sh
RUN chmod +x /usr/local/bin/run_art.sh
