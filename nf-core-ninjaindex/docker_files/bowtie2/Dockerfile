# Base Image
FROM nfcore/base
LABEL authors="Xiandong Meng" \
      description="Docker image containing all requirements for nf-core/ninjaindex pipeline" \
      software.name="Bowtie2" \
      software.description="An ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences" \
      software.website="http://bowtie-bio.sourceforge.net/bowtie2/index.shtml" \
      software.documentation="http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml" \
      container.category="tool"

# Update conda to latest version.
RUN conda update -n base -c defaults conda

# Install Bowtie2
RUN conda install -c bioconda -y bowtie2 samtools openssl=1.0 bbmap
RUN conda install -c conda-forge -y awscli

ADD run_bowtie2.sh /usr/local/bin/run_bowtie2.sh
ADD run_bowtie2_origin.sh /usr/local/bin/run_bowtie2_origin.sh
ADD run_bowtie2_header.sh /usr/local/bin/run_bowtie2_header.sh
RUN chmod +x /usr/local/bin/run_bowtie2.sh
RUN chmod +x /usr/local/bin/run_bowtie2_origin.sh
RUN chmod +x /usr/local/bin/run_bowtie2_header.sh
