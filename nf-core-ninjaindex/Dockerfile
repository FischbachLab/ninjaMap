FROM nfcore/base
LABEL authors="Xiandong Meng" \
      description="Docker image containing all requirements for nf-core/ninjaindex pipeline"

# Update conda to latest version.
RUN conda update -n base -c defaults conda

# Install samtools bamtools biopython
RUN conda install -c bioconda -y bamtools samtools
RUN conda install -c conda-forge -y awscli biopython


COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-ninjaindex-1.0dev/bin:$PATH

# Add a script to filter genomes
ADD scripts/genome_filter.py /usr/local/bin/genome_filter.py
RUN chmod +x /usr/local/bin/genome_filter.py
