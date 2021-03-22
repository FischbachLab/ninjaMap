FROM nfcore/base
LABEL authors="Xiandong Meng" \
      description="Docker image containing all requirements for nf-core/ninjaindex pipeline"

ADD environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
RUN /bin/bash -c "source activate nf-core-ninjaindex-1.0dev"
ENV PATH /opt/conda/envs/nf-core-ninjaindex-1.0dev/bin:$PATH

# Update conda to latest version.
#RUN conda update -n nf-core-ninjaindex-1.0dev -c defaults conda

# Install picard samtools bamtools biopython
RUN conda install -c bioconda picard
RUN conda install -c bioconda -y bamtools samtools
RUN conda install -c conda-forge -y awscli biopython


# Add scripts to filter genomes
ADD genome_filter.py /usr/local/bin/genome_filter.py
RUN chmod +x /usr/local/bin/genome_filter.py

ADD ninjaIndex_multiBams.sh /usr/local/bin/ninjaIndex_multiBams.sh
RUN chmod +x /usr/local/bin/ninjaIndex_multiBams.sh

ADD ninjaIndex_multiBams.py /usr/local/bin/ninjaIndex_multiBams.py
RUN chmod +x /usr/local/bin/ninjaIndex_multiBams.py
