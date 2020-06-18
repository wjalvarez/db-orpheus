# Use base image for miniconda3
FROM continuumio/miniconda3

# Set up conda environment for Snakemake with Python 3.7.3
RUN conda install -y -c conda-forge mamba && \
mamba create -y -c conda-forge -c bioconda -n snakemake snakemake

# Make RUN commands use the new environment:
SHELL ["conda", "run", "-n", "snakemake", "/bin/bash", "-c"]

#Set path
SHELL ["export", "PATH=/opt/conda/envs/snakemake/bin:$PATH"]

# Set working directory
WORKDIR /analysis
