# Use base image for miniconda3
FROM continuumio/miniconda3

# Set up arguments
ARG USER_ID
ARG GROUP_ID

RUN addgroup --gid $GROUP_ID user
RUN adduser --disabled-password --gecos '' --uid $USER_ID --gid $GROUP_ID user

# Set working directory
WORKDIR /home/user/

# Give user writing permissions
#RUN chown $USER_ID:$GROUP_ID /opt/conda/pkgs/urls.txt

# Set up conda environment for Snakemake with Python 3.7.3
RUN conda install -y -c conda-forge mamba && \
mamba create -y -c conda-forge -c bioconda -n snakemake snakemake

# Give conda privileges to user
RUN chown $USER_ID:$GROUP_ID /opt/conda/envs/snakemake/

# Set user
USER user

# Make RUN commands use the new environment:
SHELL ["conda", "run", "-n", "snakemake", "/bin/bash", "-c"]

#Set path
SHELL ["export", "PATH=/opt/conda/envs/snakemake/bin:$PATH"]
