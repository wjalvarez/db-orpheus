## Use base image from Databricks
FROM databricksruntime/standard:latest

## Set up arguments
#ARG USER_ID
#ARG GROUP_ID

#RUN addgroup --gid $GROUP_ID user
#RUN adduser --disabled-password --gecos '' --uid $USER_ID --gid $GROUP_ID user

## Set working directory
WORKDIR /home/user/

## Give user writing permissions
#RUN chown $USER_ID:$GROUP_ID /opt/conda/pkgs/urls.txt

##Set path
RUN export PATH=/databricks/conda/bin:$PATH

#RUN /databricks/conda/envs/dcs-minimal/bin/pip install mamba

## Set up conda environment for Snakemake with Python 3.7.3
RUN /databricks/conda/bin/conda install -y -c conda-forge mamba && \
/databricks/conda/bin/mamba create -y -c conda-forge -c bioconda -n snakemake snakemake

## Give conda privileges to user
#RUN chown $USER_ID:$GROUP_ID /opt/conda/envs/snakemake/

## Set user
#USER user

# Make RUN commands use the new environment:
SHELL ["/databricks/conda/bin/conda", "run", "-n", "snakemake", "/bin/bash", "-c"]

#Set path
SHELL ["export", "PATH=/databricks/conda/envs/snakemake/bin:$PATH"]
