<h1>Omics Pipeline</h1>
This workflow performs variant calling and expression quantification with
STAR, GATK, and RSEM.
<h2>Usage</h2>
<h3>Step 1: Install workflow</h3>
If you simply want to use this workflow, download and extract the latest
release. The easiest way is to [clone](https://help.github.com/en/github/creating-cloning-and-archiving-repositories/cloning-a-repository)
the repo to the desired working directory for the concrete project/run on
your machine:

```
git clone https://git.agios.local/Jeff.Alvarez/omics-pipeline.git
```
<h2>Docker</h2>
Docker can be run with the following code:

```
docker run -it --rm \
    -v /data/exploratory/Users/jeff.alvarez/omics-pipeline:/analysis \
    -v /data/exploratory/Users/jeff.alvarez/omics-pipeline/data/samples/single:/input \
    omics-pipeline:1.0 /bin/bash -c \
    "conda run -n snakemake \
    snakemake -j 6 --keep-remote --use-conda \
    --directory /analysis \
    --configfile /analysis/config/config.docker.yaml \
    -s /analysis/Snakefile"
```
