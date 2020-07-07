<h1>Omics Pipeline</h1>
This workflow performs variant calling and expression quantification with
STAR, GATK, and RSEM.
<h2>Usage</h2>
<h3>Step 0: Install Snakemake</h3>
<b>IMPORTANT</b>: Omics Pipeline was written with [Snakemake](https://snakemake.readthedocs.io/en/stable/),
a Python-based workflow manager. The recommended way to install Snakemake is via
Conda, which is not installed system-wide on HPC at the time of writing. Therefore,
<b>the easiest way to start working with Omics Pipeline is to use the Docker method
described in</b> [Step 3](#step-3-execute-workflow), 
which comes installed with Snakemake and associated dependencies.

<h3>Step 1: Get workflow</h3>
[Clone](https://help.github.com/en/github/creating-cloning-and-archiving-repositories/cloning-a-repository)
the repo to the desired working directory for the concrete project/run on
your machine:

```
git clone https://git.agios.local/Jeff.Alvarez/omics-pipeline.git
```
<h3>Step 2: Configure workflow</h3>
The workflow can be configured with the config file located in the 
<code>config/config.test.yaml</code> file. The layout is based on the config
file used for the [Array Studio RNA-seq pipeline](https://git.agios.local/Mark.Fletcher/array_studio_RNAseq_pipeline):

```
ID: 2020_04_30_test
Title: "2020-04-30 Test"
HPC_ID: jeff.alvarez
Contact_name: "Jeff Alvarez"
Organism: "Human"
fastqs: /data/exploratory/Users/jeff.alvarez/omics-pipeline/data/samples/single/
ref:
        fa: /data/exploratory/Users/jeff.alvarez/omics-pipeline/data/Human_B37.3_chr1.fasta
        gtf: /data/exploratory/Users/jeff.alvarez/omics-pipeline/data/Human_B37.3_chr1.gtf
        build: Human_B37.3
trimming:
        skip: false
known_sites: /data/exploratory/Users/jeff.alvarez/pipeline_ins/dbsnp_138.b37.chr_1.vcf.gz
```
The config file takes the following values:
* <b>ID</b>: Name of the pipeline ID--output directories will take this name.
* <b>Title</b>: String character for analysis ID.
* <b>HPC_ID</b>: User ID on HPC.
* <b>Contact_name</b>: String character of user ID.
* <b>Organism</b>: "Human" or "Mouse".
* <b>fastqs</b>: Absolute path to input FASTQ files for pipeline (.fq.gz).
* <b>ref</b>: Reference data to be used in alignment and variant calling.
     - <b>fa</b>: Absolute path to fasta file (.fasta).
     - <b>gtf</b>: Absolute path to gene annotation file (.gtf).
     - <b>build</b>: Name of reference genome--output files will use this prefix.
* <b>trimming</b>: Logical to skip trimming step [[to-do]].
* <b>known_sites</b>: Absolute path to known variants file (.vcf.gz).

<h3>Step 3: Execute workflow</h3>
<h4>Docker</h4>
The workflow may also be deployed as a Docker image, where a conda environment
is set up with Snakemake and dependencies installed. When run, the conda
environment is activated then the Snakemake directory and input files are
mounted as <code>analysis</code> and <code>input</code> respectively:

```
docker run -it --rm \
    -u $(id -u):$(id -g) \
    -v /data/exploratory/Users/jeff.alvarez/omics-pipeline:/analysis \
    -v /data/exploratory/Users/jeff.alvarez/omics-pipeline/data/samples/single:/input \
    omics-pipeline:1.0 /bin/bash -c \
    "conda run -n snakemake \
    snakemake -j 6 --keep-remote --use-conda \
    --directory /analysis \
    --configfile /analysis/config/config.docker.yaml \
    -s /analysis/Snakefile -np"
```

<h4>Local</h4>
Test your configuration by performing a dry-run:

```
snakemake --use-conda -np
```

Execute the workflow locally using <code>$N</code> cores:

```
snakemake --use-conda --cluster --cores $N
```
