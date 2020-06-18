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
