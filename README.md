```
docker run -it --rm \
    -v /data/exploratory/Users/jeff.alvarez/omics-pipeline:/analysis \
    -v /data/exploratory/Users/jeff.alvarez/omics-pipeline/data/samples/single:/input \
    omics omics-pipeline:1.0 \
    snakemake -j 6 --timestamp --verbose -np --keep-remote \
    -s /analysis/Snakefile --directory /analysis
```
