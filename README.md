```
docker run -it --rm \
    --user "$(id -u):$(id -g)" \
    -v /data/exploratory/Users/jeff.alvarez/orpheus:/home/user/analysis \
    -v /data/exploratory/Users/jeff.alvarez/orpheus/data/samples/single:/home/user/input \
    -v /data/exploratory/Users/jeff.alvarez/orpheus/data:/home/user/ref \
    orpheus:0.0 /bin/bash -c \
    "conda run -n snakemake \
    snakemake -j 6 --keep-remote --use-conda \
    --directory /home/user/analysis \
    --configfile /home/user/analysis/config/config.docker.yaml \
    -s /home/user/analysis/Snakefile -np"
```
