/databricks/conda/envs/snakemake/bin/snakemake \
--directory /databricks/db-orpheus \
--snakefile /databricks/db-orpheus/Snakefile \
--configfile /databricks/db-orpheus/config/config.CCLE_full.yaml \
--use-conda --conda-create-envs-only --cores 4
