import glob, os
import pandas as pd

configfile: "config/config.local.yaml"

sample = config["bam"].rsplit(".",1)[0].rsplit("/",1)[1]
reference="GRCh37.75"

chr_dict = {"chunk_01": "-L 1", "chunk_02": "-L 2", "chunk_03": "-L 3", "chunk_04": "-L 4",
	"chunk_05": "-L 5", "chunk_06": "-L 6", "chunk_07": "-L 7", "chunk_08": "-L 8",
	"chunk_09": "-L 9", "chunk_10": "-L 10", "chunk_11": "-L 12", "chunk_12": "-L 12 -L 13",
	"chunk_13": "-L 14 -L 15", "chunk_14": "-L 16 -L 17 -L 18", "chunk_15": "-L 19 -L 20 -L 21 -L 22", 
	"chunk_16": "-L X -L Y"}
chr_chunks = list(chr_dict.keys())
chr_params = list(chr_dict.values())

ID = config['ID']

rule all:
	input:
		"/dbfs/db-orpheus/{}/{}.vcf.gz".format(ID, sample)
#		"outs/{}/final/{}.vcf.gz".format(ID, sample)

### include rules ###
#include: 'workflow/rules/align.smk'
#include: 'workflow/rules/qc.smk'
include: 'workflow/rules/call.db.smk'
