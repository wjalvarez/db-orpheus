import glob, os
import pandas as pd

configfile: "config/config.local.yaml"

#samples, = glob_wildcards(config['fastqs'] + '/' + '{sample}_1.fq.gz')
sample = config["bam"].rsplit(".",1)[0].rsplit("/",1)[1]
reference="GRCh37.75"

chrs = [str(n) for n in range(1,23)]
chrs.append('X')
chrs.append('Y')
print(chrs)

pairs = [1, 2]
ID = config['ID']
print(sample)

rule all:
	input:
		"/dbfs/db-orpheus/{}/{}.vcf.gz".format(ID, sample)
#		"outs/{}/final/{}.vcf.gz".format(ID, sample)

### include rules ###
#include: 'workflow/rules/align.smk'
#include: 'workflow/rules/qc.smk'
include: 'workflow/rules/call.smk'
