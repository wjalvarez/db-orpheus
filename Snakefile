import glob, os
import pandas as pd

configfile: "config/config.databricks.yaml"

#samples, = glob_wildcards(config['fastqs'] + '/' + '{sample}_1.fq.gz')
sample = config["bam"].rsplit(".",1)[0].rsplit("/",1)[1]
#samples, = glob_wildcards(config['bam'])

pairs = [1, 2]
ID = config['ID']
print(sample)

rule all:
	input:
#		"outs/calls/{}.vcf.gz".format(sample)
		"/dbfs/db-orpheus/{}.vcf.gz".format(sample)

### include rules ###
#include: 'workflow/rules/align.smk'
#include: 'workflow/rules/qc.smk'
include: 'workflow/rules/call.smk'

#rule raw_counts:
#	input:
#		gtf = config['ref']['gtf'],
#		bams = expand('outs/STAR/bams/{sample}.Aligned.sortedByCoord.out.bam')
#	output:
#		'outs/counts/Pipeline.Counts.tsv'
#	threads: 
#		2
#	conda: 
#		'workflow/envs/raw_counts.yaml'
#	log: 
#		'logs/raw_counts.log'
#	script:
#		'workflow/scripts/create_counts.R'

