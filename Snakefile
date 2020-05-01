import glob, os

configfile: "configB16.yaml"

samples, = glob_wildcards(config['fastqs'] + '/' + '{sample}_1.fq.gz')
pairs = [1, 2]

rule all:
	input:
		expand('{wd}/outs/STAR_index/{build}', wd = config['wd'], build = config['ref']['build']),
		expand('outs/STAR/{sample}_pass1/SJ.out.tab', sample = samples),
		expand('outs/STAR/{sample}_pass1/{sample}_Pass1SJ.filtered.tab', sample = samples),
		expand('outs/STAR/bams/{sample}.Aligned.sortedByCoord.out.bam', sample = samples),
		'outs/counts/Pipeline.Counts.tsv',
		'outs/qc/multiqc_report.html'

### include rules ###
include: 'workflow/rules/qc.smk'

rule star_index:
	input:
		fa = config['ref']['fa'],
		gtf = config['ref']['gtf']
	output:
		directory('{wd}/outs/STAR_index/{build}')
	threads: 8
	shell:
		'rm -rf {output} && '
		'mkdir {output} && '
		'STAR --runThreadN {threads} '
		'--runMode genomeGenerate '
		'--genomeDir {output} '
		'--genomeFastaFiles {input.fa} '
		'--sjdbGTFfile {input.gtf} '
		'--sjdbOverhang 100'

rule star_pass1:
	input:
		refdir = expand('{wd}/outs/STAR_index/{build}', wd = config['wd'], build = config['ref']['build'])
#		refdir = expand('{wd}/outs/STAR_index/{build}', wd = config['wd'], build = config['ref']['build'])
	params:
		outdir = 'outs/STAR/{sample}_pass1',
		rmbam = config['wd'] + '/' + 'outs/STAR/{sample}_pass1/Aligned.out.bam',
		ID = config['fastqs'] + '/' + '{sample}',
		wd = config['wd']
	output:
		'outs/STAR/{sample}_pass1/SJ.out.tab'
	threads: 4
	shell:
		'rm -rf {params.outdir} && '
		'mkdir {params.outdir} && '
		'cd {params.outdir} && '
		'STAR --runThreadN {threads} '
		'--genomeDir {input.refdir} '
		'--readFilesIn {params.ID}_1.fq.gz {params.ID}_2.fq.gz '
		'--readFilesCommand zcat '
		'--outSAMtype BAM Unsorted && rm {params.rmbam} && cd {params.wd}'

rule star_filter:
	input:
		'outs/STAR/{sample}_pass1/SJ.out.tab'
	params:
		wd = config['wd'],
		ID = '{sample}'
	output:
		'outs/STAR/{sample}_pass1/{sample}_Pass1SJ.filtered.tab'
	shell:
		'''awk "{{ if (\$7 >= 3) print \$0 }}" {input[0]} > {input[0]}.filtered && '''
		'cd outs/STAR/{params.ID}_pass1 && '
		'mv SJ.out.tab.filtered {params.ID}_Pass1SJ.filtered.tab && cd {params.wd}'

rule pass2:
	input:
		refdir = expand('{wd}/outs/STAR_index/{build}', wd = config['wd'], build = config['ref']['build']),
		SJfiles = 'outs/STAR/{sample}_pass1/{sample}_Pass1SJ.filtered.tab'
	params:
		wd = config['wd'],
		outdir = config['wd'] + '/outs/STAR/{sample}_pass2',
		ID = config['fastqs'] + '/' + '{sample}',
		ID_2 = '{sample}',
		outdir_2 = config['wd'] + '/outs/STAR/bams'
	output:
		'outs/STAR/bams/{sample}.Aligned.sortedByCoord.out.bam',
	threads: 4
	shell:
		'rm -rf {params.outdir} && '
		'mkdir {params.outdir} && '
		'cd {params.outdir} && '
		'STAR --runThreadN {threads} '
		'--genomeDir {input.refdir} '
		'--readFilesIn {params.ID}_1.fq.gz {params.ID}_2.fq.gz '
		'--readFilesCommand zcat '
		'--outSAMtype BAM SortedByCoordinate '
		'--sjdbFileChrStartEnd {params.wd}/{input.SJfiles} '
		'--outSAMattrRGline ID:{params.ID_2} && '
		'mv Aligned.sortedByCoord.out.bam ../bams/{params.ID_2}.Aligned.sortedByCoord.out.bam && '
		'cd {params.wd}'

rule raw_counts:
	input:
		gtf = config['ref']['gtf'],
		bams = expand('outs/STAR/bams/{sample}.Aligned.sortedByCoord.out.bam', sample = samples)
	output:
		'outs/counts/Pipeline.Counts.tsv'
	threads: 
		2
	conda: 
		'workflow/envs/raw_counts.yaml'
	log: 
		'logs/raw_counts.log'
	script:
		'workflow/scripts/create_counts.R'
