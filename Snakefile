import glob, os

configfile: "config.yaml"

samples, = glob_wildcards(config['fastqs'] + '/' + '{sample}_1.fq.gz')

rule all:
	input:
#		expand('outs/STAR_index/{ref}', ref = config['ref']),
#		expand('outs/STAR/{sample}_pass1/SJ.out.tab', sample = samples),
#		expand('outs/STAR/{sample}_pass1/{sample}_Pass1SJ.filtered.tab', sample = samples),
#		'outs/STAR/bams',
#		expand('outs/STAR/bams/{sample}.Aligned.sortedByCoord.out.bam', sample = samples),
		'outs/counts/Pipeline.Counts.tsv'

rule star_index:
	input:
		fa = expand('{fa}/{ref}', fa = config['fa'], ref = config['ref']),
		gtf = expand('{gtf}/{ref}', gtf = config['gtf'], ref = config['ref'])
	output:
		directory('outs/STAR_index/{ref}')
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
		refdir = expand('{currdir}/outs/STAR_index/{ref}', currdir = config['wd'], ref = config['ref'])
	params:
		outdir = 'outs/STAR/{sample}_pass1',
		rmbam = config['wd'] + '/' + 'outs/STAR/{sample}_pass1/Aligned.out.bam',
		ID = config['wd'] + '/' + config['fastqs'] + '/' + '{sample}',
		wd = config['wd']
	output:
		'outs/STAR/{sample}_pass1/SJ.out.tab'
	threads: 2
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
		refdir = expand('{currdir}/outs/STAR_index/{ref}', currdir = config['wd'], ref = config['ref']),
		SJfiles = 'outs/STAR/{sample}_pass1/{sample}_Pass1SJ.filtered.tab'
	params:
		wd = config['wd'],
		outdir = config['wd'] + '/outs/STAR/{sample}_pass2',
		ID = config['wd'] + '/' + config['fastqs'] + '/' + '{sample}',
		ID_2 = '{sample}',
		outdir_2 = config['wd'] + '/outs/STAR/bams'
	output:
		'outs/STAR/bams/{sample}.Aligned.sortedByCoord.out.bam',
	threads: 2
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
		gtf = config['wd'] + '/' + config['gtf'],
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
