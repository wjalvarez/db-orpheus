rule replace_rg:
	input:
		'outs/STAR/bams/{sample}.Aligned.sortedByCoord.out.bam'
	output:
		temp("outs/STAR/bams/{sample}.Aligned.sortedByCoord.out.rgAligned.bam")
	log:
		"logs/picard/replace_rg/{sample}.log"
	params:
		"RGID={sample} RGLB={sample} RGPL={sample} RGPU={sample} RGSM={sample} "
		"VALIDATION_STRINGENCY=LENIENT"
	message:
		"Including read group tag in BAM for recalibrating base quality score."
	wrapper:
		"0.57.0/bio/picard/addorreplacereadgroups"

rule mark_duplicates:
	input:
		"outs/STAR/bams/{sample}.Aligned.sortedByCoord.out.rgAligned.bam"
	output:
		bam = temp("outs/STAR/bams/{sample}.Aligned.sortedByCoord.out.markedAligned.bam"),
		metrics = "outs/STAR/bams/{sample}.metrics.txt"
	log:
		"logs/picard/dedup/{sample}.log"
	params:
		"REMOVE_DUPLICATES=true"
	wrapper:
		"0.57.0/bio/picard/markduplicates"

rule split_n_cigar_reads:
	input:
		bam = "outs/STAR/bams/{sample}.Aligned.sortedByCoord.out.markedAligned.bam",
		ref = config['ref']['fa']
	output:
		temp("outs/split/{sample}.bam")
	log:
		"logs/gatk/splitNCIGARreads/{sample}.log"
	params:
		extra = "-rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS",
		java_opts = ""
	wrapper:
		"0.57.0/bio/gatk/splitncigarreads"

rule gatk_bqsr:
	input:
		bam = "outs/split/{sample}.bam",
		ref = config['ref']['fa'],
		known = "/data/exploratory/Users/jeff.alvarez/pipeline_ins/00-All.vcf.gz"
	output:
		bam = temp("outs/recal/{sample}.bam")
	log:
		"logs/gatk/bqsr/{sample}.log"
	params:
		extra = "",
		java_opts = ""
	wrapper:
		"0.57.0/bio/gatk/baserecalibrator"

rule haplotype_caller:
	input:
		bam = "outs/recal/{sample}.bam",
		ref = config['ref']['fa']
	output:
		gvcf = temp("outs/calls/{sample}.g.vcf.gz")
	log:
		"logs/gatk/haplotypecaller/{sample}.log"
	threads:
		4
	params:
		extra = "-L 19 -dontUseSoftClippedBases -stand_call_conf 20.0",
		java_opts = ""
	wrapper:
		"0.57.0/bio/gatk/haplotypecaller"

rule combine_gvcfs:
	input:
		gvcfs = expand("outs/calls/{sample}.g.vcf.gz", sample = samples),
		ref = config['ref']['fa']
	output:
		gvcf = temp("outs/calls/all.g.vcf.gz")
	log:
		"logs/gatk/combinegvcfs.log"
	params:
		extra = "",
		java_opts = ""
	wrapper:
		"0.58.0/bio/gatk/combinegvcfs"

rule genotype_gvcfs:
	input:
		gvcf = "outs/calls/all.g.vcf.gz",
		ref = config['ref']['fa']
	output:
		vcf = "outs/calls/all.vcf.gz"
	log:
		"logs/gatk/genotypegvcfs.log"
	params:
		extra = "",
		java_opts = "",
	wrapper:
		"0.58.0/bio/gatk/genotypegvcfs"

rule gatk_filter:
	input:
		vcf = "outs/calls/all.vcf.gz",
		ref = config["ref"]["fa"],
	output:
		vcf = "outs/calls/all.filtered.vcf.gz"
	log:
		"logs/gatk/filter/snvs.log"
	params:
		filters = {"FS": "FS > 30.0", "QD": "QD < 2.0"},
		extra = "-window 35 -cluster 3",
		java_opts = "",
	wrapper:
		"0.59.1/bio/gatk/variantfiltration"
