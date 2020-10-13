rule replace_rg:
	input:
		config["bam"]
	output:
		temp("/dbfs/db-orpheus/tmp/{ID}/replace_rg/{sample}.bam")
	benchmark:
		"/dbfs/db-orpheus/benchmarks/{ID}/call/00_replace_rg/{sample}.txt"
	log:
		"/dbfs/db-orpheus/logs/{ID}/00_replace_rg/{sample}.log"
#		"logs/{ID}/call/00_replace_rg/{sample}.log"
	params:
		"RGID={sample} RGLB={sample} RGPL={sample} RGPU={sample} RGSM={sample} "
		"VALIDATION_STRINGENCY=SILENT"
	wrapper:
		"0.57.0/bio/picard/addorreplacereadgroups"

rule mark_duplicates:
	input:
		"/dbfs/db-orpheus/tmp/{ID}/replace_rg/{sample}.bam"
	output:
		bam = temp("/dbfs/db-orpheus/tmp/{ID}/mark_duplicates/{sample}.bam"),
		metrics = "/dbfs/db-orpheus/tmp/{ID}/mark_duplicates/{sample}.metrics.txt"
	benchmark:
		"/dbfs/db-orpheus/benchmarks/{ID}/call/01_mark_duplicates/{sample}.txt"
	log:
		"/dbfs/db-orpheus/logs/{ID}/01_mark_duplicates/{sample}.log"
	params:
		mem = "-Xmx8g"
	wrapper:
		"0.64.0/bio/picard/markduplicates"

rule split_n_cigar_reads:
	input:
		bam = "/dbfs/db-orpheus/tmp/{ID}/mark_duplicates/{sample}.bam",
		ref = config['ref']['fa']
	output:
		temp("/dbfs/db-orpheus/tmp/{ID}/split/{sample}.bam")
	benchmark:
		"/dbfs/db-orpheus/benchmarks/{ID}/call/02_split_n_cigar_reads/{sample}.txt"
	log:
		"/dbfs/db-orpheus/logs/{ID}/02_split_n_cigar_reads/{sample}.log"
	params:
		extra = "--tmp-dir /dbfs/db-orpheus/tmp",
		java_opts = "-Xmx8g"
	wrapper:
		"0.64.0/bio/gatk/splitncigarreads"

rule gatk_baserecalibrator:
	input:
		bam = "/dbfs/db-orpheus/tmp/{ID}/split/{sample}.bam",
		ref = config['ref']['fa'],
		dict = config["ref"]["dict"],
		known = config["ref"]["known_sites"]
	output:
		recal_table = temp("/dbfs/db-orpheus/tmp/{ID}/recal/{sample}.grp")
	benchmark:
		"/dbfs/db-orpheus/benchmarks/{ID}/call/03_gatk_bqsr/{sample}.txt"
	log:
		"/dbfs/db-orpheus/logs/{ID}/03_gatk_baserecalibrator/{sample}.log"
	params:
		extra = "-DF NotDuplicateReadFilter",
		java_opts = "-Xmx8g"
	wrapper:
		"0.64.0/bio/gatk/baserecalibrator"

rule gatk_applybqsr:
	input:
		bam = "/dbfs/db-orpheus/tmp/{ID}/split/{sample}.bam",
		ref = config['ref']['fa'],
		dict = config["ref"]["dict"],
		recal_table = "/dbfs/db-orpheus/tmp/{ID}/recal/{sample}.grp"
	output:
		bam = temp("/dbfs/db-orpheus/tmp/{ID}/recal/{sample}.bam")
	log:
		"/dbfs/db-orpheus/logs/{ID}/04_gatk_bqsr/{sample}.log"
	params:
		extra = "",
		java_opts = "-Xmx8g"
	wrapper:
		"0.64.0/bio/gatk/applybqsr"

rule haplotype_caller:
	input:
		bam = "/dbfs/db-orpheus/tmp/{ID}/recal/{sample}.bam",
		ref = config['ref']['fa']
	output:
		gvcf = temp("/dbfs/db-orpheus/tmp/{ID}/gvcfs/{sample}.g.vcf.gz")
	benchmark:
		"/dbfs/db-orpheus/benchmarks/{ID}/call/05_haplotype_caller/{sample}.txt"
	log:
		"/dbfs/db-orpheus/logs/{ID}/05_haplotype_caller/{sample}.log"
	threads:
		4
	params:
		extra = "--dont-use-soft-clipped-bases true -DF NotDuplicateReadFilter "
			"--minimum-mapping-quality 0 --base-quality-score-threshold 10 -mbq 13 "
			"-L /dbfs/references/Alu.RepeatMasker.hg19.ID.bed",
		java_opts = "-Xmx8g"
	wrapper:
		"0.64.0/bio/gatk/haplotypecaller"

rule genotype_gvcfs:
	input:
		gvcf = "/dbfs/db-orpheus/tmp/{ID}/gvcfs/{sample}.g.vcf.gz",
		ref = config['ref']['fa']
	output:
		vcf = temp("/dbfs/db-orpheus/tmp/{ID}/unfiltered/{sample}.unfiltered.vcf.gz")
	benchmark:
		"/dbfs/db-orpheus/benchmarks/{ID}/call/06_genotype_gvcfs.{sample}.txt"
	log:
		"/dbfs/db-orpheus/logs/{ID}/06_genotype_gvcfs/{sample}.log"
	params:
		extra = "-stand-call-conf 0.0",
		java_opts = "",
	wrapper:
		"0.64.0/bio/gatk/genotypegvcfs"

rule gatk_filter:
	input:
		vcf = "/dbfs/db-orpheus/tmp/{ID}/unfiltered/{sample}.unfiltered.vcf.gz",
		ref = config["ref"]["fa"]
	output:
		vcf = temp("/dbfs/db-orpheus/tmp/{ID}/filtered/{sample}.vcf.gz")
	benchmark:
		"/dbfs/db-orpheus/benchmarks/{ID}/call/07_gatk_filter.{sample}.txt"
	log:
		"/dbfs/db-orpheus/logs/{ID}/07_gatk_filter/{sample}.log"
	params:
		filters = {"FS": "FS > 30.0", "QD": "QD < 2.0", "DP": "DP < 20"},
		extra = "-window 35 -cluster 3",
		java_opts = "",
	wrapper:
		"0.64.0/bio/gatk/variantfiltration"

#rule snpeff_download:
#	output:
#		directory("outs/{ID}/snpeff_ref/{reference}")
#	log:
#		"/dbfs/db-orpheus/logs/{ID}/08_snpeff_annotate/{reference}.log"
#	params:
#		reference = "GRCh37.75"
#	wrapper:
#		"0.66.0/bio/snpeff/download"

rule snpeff:
	input:
		calls = temp("/dbfs/db-orpheus/tmp/{ID}/filtered/{sample}.vcf.gz")
#		db = "outs/{ID}/snpeff_ref/GRCh37.75"
	output:
		calls = temp("/dbfs/db-orpheus/tmp/{ID}/annotated/{sample}.vcf")
	log:
		"/dbfs/db-orpheus/logs/{ID}/08_snpeff_annotate/{sample}.log"
#		"logs/{ID}/08_snpeff_annotate/{sample}.log"
	params:
		extra = "-Xmx4g -no-downstream -no-intergenic -no-intron -no-upstream",
		reference = "GRCh37.75"
	conda:
		"../envs/snpeff.yaml"
	shell:
		"snpEff {params.extra} {params.reference} "
		"{input.calls} > {output.calls}"

rule bcftools_annotate:
	input:
		calls = "/dbfs/db-orpheus/tmp/{ID}/annotated/{sample}.vcf",
		bed = "/dbfs/references/Alu.RepeatMasker.hg19.ID.bed",
		header = "/dbfs/references/Alu.RepeatMasker.hg19.ID.txt"
#		bed = "/data/exploratory/Users/jeff.alvarez/pipeline_ins/Alu.RepeatMasker.hg19.ID.bed",
#		header = "/data/exploratory/Users/jeff.alvarez/pipeline_ins/Alu.RepeatMasker.hg19.ID.txt"
	output:
		vcf = "/dbfs/db-orpheus/{ID}/{sample}.vcf.gz"
#		vcf = "outs/{ID}/final/{sample}.vcf.gz"
	log:
		"/dbfs/db-orpheus/logs/{ID}/09_bcftools_annotate/{sample}.log"
	params:
		columns = "CHROM,FROM,TO,ALU_NAME,ALU_ID,STRAND"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"bcftools annotate -a {input.bed} -h {input.header} -c {params.columns} {input.calls} | "
		"bgzip -c > {output.vcf}"
