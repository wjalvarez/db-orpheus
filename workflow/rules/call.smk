rule replace_rg:
	input:
#		'outs/star/{sample}/Aligned.sortedByCoord.out.bam',
		config["bam"]
	output:
		temp("outs/star/{sample}/Aligned.sortedByCoord.out.rgAligned.bam")
#		temp("/dbfs/db-orpheus/tmp/00_replace_rg/{sample}.rgAligned.bam")
	benchmark:
		"benchmarks/call/00_replace_rg.{sample}.txt"
	log:
#		"logs/picard/replace_rg/{sample}.log"
		"/dbfs/db-orpheus/logs/00_replace_rg/{sample}.log"
	params:
		"RGID={sample} RGLB={sample} RGPL={sample} RGPU={sample} RGSM={sample} "
		"VALIDATION_STRINGENCY=SILENT"
	wrapper:
		"0.64.0/bio/picard/addorreplacereadgroups"

rule mark_duplicates:
	input:
		"outs/star/{sample}/Aligned.sortedByCoord.out.rgAligned.bam"
		#"/dbfs/db-orpheus/tmp/00_replace_rg/{sample}.rgAligned.bam"
	output:
		bam = temp("outs/star/{sample}/Aligned.sortedByCoord.out.markedAligned.bam"),
		#bam = temp("/dbfs/db-orpheus/tmp/01_mark_duplicates/{sample}.markDup.bam"),
		metrics = "outs/star/{sample}/metrics.txt"
	benchmark:
		"benchmarks/call/01_mark_duplicates.{sample}.txt"
	log:
#		"logs/picard/dedup/{sample}.log"
		"/dbfs/db-orpheus/logs/01_mark_duplicates/{sample}.log"
	params:
		""
	wrapper:
		"0.64.0/bio/picard/markduplicates"

rule split_n_cigar_reads:
	input:
		bam = "outs/star/{sample}/Aligned.sortedByCoord.out.markedAligned.bam",
		#bam = "/dbfs/db-orpheus/tmp/01_mark_duplicates/{sample}.markDup.bam",
		ref = config['ref']['fa']
	output:
		temp("outs/split/{sample}.bam")
		#temp("/dbfs/db-orpheus/tmp/02_split_n_cigar_reads/{sample}.split.bam")
	benchmark:
		"benchmarks/call/02_split_n_cigar_reads.{sample}.txt"
	log:
#		"logs/gatk/splitNCIGARreads/{sample}.log"
		"/dbfs/db-orpheus/logs/02_split_n_cigar_reads/{sample}.log"
	params:
		extra = "",
		java_opts = ""
	wrapper:
		"0.64.0/bio/gatk/splitncigarreads"

rule gatk_baserecalibrator:
	input:
		bam = "outs/split/{sample}.bam",
#		bam = "/dbfs/db-orpheus/tmp/02_split_n_cigar_reads/{sample}.split.bam",
		ref = config['ref']['fa'],
		dict = config["ref"]["dict"],
		known = config["ref"]["known_sites"]
	output:
		recal_table = temp("outs/recal/{sample}.grp")
		#bam = temp("/dbfs/db-orpheus/tmp/03_gatk_bqsr/{sample}.bqsr.bam")
	benchmark:
		"benchmarks/call/03_gatk_bqsr.{sample}.txt"
	log:
#		"logs/gatk/bqsr/baserecalibrator.{sample}.log"
		"/dbfs/db-orpheus/logs/03_gatk_baserecalibrator/{sample}.log"
	params:
		extra = "-DF NotDuplicateReadFilter",
		java_opts = ""
	wrapper:
		"0.64.0/bio/gatk/baserecalibrator"

rule gatk_applybqsr:
	input:
		bam = "outs/split/{sample}.bam",
#		bam = "/dbfs/db-orpheus/tmp/02_split_n_cigar_reads/{sample}.split.bam",
		ref = config['ref']['fa'],
		dict = config["ref"]["dict"],
		recal_table = "outs/recal/{sample}.grp"
	output:
		bam = temp("outs/recal/{sample}.bam")
	log:
#		"logs/gatk/bqsr/gatk_applybqsr.{sample}.log"
		"/dbfs/db-orpheus/logs/04_gatk_bqsr/{sample}.log"
	params:
		extra = "",
		java_opts = ""
	wrapper:
		"0.64.0/bio/gatk/applybqsr"

#rule gatk_bqsr_old:
#	input:
#		bam = "outs/split/{sample}.bam",
#		bam = "/dbfs/db-orpheus/tmp/02_split_n_cigar_reads/{sample}.split.bam",
#		ref = config['ref']['fa'],
#		known = config["ref"]["known_sites"]
#	output:
#		bam = temp("outs/recal/{sample}.bam")
		#bam = temp("/dbfs/db-orpheus/tmp/03_gatk_bqsr/{sample}.bqsr.bam")
#	benchmark:
#		"benchmarks/call/03_gatk_bqsr.{sample}.txt"
#	log:
#		"logs/gatk/bqsr/{sample}.log"
#		"/dbfs/db-orpheus/logs/03_gatk_bqsr/{sample}.log"
#	params:
#		extra = "-DF NotDuplicateReadFilter",
#		java_opts = ""
#	wrapper:
#		"0.57.0/bio/gatk/baserecalibrator"

rule haplotype_caller:
	input:
		bam = "outs/recal/{sample}.bam",
		#bam = "/dbfs/db-orpheus/tmp/03_gatk_bqsr/{sample}.bqsr.bam",
		ref = config['ref']['fa']
	output:
		gvcf = temp("outs/gvcfs/{sample}.g.vcf.gz")
#		gvcf = temp("/dbfs/db-orpheus/tmp/04_haplotype_caller/{sample}.g.vcf.gz")
	benchmark:
		"benchmarks/call/04_haplotype_caller.{sample}.txt"
	log:
		#"logs/gatk/haplotypecaller/{sample}.log"
		"/dbfs/db-orpheus/logs/05_haplotype_caller/{sample}.log"
	threads:
		4
	params:
		extra = "--dont-use-soft-clipped-bases true -stand-call-conf 10.0 "
			"-DF NotDuplicateReadFilter --base-quality-score-threshold 10",
		java_opts = ""
	wrapper:
		"0.64.0/bio/gatk/haplotypecaller"

rule genotype_gvcfs:
	input:
		gvcf = "outs/gvcfs/{sample}.g.vcf.gz",
#		gvcf = "/dbfs/db-orpheus/tmp/04_haplotype_caller/{sample}.g.vcf.gz",
		ref = config['ref']['fa']
	output:
		vcf = temp("outs/unfiltered/{sample}.unfiltered.vcf.gz")
#		vcf = temp("/dbfs/db-orpheus/tmp/05_genotype_gvcfs/{sample}.unfiltered.vcf.gz")
	benchmark:
		"benchmarks/call/06_genotype_gvcfs.{sample}.txt"
	log:
#		"logs/gatk/genotype_gvcfs.{sample}.log"
		"/dbfs/db-orpheus/logs/06_genotype_gvcfs/{sample}.log"
	params:
		extra = "",
		java_opts = "",
	wrapper:
		"0.64.0/bio/gatk/genotypegvcfs"

rule gatk_filter:
	input:
		vcf = "outs/unfiltered/{sample}.unfiltered.vcf.gz",
#		vcf = "/dbfs/db-orpheus/tmp/05_genotype_gvcfs/{sample}.unfiltered.vcf.gz",
		ref = config["ref"]["fa"],
	output:
#		vcf = "outs/calls/{sample}.vcf.gz"
		vcf = "/dbfs/db-orpheus/{sample}.vcf.gz"
	benchmark:
		"benchmarks/call/07_gatk_filter.{sample}.txt"
	log:
	#	"logs/gatk/filter/snvs.{sample}.log"
		"/dbfs/db-orpheus/logs/07_gatk_filter/{sample}.log"
	params:
		filters = {"FS": "FS > 30.0", "QD": "QD < 2.0"},
		extra = "-window 35 -cluster 3",
		java_opts = "",
	wrapper:
		"0.64.0/bio/gatk/variantfiltration"
