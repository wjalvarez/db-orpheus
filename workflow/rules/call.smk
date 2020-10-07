rule replace_rg:
	input:
		config["bam"]
	output:
		temp("outs/{ID}/star/{sample}/Aligned.sortedByCoord.out.rgAligned.bam")
	benchmark:
		"benchmarks/{ID}/call/00_replace_rg/{sample}.txt"
	log:
		"/dbfs/db-orpheus/logs/{ID}/call/00_replace_rg/{sample}.log"
	params:
		"RGID={sample} RGLB={sample} RGPL={sample} RGPU={sample} RGSM={sample} "
		"VALIDATION_STRINGENCY=SILENT"
	wrapper:
		"0.57.0/bio/picard/addorreplacereadgroups"

rule mark_duplicates:
	input:
		"outs/{ID}/star/{sample}/Aligned.sortedByCoord.out.rgAligned.bam"
	output:
		bam = temp("outs/{ID}/star/{sample}/Aligned.sortedByCoord.out.markedAligned.bam"),
		metrics = "outs/{ID}/star/{sample}/metrics.txt"
	benchmark:
		"benchmarks/{ID}/call/01_mark_duplicates/{sample}.txt"
	log:
		"/dbfs/db-orpheus/logs/{ID}/01_mark_duplicates/{sample}.log"
	params:
		mem = "-Xmx8g"
	wrapper:
		"0.64.0/bio/picard/markduplicates"

rule split_n_cigar_reads:
	input:
		bam = "outs/{ID}/star/{sample}/Aligned.sortedByCoord.out.markedAligned.bam",
		ref = config['ref']['fa']
	output:
		temp("outs/{ID}/split/{sample}.bam")
	benchmark:
		"benchmarks/{ID}/call/02_split_n_cigar_reads/{sample}.txt"
	log:
		"/dbfs/db-orpheus/logs/{ID}/02_split_n_cigar_reads/{sample}.log"
	params:
		extra = "--tmp-dir outs/{ID}/star/{sample}",
		java_opts = "-Xmx4g"
	wrapper:
		"0.64.0/bio/gatk/splitncigarreads"

rule gatk_baserecalibrator:
	input:
		bam = "outs/{ID}/split/{sample}.bam",
		ref = config['ref']['fa'],
		dict = config["ref"]["dict"],
		known = config["ref"]["known_sites"]
	output:
		recal_table = temp("outs/{ID}/recal/{sample}.grp")
	benchmark:
		"benchmarks/{ID}/call/03_gatk_bqsr/{sample}.txt"
	log:
		"/dbfs/db-orpheus/logs/{ID}/03_gatk_baserecalibrator/{sample}.log"
	params:
		extra = "-DF NotDuplicateReadFilter",
		java_opts = "-Xmx4g"
	wrapper:
		"0.64.0/bio/gatk/baserecalibrator"

rule gatk_applybqsr:
	input:
		bam = "outs/{ID}/split/{sample}.bam",
		ref = config['ref']['fa'],
		dict = config["ref"]["dict"],
		recal_table = "outs/{ID}/recal/{sample}.grp"
	output:
		bam = temp("outs/{ID}/recal/{sample}.bam")
	log:
		"/dbfs/db-orpheus/logs/{ID}/04_gatk_bqsr/{sample}.log"
	params:
		extra = "",
		java_opts = "-Xmx4g"
	wrapper:
		"0.64.0/bio/gatk/applybqsr"

rule haplotype_caller:
	input:
		bam = "outs/{ID}/recal/{sample}.bam",
		ref = config['ref']['fa']
	output:
		gvcf = temp("outs/{ID}/gvcfs/{sample}.g.vcf.gz")
	benchmark:
		"benchmarks/{ID}/call/05_haplotype_caller/{sample}.txt"
	log:
		"/dbfs/db-orpheus/logs/{ID}/05_haplotype_caller/{sample}.log"
	threads:
		4
	params:
		extra = "--dont-use-soft-clipped-bases true -DF NotDuplicateReadFilter "
			"--minimum-mapping-quality 0 --base-quality-score-threshold 10 -mbq 13",
		java_opts = ""
	wrapper:
		"0.64.0/bio/gatk/haplotypecaller"

rule genotype_gvcfs:
	input:
		gvcf = "outs/{ID}/gvcfs/{sample}.g.vcf.gz",
		ref = config['ref']['fa']
	output:
		vcf = temp("outs/{ID}/unfiltered/{sample}.unfiltered.vcf.gz")
	benchmark:
		"benchmarks/{ID}/call/06_genotype_gvcfs.{sample}.txt"
	log:
		"/dbfs/db-orpheus/logs/{ID}/06_genotype_gvcfs/{sample}.log"
	params:
		extra = "",
		java_opts = "",
	wrapper:
		"0.64.0/bio/gatk/genotypegvcfs"

rule gatk_filter:
	input:
		vcf = "outs/{ID}/unfiltered/{sample}.unfiltered.vcf.gz",
		ref = config["ref"]["fa"],
	output:
		vcf = "/dbfs/db-orpheus/{ID}/{sample}.vcf.gz"
	benchmark:
		"benchmarks/{ID}/call/07_gatk_filter.{sample}.txt"
	log:
		"/dbfs/db-orpheus/logs/{ID}/07_gatk_filter/{sample}.log"
	params:
		filters = {"FS": "FS > 30.0", "QD": "QD < 2.0", "DP": "DP < 20"},
		extra = "-window 35 -cluster 3",
		java_opts = "",
	wrapper:
		"0.64.0/bio/gatk/variantfiltration"
