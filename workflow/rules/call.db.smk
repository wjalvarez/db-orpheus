rule replace_rg:
	input:
		config["bam"]
	output:
		temp("/dbfs/tmp-db-orpheus/{ID}.{sample}.replace_rg.bam")
	benchmark:
		"/dbfs/db-orpheus/benchmarks/{ID}/call/{sample}.00_replace_rg.txt"
	log:
		"/dbfs/db-orpheus/logs/{ID}/call/{sample}.00_replace_rg.log"
	params:
		"RGID={sample} RGLB={sample} RGPL={sample} RGPU={sample} RGSM={sample} "
		"VALIDATION_STRINGENCY=SILENT"
	wrapper:
		"0.57.0/bio/picard/addorreplacereadgroups"

rule split_n_cigar_reads:
	input:
		bam = "/dbfs/tmp-db-orpheus/{ID}.{sample}.replace_rg.bam",
		ref = config['ref']['fa']
	output:
		temp("/dbfs/tmp-db-orpheus/{ID}.{sample}.split.bam")
	benchmark:
		"/dbfs/db-orpheus/benchmarks/{ID}/call/{sample}.01_split_n_cigar_reads.txt"
	log:
		"/dbfs/db-orpheus/logs/{ID}/call/{sample}.01_split.bam"
	params:
		extra = "",
		java_opts = "-Xmx8g"
	wrapper:
		"0.64.0/bio/gatk/splitncigarreads"

rule gatk_baserecalibrator:
	input:
		bam = "/dbfs/tmp-db-orpheus/{ID}.{sample}.split.bam",
		ref = config['ref']['fa'],
		dict = config["ref"]["dict"],
		known = config["ref"]["known_sites"]
	output:
		recal_table = temp("/dbfs/tmp-db-orpheus/{ID}.{sample}.recal.grp")
	benchmark:
		"/dbfs/db-orpheus/benchmarks/{ID}/call/{sample}.02_gatk_bqsr.txt"
	log:
		"/dbfs/db-orpheus/logs/{ID}/call/{sample}.02_gatk_bqsr.txt"
	params:
		extra = "-DF NotDuplicateReadFilter",
		java_opts = "-Xmx8g"
	wrapper:
		"0.64.0/bio/gatk/baserecalibrator"

rule gatk_applybqsr:
	input:
		bam = "/dbfs/tmp-db-orpheus/{ID}.{sample}.split.bam",
		ref = config['ref']['fa'],
		dict = config["ref"]["dict"],
		recal_table = "/dbfs/tmp-db-orpheus/{ID}.{sample}.recal.grp"
	output:
		bam = temp("/dbfs/tmp-db-orpheus/{ID}.{sample}.recal.bam")
	log:
		"/dbfs/db-orpheus/logs/{ID}/call/{sample}.03_apply_bqsr.log"
	benchmark:
		"/dbfs/db-orpheus/benchmarks/{ID}/call/{sample}.03_apply_bqsr.log"
	params:
		extra = "",
		java_opts = "-Xmx8g"
	wrapper:
		"0.64.0/bio/gatk/applybqsr"

rule haplotype_caller:
	input:
		bam = "/dbfs/tmp-db-orpheus/{ID}.{sample}.recal.bam",
		ref = config['ref']['fa']
	output:
		gvcf = temp("/dbfs/tmp-HC/{ID}.{sample}.{chr_chunks}.g.vcf.gz")
	benchmark:
		"/dbfs/db-orpheus/benchmarks/{ID}/call/{sample}.{chr_chunks}.04_haplotype_caller.txt"
	log:
		"/dbfs/db-orpheus/logs/{ID}/call/{sample}.{chr_chunks}.04_haplotype_caller.txt"
	params:
		chr_intervals = lambda wildcards: chr_dict[wildcards.chr_chunks],
		extra = "--dont-use-soft-clipped-bases true -DF NotDuplicateReadFilter "
			"--minimum-mapping-quality 0 --base-quality-score-threshold 10 -mbq 13",
		java_opts = "-Xmx8g"
	conda:
		"../envs/gatk4.yaml"
	shell:
		"gatk --java-options {params.java_opts} HaplotypeCaller {params.extra} "
		"{params.chr_intervals} -R {input.ref} -I {input.bam} -ERC GVCF -O {output.gvcf}"

rule combine_gvcfs:
	input:
		gvcfs = expand("/dbfs/tmp-HC/{ID}.{sample}.{chr_chunks}.g.vcf.gz", 
			ID = ID, sample = sample, chr_chunks = chr_chunks),
		ref = config['ref']['fa']
	output:
		gvcf = temp("/dbfs/tmp-db-orpheus/{ID}.{sample}.g.vcf.gz")
	benchmark:
		"/dbfs/db-orpheus/benchmarks/{ID}/call/{sample}.05_combine_gvcfs.txt"
	log:
		"/dbfs/db-orpheus/logs/{ID}/call/{sample}.05_combine_gvcfs.txt"
	params:
		extra = "",
		java_opts = ""
	wrapper:
		"0.58.0/bio/gatk/combinegvcfs"


rule genotype_gvcfs:
	input:
		gvcf = "/dbfs/tmp-db-orpheus/{ID}.{sample}.g.vcf.gz",
		ref = config['ref']['fa']
	output:
		vcf = temp("/dbfs/tmp-db-orpheus/{ID}.{sample}.unfiltered.vcf.gz")
	benchmark:
		"/dbfs/db-orpheus/benchmarks/{ID}/call/{sample}.06_genotype_gvcfs.txt"
	log:
		"/dbfs/db-orpheus/logs/{ID}/call/{sample}.06_genotype_gvcfs.log"
	params:
		extra = "-stand-call-conf 0.0",
		java_opts = "",
	wrapper:
		"0.64.0/bio/gatk/genotypegvcfs"

rule gatk_filter:
	input:
		vcf = "/dbfs/tmp-db-orpheus/{ID}.{sample}.unfiltered.vcf.gz",
		ref = config["ref"]["fa"]
	output:
		vcf = temp("/dbfs/tmp-db-orpheus/outs/{ID}/{sample}.filtered.vcf.gz")
	benchmark:
		"/dbfs/db-orpheus/benchmarks/{ID}/call/{sample}.07_gatk_filter.txt"
	log:
		"/dbfs/db-orpheus/logs/{ID}/call/{sample}.07_gatk_filter.log"
	params:
		filters = {"FS": "FS > 30.0", "QD": "QD < 2.0", "DP": "DP < 20"},
		extra = "-window 35 -cluster 3",
		java_opts = "",
	wrapper:
		"0.64.0/bio/gatk/variantfiltration"

rule snpeff:
	input:
		calls = "/dbfs/tmp-db-orpheus/outs/{ID}/{sample}.filtered.vcf.gz"
	output:
		calls = temp("/dbfs/tmp-db-orpheus/{ID}.{sample}.annotated.vcf")
	benchmark:
		"/dbfs/db-orpheus/benchmarks/{ID}/call/{sample}.08_snpeff.txt"
	log:
		"/dbfs/db-orpheus/logs/{ID}/call/{sample}.08_snpeff.log"
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
		calls = "/dbfs/tmp-db-orpheus/{ID}.{sample}.annotated.vcf",
		bed = "/dbfs/references/Alu.RepeatMasker.hg19.ID.bed",
		header = "/dbfs/references/Alu.RepeatMasker.hg19.ID.txt"
	output:
		vcf = "/dbfs/db-orpheus/{ID}/{sample}.vcf.gz"
	benchmark:
		"/dbfs/db-orpheus/benchmarks/{ID}/call/{sample}.09_bcftools_annotate.txt"
	log:
		"/dbfs/db-orpheus/logs/{ID}/09_bcftools_annotate/{sample}.log"
	params:
		columns = "CHROM,FROM,TO,ALU_NAME,ALU_ID,STRAND"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"""bcftools annotate -a {input.bed} -h {input.header} -c {params.columns} {input.calls} | """
		"""bcftools filter -i "REF == 'A' & ALT == 'G' | REF == 'T' & ALT == 'C'" | """
		"""bcftools filter -e "FORMAT/AD[:1] < 5" | bgzip -c > {output.vcf}"""
