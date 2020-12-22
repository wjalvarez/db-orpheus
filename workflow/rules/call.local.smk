rule replace_rg:
	input:
		config["bam"]
	output:
		temp("outs/{ID}/replace/{sample}.replace_rg.bam")
	benchmark:
		"benchmarks/{ID}/call/{sample}.00_replace_rg.txt"
	log:
		"logs/{ID}/call/{sample}.00_replace_rg.log"
	params:
		"RGID={sample} RGLB={sample} RGPL={sample} RGPU={sample} RGSM={sample} "
		"VALIDATION_STRINGENCY=SILENT"
	wrapper:
		"0.57.0/bio/picard/addorreplacereadgroups"

rule split_n_cigar_reads:
	input:
		bam = "outs/{ID}/replace/{sample}.replace_rg.bam",
		ref = config['ref']['fa']
	output:
		temp("outs/{ID}/split/{sample}.split.bam")
	benchmark:
		"benchmarks/{ID}/call/{sample}.01_split_n_cigar_reads.txt"
	log:
		"logs/{ID}/call/{sample}.01_split.bam"
	params:
		extra = "",
		java_opts = "-Xmx8g"
	wrapper:
		"0.64.0/bio/gatk/splitncigarreads"

rule gatk_baserecalibrator:
	input:
		bam = "outs/{ID}/split/{sample}.split.bam",
		ref = config['ref']['fa'],
		dict = config["ref"]["dict"],
		known = config["ref"]["known_sites"]
	output:
		recal_table = temp("outs/{ID}/recal/{sample}.recal.grp")
	benchmark:
		"benchmarks/{ID}/call/{sample}.02_gatk_bqsr.txt"
	log:
		"logs/{ID}/call/{sample}.02_gatk_bqsr.txt"
	params:
		extra = "-DF NotDuplicateReadFilter",
		java_opts = "-Xmx8g"
	wrapper:
		"0.64.0/bio/gatk/baserecalibrator"

rule gatk_applybqsr:
	input:
		bam = "outs/{ID}/split/{sample}.split.bam",
		ref = config['ref']['fa'],
		dict = config["ref"]["dict"],
		recal_table = "outs/{ID}/recal/{sample}.recal.grp"
	output:
		bam = temp("outs/{ID}/recal/{sample}.recal.bam")
	log:
		"logs/{ID}/call/{sample}.03_apply_bqsr.log"
	benchmark:
		"benchmarks/{ID}/call/{sample}.03_apply_bqsr.log"
	params:
		extra = "",
		java_opts = "-Xmx8g"
	wrapper:
		"0.64.0/bio/gatk/applybqsr"

def get_intervals(wildcards):
	return {'chr_chunks': wildcards.chr_chunks,
		'chr_intervals': chr_dict[wildcards.chr_chunks]}

rule haplotype_caller:
	input:
#		unpack(get_intervals),
		bam = "outs/{ID}/recal/{sample}.recal.bam",
		ref = config['ref']['fa'],
	#	chr_intervals = lambda wildcards: chr_dict[wildcards.chr_chunks]
	output:
		gvcf = temp("outs/{ID}/gvcfs/{sample}.{chr_chunks}.g.vcf.gz")
	benchmark:
		"benchmarks/{ID}/call/{sample}.{chr_chunks}.04_haplotype_caller.txt"
	log:
		"logs/{ID}/call/{sample}.{chr_chunks}.04_haplotype_caller.txt"
	params:
		chr_intervals = lambda wildcards: chr_dict[wildcards.chr_chunks],
		extra = "--dont-use-soft-clipped-bases true -DF NotDuplicateReadFilter "
			"--minimum-mapping-quality 0 --base-quality-score-threshold 10 -mbq 13",
#			"{params.chr_intervals}",
		java_opts = "-Xmx8g"
#	wrapper:
#		"0.64.0/bio/gatk/haplotypecaller"
	conda:
		"../envs/gatk4.yaml"
	shell:
		"gatk --java-options {params.java_opts} HaplotypeCaller {params.extra} "
		"{params.chr_intervals} -R {input.ref} -I {input.bam} -ERC GVCF -O {output.gvcf}"

rule combine_gvcfs:
	input:
		gvcfs = expand("outs/{ID}/gvcfs/{sample}.{chr_chunks}.g.vcf.gz", ID = ID, sample = sample,
				chr_chunks = chr_chunks),
		ref = config['ref']['fa']
	output:
		gvcf = temp("outs/{ID}/combine_gvcfs/{sample}.g.vcf.gz")
	benchmark:
		"benchmarks/{ID}/call/{sample}.05_combine_gvcfs.txt"
	log:
		"logs/{ID}/call/{sample}.05_combine_gvcfs.txt"
	params:
		extra = "",
		java_opts = ""
	wrapper:
		"0.58.0/bio/gatk/combinegvcfs"


rule genotype_gvcfs:
	input:
		gvcf = "outs/{ID}/combine_gvcfs/{sample}.g.vcf.gz",
		ref = config['ref']['fa']
	output:
		vcf = temp("outs/{ID}/{sample}.unfiltered.vcf.gz")
	benchmark:
		"benchmarks/{ID}/call/{sample}.06_genotype_gvcfs.txt"
	log:
		"logs/{ID}/call/{sample}.06_genotype_gvcfs.log"
	params:
		extra = "-stand-call-conf 0.0",
		java_opts = "",
	wrapper:
		"0.64.0/bio/gatk/genotypegvcfs"

rule gatk_filter:
	input:
		vcf = "outs/{ID}/{sample}.unfiltered.vcf.gz",
		ref = config["ref"]["fa"]
	output:
		vcf = temp("outs/{ID}/{sample}.filtered.vcf.gz")
	benchmark:
		"benchmarks/{ID}/call/{sample}.07_gatk_filter.txt"
	log:
		"logs/{ID}/call/{sample}.07_gatk_filter.log"
	params:
		filters = {"FS": "FS > 30.0", "QD": "QD < 2.0", "DP": "DP < 20"},
		extra = "-window 35 -cluster 3",
		java_opts = "",
	wrapper:
		"0.64.0/bio/gatk/variantfiltration"

rule snpeff:
	input:
		calls = "outs/{ID}/{sample}.filtered.vcf.gz"
	output:
		calls = temp("outs/{ID}/annotated/{sample}.vcf")
	benchmark:
		"benchmarks/{ID}/call/{sample}.08_snpeff.txt"
	log:
		"logs/{ID}/call/{sample}.08_snpeff.log"
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
		calls = "outs/{ID}/annotated/{sample}.vcf",
		bed = "/data/exploratory/Users/jeff.alvarez/pipeline_ins/Alu.RepeatMasker.hg19.ID.bed",
		header = "/data/exploratory/Users/jeff.alvarez/pipeline_ins/Alu.RepeatMasker.hg19.ID.txt"
	output:
		vcf = "outs/{ID}/final/{sample}.vcf.gz"
	benchmark:
		"benchmarks/{ID}/call/{sample}.09_bcftools_annotate.txt"
	log:
		"logs/{ID}/09_bcftools_annotate/{sample}.log"
	params:
		columns = "CHROM,FROM,TO,ALU_NAME,ALU_ID,STRAND"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"bcftools annotate -a {input.bed} -h {input.header} -c {params.columns} {input.calls} | "
		"bgzip -c > {output.vcf}"
