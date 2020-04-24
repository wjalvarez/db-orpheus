rule fastqc:
	input:
#		expand('{fastqs}/{sample}_{pair}.fq.gz', fastqs = config['fastqs'], sample = samples, pair = pairs)
		config['fastqs'] + '/{sample}_{pair}.fq.gz' 
	output:
		html="outs/qc/{sample}_{pair}.html",
		zip="outs/qc/{sample}_{pair}_fastqc.zip"
	wrapper:
		"0.31.1/bio/fastqc"

rule multiqc:
	input:
		expand("outs/qc/{sample}_{pair}_fastqc.zip", sample = samples, pair = pairs)
	output:
		"outs/qc/multiqc_report.html"
	wrapper:
		"0.51.3/bio/multiqc"
