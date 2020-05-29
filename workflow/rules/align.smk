rule star_index_new:
	input:
		fasta = config["ref"]["fa"],
		gtf = config["ref"]["gtf"],
	message:
		"Testing STAR index"
	threads:
		12
	params:
		extra = "",
		build = config["ref"]["build"],
		analysis = config["Title"],
	output:
		directory("outs/{}/{}".format(config["Title"], config["ref"]["build"])),
	log:
		"logs/star_index_{input.build}}.log"
	wrapper:
		"0.59.1/bio/star/index"

rule star_pe_multi:
	input:
		fq1 = ["reads/{sample}_R1.fq.gz"],
		fq2 = ["reads/{sample}_R2.fq.gz"]
	output:
		"outs/star/pe/{sample}.Aligned.sortedByCoord.out.bam"
	log:
		"logs/star/pe/{sample}.log"
	params:
		index = "index",
		extra = "--twopassMode Basic --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate"
	threads:
		4
	wrapper:
		"0.59.1/bio/star/align"
