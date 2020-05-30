rule star_index_new:
	input:
		fasta = config["ref"]["fa"],
		gtf = config["ref"]["gtf"],
	threads:
		12
	params:
		extra = "",
		build = config["ref"]["build"],
		analysis = config["ID"]
	output:
		directory("outs/{}/{}".format(config["ID"], config["ref"]["build"])),
	benchmark:
		"benchmarks/align/00_star_index.txt"
	log:
		"logs/star_index_{}.log".format(config["ref"]["build"])
	wrapper:
		"0.59.1/bio/star/index"

rule star_pe_multi:
	input:
		directory("outs/{}/{}".format(config["ID"], config["ref"]["build"])),
		fq1 = [config["fastqs"] + "{sample}_1.fq.gz"],
		fq2 = [config["fastqs"] + "{sample}_2.fq.gz"]
	output:
		"outs/star/{sample}/Aligned.sortedByCoord.out.bam"
	benchmark:
		"benchmarks/align/01_star_align.{sample}.txt"
	log:
		"logs/star/{sample}.log"
	params:
		index = "outs/{}/{}".format(config["ID"], config["ref"]["build"]),
		extra = "--twopassMode Basic --outSAMtype BAM SortedByCoordinate"
	threads:
		12
	wrapper:
		"0.59.1/bio/star/align"
