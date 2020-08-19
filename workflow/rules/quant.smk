rule rsem_prepare_reference:
	input:
		fasta = config["ref"]["fa"],
		gtf = config["ref"]["gtf"],
	threads:
		12
	params:
		build = config["ref"]["build"],
		analysis = config["ID"]
	output:
		directory("outs/{}/rsem/{}".format(config["ID"], config["ref"]["build"])),
	benchmark:
		"benchmarks/quant/00_rsem_index.txt"
	log:
		"logs/rsem_index_{}.log".format(config["ref"]["build"])
	conda:
		"workflow/envs/quant.yaml"
	shell:
		"rsem-prepare-reference --gtf {input.gtf} -p {threads} {input.fasta} {params.build}"

rule rsem_calculate_expression:
	input:
		bam = "outs/star/{sample}/Aligned.sortedByCoord.out.markedAligned.bam"
	threads:
		12
	output:
		directory("outs/RSEM/{sample}/{sample}")
