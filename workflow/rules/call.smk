rule replace_rg:
	input:
		'outs/STAR/bams/{sample}.Aligned.sortedByCoord.out.bam'
	output:
		"outs/STAR/bams/{sample}.Aligned.sortedByCoord.out.rgAligned.bam"
	log:
		"logs/picard/replace_rg/{sample}.log"
	params:
		"RGID={sample} RGLB={sample} RGPL={sample} RGPU={sample} RGSM={sample} "
		"VALIDATION_STRINGENCY=LENIENT"
	message:
		"Including read group tag in BAM for recalibrating base quality score."
	wrapper:
		"0.57.0/bio/picard/addorreplacereadgroups"
