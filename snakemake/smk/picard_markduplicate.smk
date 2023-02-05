rule picard_markduplicate:
	input:
		"picard_reads/{sample}.clean.bam"
	output:
		"picard_reads/{sample}.mkdup.bam",
		"picard_reads/{sample}.dup.metrics.txt"
	shell:
		"""
		picard -Xmx64g MarkDuplicates \
		--INPUT {input} \
		--OUTPUT {output[0]} \
		--METRICS_FILE {output[1]} \
		--VALIDATION_STRINGENCY LENIENT \
		--CREATE_INDEX true
		"""
